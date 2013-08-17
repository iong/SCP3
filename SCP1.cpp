//
//  SCP1.cpp
//  SCP3
//
//  Created by Ionut Georgescu on 7/15/13.
//
//

#include <iostream>

#include <mpi.h>

#include "DiskIO.h"
#include "Hessian.h"
#include "SCP1.h"
#include "sobol.hpp"

vec SCP1::disp_disp_corr_diag(const vec& omega, double kT)
{
    vec d(omega.n_rows);
    if (kT == 0.0) {
        d = 0.5 / omega;
    }
    else {
        d = 0.5 / (omega * tanh(0.5/kT * omega) );
    }

    d ( find(omega < (omega_cutoff)) ).fill(0.0);

    return d;
}

double SCP1::average(const vec& q, const vec& d, const mat& isqrtM_U, vec& UX, mat& UX_q)
{
    double Upot;
    int dim = q.n_rows;
    int nw = dim / 9;
    
    Upot = 0.0;
    UX.zeros(dim);
    UX_q.zeros(dim, dim);

    {
        h2o::PES *pot = getPES();
        double Upot_private;
        vec UX_private, UX_(dim), y(dim), dq(dim);
        mat UX_q_private;
        
        Upot_private = 0.0;
        UX_private.zeros(dim);
        UX_q_private.zeros(dim, dim);
        
        for (int i=(rank * NSobol) / nprocs;
                i < ((rank+1) * NSobol) / nprocs; i++) {
            long long current_skip = sobol_skip + i;
            sobol::std_normal(y.n_rows, &current_skip, y.memptr());
            //cout << omp_get_thread_num() <<" "<< current_skip << endl;
            
            y %= sqrt(d);
            
            dq = isqrtM_U * y;
            y = bohr * (q + dq);
            Upot_private += (*pot)(nw, y.memptr(), UX_.memptr()) / autokcalpmol;
            UX_ *= bohr/autokcalpmol;
            
            UX_private += UX_;
            UX_q_private += UX_ * dq.t();
        }
        delete pot;

        MPI::COMM_WORLD.Allreduce(&Upot_private, &Upot, 1,
                MPI::DOUBLE, MPI::SUM);

        MPI::COMM_WORLD.Allreduce(UX_private.memptr(), UX.memptr(), UX.n_elem,
                MPI::DOUBLE, MPI::SUM);

        MPI::COMM_WORLD.Allreduce(UX_q_private.memptr(), UX_q.memptr(), UX_q.n_elem,
                MPI::DOUBLE, MPI::SUM);
    }
    
    Upot /= NSobol;
    UX /= NSobol;
    UX_q /= NSobol;
    
    return Upot;
}


double SCP1::operator()(vec& q, double kT, mat& Ks, int max_iterations)
{
    double Upot, F;
    mat U;
    vec omega, omegasq, d, UX;

    if (Ks.is_empty()) {
        h2o::PES *pot = getPES();
        getHessian(*pot, q, 1e-3, Ks);
        delete pot;

        massScaleHessian(mass, Ks);
        regtransrot(transrotBasis(q, mass), Ks);    
    }
    

    eig_sym(omegasq, U, Ks);
    omega = sqrt(abs(omegasq));
    d = disp_disp_corr_diag(omega, kT);

    mat isqrtM_U = U;
    isqrtM_U.each_col()  /= sqrt(mass);

    bool finished = false;
    for (int niter = 0; niter < max_iterations; niter++) {
        vec d_old = d;
        vec q_old = q;
 
        mat UX_q;
        Upot = average(q, d, isqrtM_U, UX, UX_q);
        
        for (int i=0; i < UX_q.n_cols; i++) {
            for (int j=0; j < UX_q.n_rows; j++) {
                UX_q(j, i) *= sqrt(mass[i]/mass[j]);
            }
        }
    
        mat KD = UX_q * U;
        for (int i=0; i<Ks.n_cols; i++) {
            if (d[i] != 0.0) {
                KD.col(i) /= d[i];
            }
            else {
                KD.col(i).fill(0.0);
            }
        }
        mat Ks_old = Ks;
        Ks = KD * U.t();
        
        KD = Ks + Ks.t();

        if (niter < min(50, max_iterations) ) {
            Ks = 0.5 * KD;
        }
        else {
            Ks = 0.5 * KD;
        }

        regtransrot(transrotBasis(q, mass), Ks);

        eig_sym(omegasq, U, Ks);
        isqrtM_U = U;
        isqrtM_U.each_col()  /= sqrt(mass);
        
        omega = sqrt(abs(omegasq));
        //omega((omega > 1.0/autocm) % (omega < omega_cutoff)).fill(omega_cutoff);
        //omega(find((omega < omega_cutoff) % (omega > 1e-8) ).fill(omega_cutoff);)

        
        d = disp_disp_corr_diag(omega, kT);
        
        q -= qeps * isqrtM_U * (d % (isqrtM_U.t() * UX));

        finished = prod(abs(q - q_old) <= abs(q_old) * epsrel)
                    * prod(abs(d - d_old) <= abs(d_old) * epsrel);
        
        F  = Upot;
        vec omega_eckart = omega.subvec(6,omega.n_elem - 1);
        if (kT > 0.0) {
            F += kT * sum(log(2.0*sinh(0.5/kT * omega_eckart)))
            - 0.25*sum(omega_eckart/tanh(0.5/kT * omega_eckart)); // - 3*kT;
        }
        else {
            F += 0.25 * sum(omega_eckart);
        }
        if (rank == 0) {
            omega_out << omega.t();
            free_energy_out << F << endl;
            save_for_vladimir("coord.xyz", F, q, Ks);

            cout << "iter: " << niter << endl;
        }
    }
    if (rank == 0) {
        cout << "omega = " << omega.t()*autocm << endl;
        cout << "<U> = " << Upot << endl;
        cout << "<UX> = " << UX.t() << endl;
        cout << "F = " << F * autokcalpmol <<" kcal/mol"<< endl;
    }
    
    return F;
}
