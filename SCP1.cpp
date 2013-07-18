//
//  SCP1.cpp
//  SCP3
//
//  Created by Ionut Georgescu on 7/15/13.
//
//

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

    d ( find(abs(omega) < omega_cutoff) ).fill(0.0);

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

#pragma omp parallel
    {
        h2o::PES *pot = getPES();
        double Upot_private;
        vec UX_private, UX_(dim), y(dim), dq(dim);
        mat UX_q_private;
        
        Upot_private = 0.0;
        UX_private.zeros(dim);
        UX_q_private.zeros(dim, dim);
        
#pragma omp for schedule(static) nowait
        for (int i=0; i < NSobol; i++) {
            long long current_skip = sobol_skip + i;
            sobol::std_normal(y.n_rows, &current_skip, y.memptr());
            
            y %= sqrt(d);
            
            dq = isqrtM_U * y;
            y = bohr * (q + dq);
            Upot_private += (*pot)(nw, y.memptr(), UX_.memptr()) / autokcalpmol;
            UX_ *= bohr/autokcalpmol;
            
            UX_private += UX_;
            UX_q_private += UX_ * dq.t();
        }
        
#pragma omp critical
        {
            Upot += Upot_private;
            UX += UX_private;
            UX_q += UX_q_private;
        }
        delete pot;
    }
    
    Upot /= NSobol;
    UX /= NSobol;
    UX_q /= NSobol;
    
    return Upot;
}


double SCP1::operator()(vec& q, double kT, mat& Ks)
{
    double Upot, F;
    mat U;
    vec omega, omegasq, d, UX;

    h2o::PES *pot = getPES();
    getHessian(*pot, q, 1e-3, Ks);
    delete pot;
    
    massScaleHessian(mass, Ks);
    regtransrot(transrotBasis(q, mass), Ks);    

    eig_sym(omegasq, U, Ks);
    omega = sqrt(abs(omegasq));
    d = disp_disp_corr_diag(omega, kT);

    mat isqrtM_U = U;
    isqrtM_U.each_col()  /= sqrt(mass);

    int niter = 0;
    for (bool finished = false; !finished; niter++) {
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
        Ks = KD * U.t();
        
        KD = Ks + Ks.t();
        Ks = 0.5 * KD;

        regtransrot(transrotBasis(q, mass), Ks);

        eig_sym(omegasq, U, Ks);
        isqrtM_U = U;
        isqrtM_U.each_col()  /= sqrt(mass);
        
        omega = sqrt(abs(omegasq));
        
        omega_out << omega.t();
        
        d = disp_disp_corr_diag(omega, kT);
        
        q -= qeps * isqrtM_U * (d % (isqrtM_U.t() * UX));

        finished = prod(abs(q - q_old) <= abs(q_old) * epsrel)
                    * prod(abs(d - d_old) <= abs(d_old) * epsrel);
        
         F  = Upot;
        if (kT > 0.0) {
            F += kT * sum(log(2.0*sinh(0.5/kT * omega)))
            - 0.25*sum(omega/tanh(0.5/kT * omega)); // - 3*kT;
        }
        else {
            F += 0.25 * sum(omega);
        }
        free_energy_out << F << endl;
    }
    omega_out << endl;
    save_for_vladimir("coord.xyz", F, q, Ks);
    
    omega.shed_rows(0, 5);
    

    cout << "omega = " << omega.t()*autocm << endl;
    cout << "<U> = " << Upot << endl;
    cout << "<UX> = " << UX.t() << endl;
    cout << "F = " << F * autokcalpmol <<" kcal/mol"<< endl;
    
    return F;
}
