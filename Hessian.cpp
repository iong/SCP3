//
//  Hessian.cpp
//  SCP3
//
//  Created by Ionut Georgescu on 7/15/13.
//
//

#include "Constants.h"
#include "Hessian.h"


void getHessian(h2o::PES& pot, const vec& r_au, double s, mat& H)
{
    int N = r_au.n_rows;
    int Nw = N/9;
    
    vec Vrp(N), Vrm(N);
    H.set_size(N,N);
    
    vec r(r_au * bohr);
    
    for (int i=0; i<N; i++) {
        double ri0 = r[i];
        double V;
        
        r[i] = ri0 + s * bohr;
        V = pot(Nw, r.memptr(), Vrp.memptr());
        
        r[i] = ri0 - s * bohr;
        V = pot(Nw, r.memptr(), Vrm.memptr());
        
        H.col(i) = (Vrp - Vrm) / (2.0*s*autokcalpmol);
        
        r[i] = ri0;
    }
    
    for (int i=1; i<N; i++) {
        H.col(i).subvec(0, i-1) = 0.5 * (H.col(i).subvec(0, i-1)
                                         + H.row(i).subvec(0, i-1).t() );
    }
}


void massScaleHessian(const vec& mass, mat& H)
{
    vec isqrt_mass = sqrt(mass);
    for (int i=0; i<isqrt_mass.n_rows; i++) {
	    isqrt_mass[i] = 1.0/isqrt_mass[i];
    }
    
    for (int i=0; i<H.n_cols; i++) {
        H.col(i) %= isqrt_mass * isqrt_mass[i];
    }
}

mat transrotBasis(const vec& q)
{
    mat U(q.n_rows, 6);
    U.fill(0.0);
    
    for (int i=0; i<q.n_rows; i += 3) {
        U(i  , 0) = 1.0/sqrt(q.n_rows / 3.0);
        U(i+1, 1) = 1.0/sqrt(q.n_rows / 3.0);
        U(i+2, 2) = 1.0/sqrt(q.n_rows / 3.0);
        
        U(i+1, 3) = -q[i+2];
        U(i+2, 3) =  q[i+1];
        
        U(i  , 4) =  q[i+2];
        U(i+2, 4) = -q[i  ];
        
        U(i  , 5) = -q[i+1];
        U(i+1, 5) =  q[i  ];
    }
    
    for (int i=3; i<6; i++) {
        for (int j = 0; j<i; j++) {
            U.col(i) -= U.col(j) * dot(U.col(i), U.col(j));
        }
        
        U.col(i) /= norm(U.col(i), 2);
    }
    
    return U;
}


mat transrotBasis(const vec& q, const vec& mass)
{
    mat U = zeros(q.n_rows, 6);
    vec qs = q % sqrt(mass);
    
    double total_mass = sum(mass)/3.0;
    
    for (int i=0; i<q.n_rows; i += 3) {
        U(i  , 0) = sqrt(mass[i  ] / total_mass);
        U(i+1, 1) = sqrt(mass[i+1] / total_mass);
        U(i+2, 2) = sqrt(mass[i+2] / total_mass);
        
        U(i+1, 3) = -qs[i+2];
        U(i+2, 3) =  qs[i+1];
        
        U(i  , 4) =  qs[i+2];
        U(i+2, 4) = -qs[i  ];
        
        U(i  , 5) = -qs[i+1];
        U(i+1, 5) =  qs[i  ];
    }
    
    for (int i=3; i<6; i++) {
        for (int j = 0; j<i; j++) {
            U.col(i) -= U.col(j) * dot(U.col(i), U.col(j));
        }
        
        U.col(i) /= norm(U.col(i), 2);
    }
    
    return U;
}


void regtransrot(const mat& U, mat &H)
{
        //mat P = eye(q.n_rows, q.n_rows) - U * U.t();
    
    mat HP = H - (H*U) * U.t();
    
    H = HP - U*(U.t() * HP);
}



