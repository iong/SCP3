/*
 *  ScaledMatrixElements.h.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 3/7/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <iostream>
#include <algorithm>

#include "F90.h"
#include "HarmonicOscillatorBasis.h"
#include "ScaledMatrixElements.h"
#include "SimplexIterator.h"


int ScaledMatrixElements::k_index(int k)
{
    return k_fence + k;
}

int ScaledMatrixElements::k2_index(int k)
{
    return k2_fence + k;
}

/**
 1 + 2*Nmodes + (Nmodes - 1 + Nmodes - k)*k/2 + l - k - 1
 */
int ScaledMatrixElements::kl_index(int k, int l)
{
    return kl_fence + k*Nmodes2 - ((k + 2) * (k + 1))/2 + l;
}

int ScaledMatrixElements::k3_index(int k)
{
    return k3_fence + k;
}

int ScaledMatrixElements::k2l_index(int k,int l)
{
    return k2l_fence + 2*Nmodes3*k + 2*l - (k+2)*(k + 1);
}

/**
 \f[
 \sum_{n=Nmodes3-(k-1)-1}^{Nmodes3-1} n*(n-1)/2 +l*(Nmodes3-k-1) - ((l + 2) * (l + 1))/2 + j;
 */
int ScaledMatrixElements::klj_index(int k, int l, int j)
{
    int kth_plane_offset, jl_pos;
    
    kth_plane_offset = (k * (2 + k*k - 3*k*(Nmodes3 - 1) 
                               + 3*(Nmodes3 - 2) * Nmodes3) ) / 6;
    /*
     int Np = Nmodes3 - k - 1;
     int lp = l - k - 1;
     int jp = j - k - 1;
     return klj_fence + jl_offset + lp*Np - ((lp + 2) * (lp + 1))/2 + jp;
     */
    jl_pos = j + (k*k - l*(3 + l) + k*(3 - 2*Nmodes3) + 2*(l - 1)*Nmodes3)/2;
    
    return klj_fence + kth_plane_offset + jl_pos;
}

size_t ScaledMatrixElements::getBasisSize()
{
    size_t Nstates2 = 1 + Nmodes + Nmodes2*(Nmodes2 + 1)/2;

    return Nstates2 + (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6;
}

size_t ScaledMatrixElements::getSubBasisSize(int n)
{
    size_t Nstates = 0;

    switch (n) {
        case 3:
            Nstates += (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6;
        case 2:
            Nstates += Nmodes2*(Nmodes2 + 1)/2;
        case 1:
            Nstates += Nmodes;
        case 0:
            Nstates += 1;
        default:
            break;
    }

    return Nstates;
}

void ScaledMatrixElements::addEpotSingles(mat &M)
{
    int m, k;
    
    vec v(M.n_rows);
    
#pragma omp for schedule(dynamic) nowait
    for (k=0; k<Nmodes; k++) {
        int C = k_index(k);
        int R = 0;
        
        v[R++] = f_k1(k);
        
        for (m=0; m<k; m++) {
            v[R++] = m1_f_k1(m, k);
        }
        
        if (R != C) {
            cerr << "singles: R != k" << endl;
            exit(EXIT_FAILURE);
        }
        
        v[R] = m1_f_m1(k);
        
        M.col(C).subvec(0,R) += v.subvec(0,R);
    }
    
    //cout << "+ " << m1_f_m1(Nmodes-2) << " " << m1_f_m1(Nmodes-1) << endl;
    //cout << "  " << M(Nmodes-1,Nmodes-1) << " " << M(Nmodes,Nmodes) << endl;
    
}

void ScaledMatrixElements::addEpotDoubles(mat &M)
{
    int m, n, k, l;
    
    vec v(M.n_rows);
#pragma omp for schedule(dynamic) nowait
    for (k=0; k<Nmodes2; k++) {
        int R = 0;
        int C = k2_index(k);
        
        v.subvec(0,C).fill(0.0);
        
        v[R++] = f_k2(k);
        
        for (m=0; m<Nmodes; m++) {
            v[R++] = m1_f_k2(m, k);
        }
        m=k;
        v[k_index(m)] = m1_f_m2(k);
        
        for (m=0; m<k; m++) {
            v[R++] = m2_f_k2(m, k);
        }
        
        if (R != C) {
            cerr << "doubles: R != C" << endl;
            exit(EXIT_FAILURE);
        }
        
        m=k;
        v[R] = m2_f_m2(k);
        
        M.col(C).subvec(0,R) += v.subvec(0,R);
    }
    
    SimplexIterator<2> s2(Nmodes2);
#pragma omp for schedule(dynamic) nowait
    for (int i=0; i < (s2.volume() + 1) / 2; i++) {
        s2 = i;
        for (int i2=0; i2 < 2; i2++) {
            k = s2.index[0];
            l = s2.index[1];

            int C = kl_index(k,l);
            v.subvec(0, C).fill(0.0);
            
            int R = 0;
            v[R++] = f_k1_l1(k, l);
            
            for (m=0; m<Nmodes; m++) {
                v[R++] = m1_f_k1_l1(m, k, l);
            }
            m=k;
            v[k_index(m)] = m1_f_m1_l1(k, l);
            
            m=l;
            v[k_index(m)] = m1_f_m1_l1(l, k);
            
            for (m=0; m<Nmodes2; m++) {
                v[R++] = m2_f_k1_l1(m, k, l);
            }
            m=k;
            v[k_index(m)] = m2_f_m1_l1(k, l);
            
            m=l;
            v[k_index(m)] = m2_f_m1_l1(l, k);
            
            for (m=0; m<k; m++) {
                for (n=m+1; n<Nmodes2; n++) {
                    v[R++] = m1_n1_f_k1_l1(m, n, k, l);
                }
            }
            
            n=k;
            for (m=0; m<k; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1(k, m, l);
            }
            
            n=l;
            for (m=0; m<k; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1(l, m, k);
            }
            
            m=k;
            for (n=m+1; n<l; n++) {
                v[R++] = m1_n1_f_m1_l1(k, n, l);
            }
            
            if (R != C) {
                cerr << "doubles_mn: R != C" << endl;
                exit(EXIT_FAILURE);
            }
            
            v[R] = m1_n1_f_m1_n1(k, l);
            
            M.col(C).subvec(0, R) += v.subvec(0,R);

            if ( (2*i + 1) == s2.volume() ) break;

            s2 = s2.volume() - 1 - i;
        }
    }
}

void ScaledMatrixElements::addEpotTriples(mat &M)
{
    int m, n, p, k, l, j;
    
    if (Nmodes3 == 0) return;
    
    vec v(M.n_rows), v2(M.n_rows);
    
#pragma omp for schedule(dynamic) nowait
    for (k=0; k<Nmodes3; k++) {
        int C = k3_index(k);
        int R = 0;
        
        v.subvec(0, C).fill(0.0);
        
        v[R++] = f_k3(k);
        
        for (m=0; m<Nmodes; m++) {
            v[R++] = m1_f_k3(m, k);
        }
        m=k;
        v[k_index(m)] = m1_f_m3(m);
        
        for (m=0; m<Nmodes2; m++) {
            v[R++] = m2_f_k3(m, k);
        }
        m=k;
        v[k_index(m)] = m2_f_m3(m);
        
        for (m=0; m<Nmodes2; m++) {
            for (n=m+1; n<Nmodes2; n++) {
                v[R++] = m3_f_k1_l1(k, m, n);
            }
        }
        
        m=k;
        for (n=m+1; n<Nmodes2; n++) {
            v[kl_index(m, n)] = m3_f_m1_l1(m, n);
        }
        
        n=k;
        for (m=0; m<n; m++) {
            v[kl_index(m, n)] = m3_f_m1_l1(n, m);
        }
        
        for (m=0; m<k; m++) {
            v[R++] = m3_f_k3(m, k);
        }
        
        v[R] = m3_f_m3(k);
        
        if (R != k3_index(k)) {
            cerr << "R != k3_index(k)" << endl;
            exit(EXIT_FAILURE);
        }
        
        M.col(R).subvec(0,R) += v.subvec(0, R);
    }

    SimplexIterator<2> s2(Nmodes3);
#pragma omp for schedule(dynamic) nowait
    for (int i=0; i< (s2.volume() + 1) / 2; i++) {
        s2 = i;
        for (int i2=0; i2 < 2; i2++) {
            k = s2.index[0];
            l = s2.index[1];
            
            int C = k2l_index(k, l);
            
            v.subvec(0, C).fill(0.0);
            v2.subvec(0, C).fill(0.0);
            
            int R = 0;
            v[R] = f_k2_l1(k, l);
            v2[R++] = f_k2_l1(l, k);
            
            for (m=0; m<Nmodes; m++) {
                v[R] = m1_f_k2_l1(m, k, l);
                v2[R++] = m1_f_k2_l1(m, l, k);
            }
            m=k;
            v[k_index(m)] = m1_f_m2_l1(m, l);
            v2[k_index(m)] = m1_f_m1_l2(m, l);//unique, no coordinate swap possible!
            
            m=l;
            v[k_index(m)] = m1_f_m1_l2(m, k);
            v2[k_index(m)] = m1_f_m2_l1(m, k); // unique again.
            
            for (m=0; m<Nmodes2; m++) {
                v[R] = m2_f_k2_l1(m, k, l);
                v2[R++] = m2_f_k1_l2(m, k, l);
            }
            m=k;
            v[k2_index(m)] = m2_f_m2_l1(m, l);
            v2[k2_index(m)] = m2_f_m1_l2(m, l);
            
            m=l;
            v[k2_index(m)] = m2_f_m1_l2(m, k);
            v2[k2_index(m)] = m2_f_m2_l1(m, k);
            
            R = kl_fence;
            for (m=0; m<Nmodes2; m++) {
                for (n=m+1; n<Nmodes2; n++) {
                    v[R] = m1_n1_f_k2_l1(m, n, k, l);
                    v2[R++] = m1_n1_f_k1_l2(m, n, k, l);
                }
            }
            
            m=k;
            for (n=m+1; n<Nmodes2; n++) {
                R = kl_index(m, n);
                v[R] = m1_n1_f_m2_l1(m, n, l);
                v2[R] = m1_n1_f_m1_l2(m, n, l);
            }
            
            m=l;
            for (n=l+1; n<Nmodes2; n++) {
                R = kl_index(l, n);
                v[R] = m1_n1_f_m1_l2(l, n, k);
                v2[R] = m1_n1_f_m2_l1(l, n, k);
            }
            
            n=k;
            for (m=0; m<n; m++) {
                R = kl_index(m, n);
                v[R] = m1_n1_f_m2_l1(n, m, l);
                v2[R] = m1_n1_f_m1_l2(n, m, l);
            }
            
            n=l;
            for (m=0; m<n; m++) {
                R = kl_index(m, n);
                v[R] = m1_n1_f_m1_l2(n, m, k);
                v2[R] = m1_n1_f_m2_l1(n, m, k);
            }
            m=k; n=l;
            R = kl_index(m, n);
            v[R] = m1_n1_f_m2_n1(m, n);
            v2[R] = m1_n1_f_m1_n2(m, n);
            
            R = k3_fence;
            for (m=0; m<Nmodes3; m++) {
                v[R] = m3_f_k2_l1(m, k, l);
                v2[R++] = m3_f_k1_l2(m, k, l);
            }
            
            m=k;
            v[k3_index(m)] = m3_f_m2_l1(k, l);
            v2[k3_index(m)] = m3_f_m1_l2(k, l);
            
            m=l;
            v[k3_index(m)] = m3_f_m1_l2(l, k);
            v2[k3_index(m)] = m3_f_m2_l1(l, k);
            
            R=k2l_fence;
            for (m=0; m<k; m++) {
                for (n=m+1; n<Nmodes3; n++) {
                    v[R] = m2_n1_f_k2_l1(m, n, k, l);
                        //v2[R++] = m2_n1_f_k1_l2(m, n, k, l);
                    v2[R++] = m1_n2_f_k2_l1(k, l, m, n);
                    
                    v[R] = m1_n2_f_k2_l1(m, n, k, l);
                    v2[R++] = m1_n2_f_k1_l2(m, n, k, l);
                }
            }
            
            n=k;
            for (m=0; m<k; m++) {
                R = k2l_index(m, n);
                v[R] = m1_n2_f_m2_l1(k, m, l);
                v2[R++] = m1_n2_f_m1_l2(k, m, l);
                
                
                v[R] = m2_n1_f_m2_l1(k, m, l);
                v2[R] = m1_n2_f_m2_l1(k, l, m);
            }
            
            n=l;
            for (m=0; m<k; m++) {
                R = k2l_index(m, n);
                v[R] = m1_n2_f_m1_l2(l, m, k);
                v2[R++] = m1_n2_f_m2_l1(l, m, k);
                
                
                    //v[R] = m2_n1_f_m1_l2(l, m, k);
                v[R] = m1_n2_f_m2_l1(l, k, m);
                v2[R] = m2_n1_f_m2_l1(l, m, k);
            }
            
            
            m=k;
            for (n=m+1; n<l; n++) {
                R = k2l_index(m, n);
                v[R] = m2_n1_f_m2_l1(k, n, l);
                v2[R++] = m1_n2_f_m2_l1(k, l, n);
                
                
                v[R] = m1_n2_f_m2_l1(k, n, l);
                v2[R] = m1_n2_f_m1_l2(k, n, l);
            }
            
            n=l;
            R = k2l_index(m, n);
            v[R] = m2_n1_f_m2_n1(m, n);
            v2[R++] = m1_n2_f_m2_n1(m, n);
            
            v2[R] = m1_n2_f_m1_n2(m, n);
            
            if ( R!= C+1) {
                cerr << "R != C\n";
                exit(EXIT_FAILURE);
            }
            M.col(C).subvec(0,C) += v.subvec(0, C);
            M.col(C+1).subvec(0,C+1) += v2.subvec(0, C+1);

            if ( (2*i + 1) == s2.volume() ) break;

            s2 = s2.volume() - 1 - i;
        }
    }
    
    SimplexIterator<3> s3(Nmodes3);
#pragma omp for schedule(dynamic)
    for (int i=0; i< (s3.volume() + 1) / 2; i++) {
        s3 = i;
        for (int i2=0; i2 < 2; i2++) {
            k = s3.index[0];
            l = s3.index[1];
            j = s3.index[2];

            int C = klj_index(k, l, j);
            int R=0;
            
            v.subvec(0, C).fill(0.0);
            
            v[R++] = f_k1_l1_j1(k, l, j);
            
            for (m=0; m<Nmodes; m++) {
                v[R++] = m1_f_k1_l1_j1(m, k, l, j);
            }
            v[k+1] = m1_f_m1_l1_j1(k, l, j);
            v[l+1] = m1_f_m1_l1_j1(l, k, j);
            v[j+1] = m1_f_m1_l1_j1(j, k, l);
            
            for (m=0; m<Nmodes2; m++) {
                v[R++] = m2_f_k1_l1_j1(m, k, l, j);
            }
            v[k2_index(k)] = m2_f_m1_l1_j1(k, l, j);
            v[k2_index(l)] = m2_f_m1_l1_j1(l, k, j);
            v[k2_index(j)] = m2_f_m1_l1_j1(j, k, l);
            
            for (m=0; m<Nmodes2; m++) {
                for (n=m+1; n<Nmodes2; n++) {
                    v[R++] = m1_n1_f_k1_l1_j1(m, n, k, l, j);
                }
            }
            
            m=k;
            for (n=m+1; n<Nmodes2; n++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, l, j);
            }
            
            m=l;
            for (n=m+1; n<Nmodes2; n++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, k, j);
            }
            
            m=j;
            for (n=m+1; n<Nmodes2; n++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, k, l);
            }
            
            n=k;
            for (m=0; m<n; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(n, m, l, j);
            }
            
            n=l;
            for (m=0; m<n; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(n, m, k, j);
            }
            
            n=j;
            for (m=0; m<n; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(n, m, k, l);
            }
            
            m=k; n=l;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(k, l, j);
            n=j;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(k, j, l);
            
            
            m=l; n=j;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(l, j, k);
            
            
            for (m=0; m<Nmodes3; m++) {
                v[R++] = m3_f_k1_l1_j1(m, k, l, j);
            }
            v[k3_index(k)] = m3_f_m1_l1_j1(k, l, j);
            v[k3_index(l)] = m3_f_m1_l1_j1(l, k, j);
            v[k3_index(j)] = m3_f_m1_l1_j1(j, k, l);
            
            
            for (m=0; m<Nmodes3; m++) {
                for (n = m+1; n<Nmodes3; n++) {
                    v[R++] = m2_n1_f_k1_l1_j1(m, n, k, l, j);
                    v[R++] = m1_n2_f_k1_l1_j1(m, n, k, l, j);
                }
            }
            
            m=k;
            for (n=m+1; n<Nmodes3; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(m, n, l, j);
                v[R++] = m1_n2_f_m1_l1_j1(m, n, l, j);
            }
            
            m=l;
            for (n=m+1; n<Nmodes3; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(m, n, k, j);
                v[R++] = m1_n2_f_m1_l1_j1(m, n, k, j);
            }
            
            m=j;
            for (n=m+1; n<Nmodes3; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(m, n, k, l);
                v[R++] = m1_n2_f_m1_l1_j1(m, n, k, l);
            }
            
            n=k;
            for (m=0; m<n; m++) {
                R = k2l_index(m, n);
                v[R++] = m1_n2_f_m1_l1_j1(n, m, l, j);
                v[R++] = m2_n1_f_m1_l1_j1(n, m, l, j);
            }
            
            n=l;
            for (m=0; m<n; m++) {
                R = k2l_index(m, n);
                v[R++] = m1_n2_f_m1_l1_j1(n, m, k, j);
                v[R++] = m2_n1_f_m1_l1_j1(n, m, k, j);
            }
            
            n=j;
            for (m=0; m<n; m++) {
                R = k2l_index(m, n);
                v[R++] = m1_n2_f_m1_l1_j1(n, m, k, l);
                v[R++] = m2_n1_f_m1_l1_j1(n, m, k, l);
            }
            
            m=k; n=l;
            R = k2l_index(m, n);
            v[R++] = m2_n1_f_m1_n1_j1(m, n, j);
            v[R  ] = m1_n2_f_m1_n1_j1(m, n, j);
            
            n=j;
            R = k2l_index(m, n);
            v[R++] = m2_n1_f_m1_n1_j1(m, n, l);
            v[R  ] = m1_n2_f_m1_n1_j1(m, n, l);
            
            m=l; n=j;
            R = k2l_index(m, n);
            v[R++] = m2_n1_f_m1_n1_j1(m, n, k);
            v[R  ] = m1_n2_f_m1_n1_j1(m, n, k);
            
            
            R = klj_fence;
            for (m=0; m<k; m++) {
                for (n=m+1; n<Nmodes3; n++) {
                    for (p=n+1; p<Nmodes3; p++) {
                        v[R++] = m1_n1_p1_f_k1_l1_j1(m, n, p, k, l, j);
                    }
                }
            }
            
            /* m<k, n=k */
            n=k;
            for (m=0; m<n; m++) {
                for (p=n+1; p<Nmodes3; p++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, l, j);
                }
            }
            n=l;
            for (m=0; m<n; m++) {
                for (p=n+1; p<Nmodes3; p++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, k, j);
                }
            }
            n=j;
            for (m=0; m<n; m++) {
                for (p=n+1; p<Nmodes3; p++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, k, l);
                }
            }
            
            p=k;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, l, j);
                }
            }
            p=l;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, k, j);
                }
            }
            p=j;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, k, l);
                }
            }
            
            n=k; p=l;
            for (m=0; m<k; m++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, j);
            }
            p=j;
            for (m=0; m<k; m++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, l);
            }
            n=l; p=j;
            for (m=0; m<k; m++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, k);
            }
            
            m=k;
            for (n=m+1; n<l; n++) {
                R = klj_index(m, n, n+1);
                for (p=n+1; p<Nmodes3; p++) {
                    v[R++] = m1_n1_p1_f_m1_l1_j1(m, n, p, l, j);
                }
            }
            
            p=l;
            for (n=m+1; n<p; n++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(m, p, n, j);
            }
            p=j;
            for (n=m+1; n<p; n++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(m, p, n, l);
            }
            
            n=l;
            R = klj_index(m, n, n+1);
            for (p=n+1; p<j; p++) {
                v[R++] = m1_n1_p1_f_m1_n1_j1(m, n, p, j);
            }
            
            p=j;
            v[R] = m1_n1_p1_f_m1_n1_p1(m, n, p);
            
            if ( R!= C) {
                cerr << "R != C\n";
                exit(EXIT_FAILURE);
            }
            
            M.col(C).subvec(0, C) += v.subvec(0, C);

            if ( (2*i + 1) == s3.volume() ) break;

            s3 = s3.volume() - 1 - i;
        }
    }
}

void ScaledMatrixElements::addEpot(const vec &q_, double V_, const vec& Vq_, mat &M)
{
    q = q_;
    V = V_;
    Vq = Vq_;
    
    M(0,0) += V;
    
#pragma omp parallel
    {
        addEpotSingles(M);
        addEpotDoubles(M);
        addEpotTriples(M);
    }
}

void ScaledMatrixElements::test_index()
{
    cout << klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1) << " ";
    cout << 1 + Nmodes + Nmodes2*(Nmodes2 + 1)/2 
                + (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6;
    cout << endl;
}

void ScaledMatrixElements::addHODiagonal(mat &M)
{
    int k, l, j, C;
    double E0 = 0.5*sum(omega);
    
    M(0,0) += E0;
    C = k_fence;
    for (int k=0; k<Nmodes; k++, C++) {
        M(C,C) += E0+omega(k);
    }
    
    for (k=0; k<Nmodes2; k++, C++) {
        M(C,C) += E0+2*omega(k);
    }
    
    for (k=0; k<Nmodes2; k++) {
        for (l=k+1; l<Nmodes2; l++, C++) {
            M(C,C) += E0 + omega(k) + omega(l);
        }
    }
    
    for (k=0; k<Nmodes3; k++, C++) {
        M(C,C) += E0 + 3.0*omega(k);
    }
    
    for (k=0; k<Nmodes3; k++) {
        for (l=k+1; l<Nmodes3; l++) {
            M(C,C) += E0 + 2.0*omega(k) +       omega(l);
            C++;
            M(C,C) += E0 +     omega(k) + 2.0 * omega(l);
            C++;
        }
    }
    
    for (k=0; k<Nmodes3; k++) {
        for (l=k+1; l<Nmodes3; l++) {
            for (j=l+1; j<Nmodes3; j++, C++) {
                M(C,C) += E0 + omega(k) + omega(l) + omega(j);
            }
        }
    }
    /*
     if (C-1 != klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1) ) {
     cerr << "HO diagonal term ended in wrong col number." << endl;
     cerr << C-1 << " != " << klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1);
     cerr << endl;
     exit(EXIT_FAILURE);
     }
     */
}


mat ScaledMatrixElements::transitionDipole(vec& charges, mat& MUa, mat& C)
{
    mat Q(3,MUa.n_rows);
    Q.fill(0);
    for (int i=0; i<Q.n_cols; i += 3) {
        Q(0, i  ) = charges[i/3];
        Q(1, i+1) = charges[i/3];
        Q(2, i+2) = charges[i/3];
    }

    mat mu_m_m2 = Q*MUa;
    mat mu_0_m = mu_m_m2 / sqrt(2.0);
    
    mat mu(3,C.n_cols-1);    /*
    for (int f=1; f<C.n_cols; f++) {
        vec d;
        
        d = mu_0_m* C(0,0)*C(span(1,Nmodes),f);
        d += mu_0_m* C(span(1,Nmodes), 0)*C(0,f);

        for (m=0; m<Nmodes2; m++) {
            
            d+= mu_0_m2.col(m) * ( C(k_index(m), 0) % C(k2_index(m), f) 
                                  + C(k2_index(m), 0) % C(k_index(m), f) );
            double t0 = 0.0
            for (n=0; n<m; n++) {
                t0 += C(k_index(n), 0) * C(kl_index(n, m), f);
                t0 += C(kl_index(n, m), 0) * C(k_index(n), f);
            }
            for (n=m+1; n<Nmodes2; n++) {
                t0 += C(k_index(n), 0) * C(kl_index(m, n), f);
                t0 += C(kl_index(m, n), 0) * C(k_index(n), f);
            }
            
            d += mu_0_m.col(m)*t0;
        }
        
        for (m=0; m<Nmodes3; m++) {
            double t1 = 0.0, t2 = 0.0, t3 = 0.0;
            t3 = C(k2_index(m), 0) * C(k3_index(m), f)
                + C(k3_index(m), 0) * C(k2_index(m), f);
            
            for (n=0; n<m; n++) {
                t1 += C(k2_index(n), 0) * C(k2l_index(n, m), f);
                t2 += C(kl_index(n, m), 0) * C(k2l_index(n, m)+1, f);
                
                t2 += C(k2l_index(n, m), 0) * C(k2_index(n), f);
                t2 += C(k2l_index(n, m)+1, 0) * C(kl_index(n, m), f);
                
                for (p=m+1; p<Nmodes3; p++) {
                    t1 += C(kl_index(n, p), 0) * C(klj_index(n, m, p), f);
                    t1 += C(klj_index(n, m, p), 0) * C(kl_index(n, p), f);
                }
                for (p=n+1; p<m; p++) {
                    t1 += C(kl_index(n, p), 0) * C(klj_index(n, p, m), f);
                    t1 += C(klj_index(n, p, m), 0) * C(kl_index(n, p), f);
                }
            }
            
            for (n=m+1; n<Nmodes3; n++) {
                t1 += C(k2_index(n), 0) * C(k2l_index(m, n)+1, f);
                
                t2 += C(kl_index(m, n), 0) * C(k2l_index(m, n), f);
                t2 += C(k2l_index(m, n), 0) * C(kl_index(m, n), f);
                
                t1 += C(k2l_index(m, n)+1, 0) * C(k2_index(n), f);
                
                for (p=n+1; p<Nmodes3; p++) {
                    t1 += C(kl_index(n, p), 0) * C(klj_index(m, n, p), f);
                    t1 += C(klj_index(m, n, p), 0) * C(kl_index(n, p), f);
                }
            }
            
            d += mu_0_m.col(m)*t1 + mu_m_m2.col(m)*t2 + mu_m2_m3.col(m)*t3;
        }
        
        mu.col(f-1) = d;
    }
     */
    
    return mu;
}

vec ScaledMatrixElements::get_bra(const vec& q)
{
    int m, n, p;
    
    int R = 0;
    
    vec bra(getBasisSize());
    
    bra[R++] = 1.0;
    
    
    
    for (m=0; m<Nmodes; m++) {
        bra[R++] = ho_basis[1] (q[m]);
    }
    
    
    for (m=0; m<Nmodes2; m++) {
        bra[R++] = ho_basis[2] (q[m]);
    }
    
    
    for (m=0; m<Nmodes2; m++) {
        double bra_m = ho_basis[1](q[m]);
        
        for (n=m+1; n<Nmodes2; n++) {
            bra[R++] = bra_m * ho_basis[1] (q[n]);
        }
    }
    
    
    for (m=0; m<Nmodes3; m++) {
        bra[R++] = ho_basis[3] (q[m]);
    }
    
    
    for (m=0; m<Nmodes3; m++) {
        double bra_m  = ho_basis[1](q[m]);
        double bra_m2 = ho_basis[2](q[m]);
        
        for (n=m+1; n<Nmodes3; n++) {
            bra[R++] = bra_m2 * ho_basis[1](q[n]);
            bra[R++] = bra_m  * ho_basis[2](q[n]);
        }
    }
    
    
    for (m=0; m<Nmodes3; m++) {
        double bra_m = ho_basis[1]( q[m] );
        
        for (n=m+1; n<Nmodes3; n++) {
            double bra_mn = bra_m * ho_basis[1]( q[n] );
            
            for (p=n+1; p<Nmodes3; p++) {
                bra[R++] = bra_mn * ho_basis[1]( q[p] );
            }
        }
    }
    
    if ( R != bra.n_rows) {
        cerr << "R != M.n_rows";
        exit(EXIT_FAILURE);
    }
    
    return bra;
}