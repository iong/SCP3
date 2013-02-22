/*
 *  ScaledMatrixElements.h
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 2/5/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */


#ifndef __SCALED_MATRIX_ELEMENTS_H__
#define __SCALED_MATRIX_ELEMENTS_H__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

double sqr(double x)
{
    return x*x;
}

class ScaledMatrixElements {
    double V, *Vq, *q;
    int Nmodes, Nmodes2, Nmodes3;
    vec &omega;
    
    int k_fence, k2_fence, kl_fence, k2l_fence, k3_fence, klj_fence;
    
public:
    ScaledMatrixElements(vec &omega_, int Nmodes2_ = 0, int Nmodes3_ = 0) :
    omega(omega_), Nmodes2(Nmodes2_), Nmodes3(Nmodes3_)
    {
        Nmodes = omega.n_rows;
        k_fence = 1;
        k2_fence = k_fence + Nmodes;
        kl_fence = k2_fence + Nmodes2;
        k3_fence = kl_fence + (Nmodes2*(Nmodes2-1)) / 2;
        k2l_fence = k3_fence + Nmodes3;
        klj_fence = k2l_fence + Nmodes3*(Nmodes3-1);
        
        test_index();
    }
    
private:
//#include "DoubleExcitations.h"
#include "TripleExcitations.h"
    
    int k_index(int k)
    {
        return k_fence + k;
    }
    
    int k2_index(int k)
    {
        return k2_fence + k;
    }
    
    /**
     1 + 2*Nmodes + (Nmodes - 1 + Nmodes - k)*k/2 + l - k - 1
     */
    int kl_index(int k, int l)
    {
        return kl_fence + k*Nmodes2 - ((k + 2) * (k + 1))/2 + l;
    }
    
    int k3_index(int k)
    {
        return k3_fence + k;
    }
    
    int k2l_index(int k,int l)
    {
        return k2l_fence + 2*Nmodes3*k + 2*l - (k+2)*(k + 1);
    }
    
    /**
     \f[
     \sum_{n=Nmodes3-(k-1)-1}^{Nmodes3-1} n*(n-1)/2 +l*(Nmodes3-k-1) - ((l + 2) * (l + 1))/2 + j;
     */
    int klj_index(int k, int l, int j)
    {
        int jl_offset = (k*(2 + k*k - 3*k*(Nmodes3 - 1) + 3*(Nmodes3 - 2)*Nmodes3))/6;
        return klj_fence + jl_offset + l*(Nmodes3-k-1) - ((l + 2) * (l + 1))/2 + j;
    }
    
    
public:
    void addEpotSingles(vec &q_, double V_, vec& Vq_, mat &M)
    {
        int m, k;
        
        q = q_.memptr();
        V = V_;
        Vq = Vq_.memptr();
        
        M(0,0) += V;
        
        for (k=0; k<Nmodes; k++) {
            double v = f_k1(k);
            M(0,k+1) += v;
            M(k+1,0) += v;
        }
        
        for (m=0; m<Nmodes; m++) {
            M(m+1,m+1) += m1_f_m1(m);
            for (k=m+1; k<Nmodes; k++) {
                double v = m1_f_k1(m, k);
                M(m+1,k+1) += v;
                M(k+1,m+1) += v;
            }
        }
        
        //cout << "+ " << m1_f_m1(Nmodes-2) << " " << m1_f_m1(Nmodes-1) << endl;
        //cout << "  " << M(Nmodes-1,Nmodes-1) << " " << M(Nmodes,Nmodes) << endl;
        
    }
    
    void addEpotDoubles(vec &q_, double V_, vec& Vq_, mat &M)
    {
        int m, n, k, l;
        
        addEpotSingles(q_, V_, Vq_, M);
#pragma omp parallel private(m, n, k, l)
        {
#pragma omp for schedule(static) nowait
            for (k=0; k<Nmodes2; k++) {
                double v = f_k2(k);
                int k2pos = k2_index(k);
                M(0, k2pos) += v;
                M(k2pos,0) += v;
            }
            
#pragma omp for schedule(static) nowait
            for (m=0; m<Nmodes; m++) {
                for (k=0; k<Nmodes2; k++) {
                    int k2pos = k2_index(k);
                    double v = m1_f_k2(m, k);
                    
                    if (m==k) {
                        v = m1_f_m2(m);
                    }
                    
                    M(m+1, k2pos) += v;
                    M(k2pos, m+1) += v;
                }
            }
            
#pragma omp for schedule(static) nowait
            for (m=0; m<Nmodes2; m++) {
                int m2pos = k2_index(m);
                M(m2pos, m2pos) += m2_f_m2(m);
                for (k=m+1; k<Nmodes2; k++) {
                    int k2pos = k2_index(k);
                    double v = m2_f_k2(m, k);
                    
                    M(m2pos, k2pos) += v;
                    M(k2pos, m2pos) += v;
                }
            }
            
#pragma omp for schedule(static) nowait
            for (k=0; k<Nmodes2; k++) {
                for (l=k+1; l<Nmodes2; l++) {
                    int klpos = kl_index(k,l);
                    
                    double v = f_k1_l1(k, l);
                    M(0, klpos) += v;
                    M(klpos, 0) += v;
                }
            }
            
#pragma omp for schedule(static) nowait
            for (m=0; m<Nmodes; m++) {
                for (k=0; k<Nmodes2; k++) {
                    for (l=k+1; l<Nmodes2; l++) {
                        int klpos = kl_index(k,l);
                        double v = m1_f_k1_l1(m, k, l);
                        
                        if (m==l) {
                            v = m1_f_m1_l1(m, k);
                        }
                        else if (m==k) {
                            v = m1_f_m1_l1(m, l);
                        }
                        
                        M(m+1, klpos) += v;
                        M(klpos, m+1) += v;
                    }
                }
            }
            
#pragma omp for schedule(static) nowait
            for (m=0; m<Nmodes; m++) {
                int m2pos = k2_index(m);
                
                for (k=0; k<Nmodes2; k++) {
                    for (l=k+1; l<Nmodes2; l++) {
                        int klpos = kl_index(k,l);
                        double v = m2_f_k1_l1(m, k, l);
                        
                        if (m==l) {
                            v = m2_f_m1_l1(m, k);
                        }
                        else if (m==k) {
                            v = m2_f_m1_l1(m, l);
                        }
                        
                        M(m2pos, klpos) += v;
                        M(klpos, m2pos) += v;
                    }
                }
            }
            
#pragma omp for schedule(dynamic,2)
            for (m=0; m<Nmodes2; m++) {
                for (n=m+1; n<Nmodes2; n++) {
                    k = m;
                    l = n;
                    
                    int mnpos = kl_index(m, n);
                    M(mnpos, mnpos) += m1_n1_f_m1_n1(m, n);
                    
                    for (l=n+1; l<Nmodes2; l++) {
                        int klpos = kl_index(k, l);
                        double v = m1_n1_f_m1_l1(m, n, l);
                        
                        M(mnpos, klpos) += v;
                        M(klpos, mnpos) += v;
                    }
                    
                    for (k=m+1; k<Nmodes2; k++) {
                        for (l=k+1; l<Nmodes2; l++) {
                            double v =  m1_n1_f_k1_l1(m, n, k, l);
                            int klpos = kl_index(k, l);
                            
                            if (l==n) {
                                v = m1_n1_f_m1_l1(n, m, k);
                            }
                            
                            M(mnpos, klpos) += v;
                            M(klpos, mnpos) += v;
                        }
                    }
                }
            }
        }
    }
    
    void addEpotTriples(vec &q_, double V_, vec& Vq_, mat &M)
    {
        int m, n, p, k, l, j;
        
        addEpotDoubles(q_, V_, Vq_, M);
        
        vec v(M.n_rows);
        for (k=0; k<Nmodes3; k++) {
            int R = 0;
            v[R] = f_k3(k);
            
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
        
            R = k3_fence;
            for (m=0; m<k; m++) {
                v[R++] = m3_f_k3(m, k);
            }
        
            v[R] = m3_f_m3(k);
        
            if (R != k3_index(k)) {
                cerr << "R != k3_index(k)" << endl;
                exit(EXIT_FAILURE);
            }

            M.col(R).subvec(0,R-1) += v.subvec(0, R-1);
        }
    
        vec v2(M.n_rows);
        for (k=0; k<Nmodes3; k++) {
            for (l=k+1; l<Nmodes3; l++) {
                int C = k2l_index(k, l);
                
                v .fill(0.0);
                v2.fill(0.0);
                
                int R = 0;
                v[R] = f_k2_l1(k, l);
                v2[R++] = f_k1_l2(k, l);
                
                for (m=0; m<Nmodes; m++) {
                    v[R] = m1_f_k2_l1(m, k, l);
                    v2[R++] = m1_f_k1_l2(m, k, l);
                }
                m=k;
                v[k_index(m)] = m1_f_m2_l1(m, l);
                v2[k_index(m)] = m1_f_m1_l2(m, l);
                
                m=l;
                v[k_index(m)] = m1_f_m2_l1(m, k);
                v2[k_index(m)] = m1_f_m1_l2(m, k);
                
                for (m=0; m<Nmodes2; m++) {
                    v[R] = m2_f_k2_l1(m, k, l);
                    v2[R++] = m2_f_k1_l2(m, k, l);
                }
                m=k;
                v[k2_index(m)] = m2_f_m2_l1(m, l);
                v2[k2_index(m)] = m2_f_m1_l2(m, l);
                
                m=l;
                v[k2_index(m)] = m2_f_m2_l1(m, k);
                v2[k2_index(m)] = m2_f_m1_l2(m, k);
                
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
                for (n=m+1; n<Nmodes2; n++) {
                    R = kl_index(m, n);
                    v[R] = m1_n1_f_m2_l1(m, n, k);
                    v2[R] = m1_n1_f_m1_l2(m, n, k);
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
                    v[R] = m1_n1_f_m2_l1(n, m, k);
                    v2[R] = m1_n1_f_m1_l2(n, m, k);
                }
                m=k; n=l;
                R = kl_index(m, n);
                v[R] = m1_n1_f_m2_n1(m, n);
                v2[R] = m1_n1_f_m1_n2(m, n);
                
                for (m=0; m<Nmodes3; m++) {
                    v[R] = m3_f_k2_l1(m, k, l);
                    v2[R++] = m3_f_k1_l2(m, k, l);
                }
                m=k;
                v[k3_index(m)] = m3_f_m2_l1(m, l);
                v2[k3_index(m)] = m3_f_m1_l2(m, l);
                
                m=l;
                v[k3_index(m)] = m3_f_m2_l1(m, k);
                v2[k3_index(m)] = m3_f_m1_l2(m, k);
                
                R=k2l_fence;
                for (m=0; m<k; m++) {
                    for (n=m+1; n<Nmodes3; n++) {
                        v[R] = m2_n1_f_k2_l1(m, n, k, l);
                        v2[R++] = m2_n1_f_k1_l2(m, n, k, l);
                        
                        v[R] = m1_n2_f_k2_l1(m, n, k, l);
                        v2[R++] = m1_n2_f_k1_l2(m, n, k, l);
                    }
                    
                    n=k;
                    R = k2l_index(m, n);
                    v[R] = m2_n1_f_m2_l1(n, m, l);
                    v2[R++] = m2_n1_f_m1_l2(n, m, l);
                    
                    
                    v[R] = m1_n2_f_m2_l1(n, m, l);
                    v2[R] = m1_n2_f_m1_l2(n, m, l);
                    
                    n=l;
                    R = k2l_index(m, n);
                    v[R] = m2_n1_f_m2_l1(n, m, k);
                    v2[R++] = m2_n1_f_m1_l2(n, m, k);
                    
                    
                    v[R] = m1_n2_f_m2_l1(n, m, k);
                    v2[R] = m1_n2_f_m1_l2(n, m, k);
                }
                
                
                m=k;
                for (n=m+1; n<Nmodes3; n++) {
                    R = k2l_index(m, n);
                    v[R] = m2_n1_f_m2_l1(m, n, l);
                    v2[R++] = m2_n1_f_m1_l2(m, n, l);
                    
                    
                    v[R] = m1_n2_f_m2_l1(m, n, l);
                    v2[R] = m1_n2_f_m1_l2(m, n, l);
                }
                
                n=l;
                R = k2l_index(m, n);
                v[R] = m2_n1_f_m2_n1(m, n);
                v2[R++] = m2_n1_f_m1_n2(m, n);
                
                v2[R] = m1_n2_f_m1_n2(m, n);
            }
        }
        
        for (k=0; k<Nmodes3; k++) {
            for (l=k+1; l<Nmodes3; l++) {
                for (j=l+1; j<Nmodes3; j++) {
                    int C = klj_index(k, l, j);
                    int R=0;
                    
                    v.fill(0.0);
                    
                    v[0] = f_k1_l1_j1(k, l, j);
                    
                    for (m=0; m<Nmodes; m++) {
                        v[m+1] = m1_f_k1_l1_j1(m, k, l, j);
                    }
                    v[k+1] = m1_f_m1_l1_j1(k, l, j);
                    v[l+1] = m1_f_m1_l1_j1(l, k, j);
                    v[j+1] = m1_f_m1_l1_j1(j, k, l);
                    
                    R = k2_fence;
                    for (m=0; m<Nmodes2; m++) {
                        v[R++] = m2_f_k1_l1_j1(m, k, l, j);
                    }
                    v[k2_index(k)] = m2_f_m1_l1_j1(k, l, j);
                    v[k2_index(l)] = m2_f_m1_l1_j1(l, k, j);
                    v[k2_index(j)] = m2_f_m1_l1_j1(j, k, l);
                    
                    for (m=0; m<Nmodes2; m++) {
                        R = kl_fence;
                        for (n=m+1; n<Nmodes2; n++) {
                            v[R++] = m1_n1_f_k1_l1_j1(m, n, k, l, j);
                        }
                    }
                    
                    m=k;
                    for (n=m+1; n<Nmodes2; n++) {
                        v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, l, j);
                    }
                    n=l;
                    v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(m, n, j);
                    n=j;
                    v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(m, n, l);
                    
                    
                    m=l;
                    for (n=m+1; n<Nmodes2; n++) {
                        v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, k, j);
                    }
                    n=j;
                    v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(m, n, k);
                    
                    m=j;
                    for (n=m+1; n<Nmodes2; n++) {
                        v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(m, n, k, l);
                    }
                    
                    R = k2l_fence;
                    for (m=0; m<Nmodes3; m++) {
                        for (n = m+1; n<Nmodes3; n++) {
                            v[R++] = m2_n1_f_k1_l1_j1(m, n, k, l, j);
                            v[R++] = m2_n1_f_k1_l1_j1(n, m, k, l, j);
                        }
                    }
                    
                    m=k;
                    for (n=m+1; n<k; n++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, l, j);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, l, j);
                    }
                    m=l;
                    for (n=m+1; n<Nmodes3; n++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, k, j);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, k, j);
                    }
                    m=j;
                    for (n=m+1; n<Nmodes3; n++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, k, l);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, k, l);
                    }
                    
                    n=k;
                    for (m=0; m<n; m++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, l, j);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, l, j);
                    }
                    n=l;
                    for (m=0; m<n; m++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, k, j);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, k, j);
                    }
                    n=j;
                    for (m=0; m<n; m++) {
                        R = k2l_index(m, n);
                        v[R++] = m2_n1_f_m1_l1_j1(m, n, k, l);
                        v[R++] = m2_n1_f_m1_l1_j1(n, m, k, l);
                    }
                    
                    m=k; n=l;
                    R = k2l_index(m, n);
                    v[R++] = m2_n1_f_m1_n1_j1(m, n, j);
                    v[R  ] = m2_n1_f_m1_n1_j1(n, m, j);
                    
                    n=j;
                    R = k2l_index(m, n);
                    v[R++] = m2_n1_f_m1_n1_j1(m, n, l);
                    v[R  ] = m2_n1_f_m1_n1_j1(n, m, l);
                    
                    m=l; n=j;
                    R = k2l_index(m, n);
                    v[R++] = m2_n1_f_m1_n1_j1(m, n, k);
                    v[R  ] = m2_n1_f_m1_n1_j1(n, m, k);
                    
                    R = k3_index(0);
                    for (m=0; m<Nmodes3; n++) {
                        v[R++] = m3_f_k1_l1_j1(m, k, l, j);
                    }
                    m=k;
                    v[k3_index(m)] = m3_f_m1_l1_j1(m, l, j);
                    m=l;
                    v[k3_index(m)] = m3_f_m1_l1_j1(m, k, j);
                    m=j;
                    v[k3_index(m)] = m3_f_m1_l1_j1(m, k, l);

                    
                    R = klj_index(0, 1, 2);
                    for (m=0; m<k; m++) {
                        for (n=m+1; n<Nmodes3; n++) {
                            for (p=n+1; p<Nmodes3; p++) {
                                
                                v[R++] = m1_n1_p1_f_k1_l1_j1(m, n, p, k, l, j);
                            }
                        }
                        
                        /* m<k, n=k */
                        n=k;
                        for (p=n+1; p<Nmodes3; p++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, l, j);
                        }
                        n=l;
                        for (p=n+1; p<Nmodes3; p++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, k, j);
                        }
                        n=j;
                        for (p=n+1; p<Nmodes3; p++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(n, m, p, k, l);
                        }
                        
                        p=k;
                        for (n=m+1; n<p; n++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, l, j);
                        }
                        p=l;
                        for (n=m+1; n<p; n++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, k, j);
                        }
                        p=j;
                        for (n=m+1; n<p; n++) {
                            v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(p, m, n, k, l);
                        }
                        
                        n=k; p=l;
                        v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, j);
                        p=j;
                        v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, l);

                        n=l; p=j;
                        v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(n, p, m, k);
                    }
                    
                    m=k;
                    for (n=m+1; n<l; n++) {
                        R = klj_index(m, n, n+1);
                        for (p=n+1; p<Nmodes3; p++) {
                            v[R++] = m1_n1_p1_f_m1_l1_j1(m, n, p, l, j);
                        }
                        
                        p=l;
                        v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(m, p, n, j);
                        p=j;
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
                    
                    M.col(C).subvec(0, C-1) += v.subvec(0, C-1);
                }
            }
        }
    }
    
    void test_index()
    {
        cout << klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1) << " ";
        cout << (6*Nmodes + 3*Nmodes2*(1 + Nmodes2) + (3 + Nmodes3)*(2 + Nmodes3*Nmodes3))/6;
        cout << endl;
    }

    void addHODiagonal(mat &M)
    {
        double E0 = 0.5*sum(omega);
        
        M(0,0) += E0;
        for (int i=0; i<Nmodes; i++) {
            M(i+1,i+1) += E0+omega(i);
        }
        
        for (int i=0; i<Nmodes2; i++) {
            M(k2_index(i),k2_index(i)) += E0+2*omega(i);
        }
        
        for (int k=0; k<Nmodes2; k++) {
            for (int l=k+1; l<Nmodes2; l++) {
                int klpos = kl_index(k, l);
                
                M(klpos,klpos) += E0 + omega(k) + omega(l);
            }
        }
    }
    
    
    mat transitionDipole(vec& charges, mat& MUa, mat& C)
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
        
        mat mu(3,C.n_cols-1);
        for (int f=1; f<C.n_cols; f++) {
            mat d(3,3);
            d.col(0) = mu_0_m* C(0,0)*C(span(1,Nmodes),f);
            d.col(1) = mu_0_m* C(span(1,Nmodes), 0)*C(0,f);
            
            if (C.n_cols == Nmodes + 1) continue;
            
            //vec d1 = mu_m_m2* ( C(span(1,Nmodes),0) % C(span(Nmodes+1,Nmodes+Nmodes2), f)
            //                       + C(span(Nmodes+1,Nmodes+Nmodes2),0) % C(span(1,Nmodes),f) );
            d.col(2).fill(0.0);
            for (int m=0; m<Nmodes2-1; m++) {
                span s(kl_index(m,m+1), kl_index(m,Nmodes2-1));
                d.col(2) += mu_0_m.cols(m+1,Nmodes2-1) * (C(m, 0)*C(s,f) + C(s,0)*C(m,f));
            }
            //cout << reshape(d, 1, 9);
            mu.col(f-1) = d.col(0)+d.col(2);// + d1 + d2;
        }
        
        return mu;
    }
};

#endif
