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

#include "SimplexIterator.h"

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
        int kth_plane_offset = (k*(2 + k*k - 3*k*(Nmodes3 - 1) + 3*(Nmodes3 - 2)*Nmodes3))/6;
        /*
        int Np = Nmodes3 - k - 1;
        int lp = l - k - 1;
        int jp = j - k - 1;
        return klj_fence + jl_offset + lp*Np - ((lp + 2) * (lp + 1))/2 + jp;
         */
        int jl_pos = j + (k*k - l*(3 + l) + k*(3 - 2*Nmodes3) + 2*(l - 1)*Nmodes3)/2;
        return klj_fence + kth_plane_offset + jl_pos;
    }


    void addEpotSingles(mat &M)
    {
        int m, k;

        vec v(M.n_rows);

#pragma omp for schedule(dynamic)
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

    void addEpotDoubles(mat &M)
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
        for (int i=0; i<s2.end(); i++) {
            s2 = i;

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
        }
    }

    void addEpotTriples(mat &M)
    {
        int m, n, p, k, l, j;

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
        for (int i=0; i<s2.end(); i++) {
            s2 = i;

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
            v2[k_index(m)] = m1_f_m1_l2(m, l); // unique, no coordinate swap possible!

            m=l;
            v[k_index(m)] = m1_f_m2_l1(m, k);
            v2[k_index(m)] = m1_f_m1_l2(m, k); // unique again.

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
        }

        SimplexIterator<3> s3(Nmodes3);
#pragma omp for schedule(dynamic,2)
        for (int i=0; i<s3.end(); i++) {
            s3 = i;

            k = s3.index[0];
            l = s3.index[1];
            j = s3.index[2];

            int C = klj_index(k, l, j);
            int R=0;

            v.subvec(0, C).fill(0.0);

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

            R = kl_fence;
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
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(k, m, l, j);
            }

            n=l;
            for (m=0; m<n; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(l, m, k, j);
            }

            n=j;
            for (m=0; m<n; m++) {
                v[kl_index(m, n)] = m1_n1_f_m1_l1_j1(j, m, k, l);
            }

            m=k; n=l;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(k, l, j);
            n=j;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(k, j, l);


            m=l; n=j;
            v[kl_index(m, n)] = m1_n1_f_m1_n1_j1(l, j, k);


            R = k3_fence;
            for (m=0; m<Nmodes3; m++) {
                v[R++] = m3_f_k1_l1_j1(m, k, l, j);
            }
            v[k3_index(k)] = m3_f_m1_l1_j1(k, l, j);
            v[k3_index(l)] = m3_f_m1_l1_j1(l, k, j);
            v[k3_index(j)] = m3_f_m1_l1_j1(j, k, l);


            R = k2l_fence;
            for (m=0; m<Nmodes3; m++) {
                for (n = m+1; n<Nmodes3; n++) {
                    v[R++] = m2_n1_f_k1_l1_j1(m, n, k, l, j);
                    v[R++] = m1_n2_f_k1_l1_j1(m, n, k, l, j);
                }
            }

            m=k;
            for (n=m+1; n<k; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(k, n, l, j);
                v[R++] = m1_n2_f_m1_l1_j1(k, n, l, j);
            }

            m=l;
            for (n=m+1; n<Nmodes3; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(l, n, k, j);
                v[R++] = m1_n2_f_m1_l1_j1(l, n, k, j);
            }

            m=j;
            for (n=m+1; n<Nmodes3; n++) {
                R = k2l_index(m, n);
                v[R++] = m2_n1_f_m1_l1_j1(j, n, k, l);
                v[R++] = m1_n2_f_m1_l1_j1(j, n, k, l);
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
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(k, m, p, l, j);
                }
            }
            n=l;
            for (m=0; m<n; m++) {
                for (p=n+1; p<Nmodes3; p++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(l, m, p, k, j);
                }
            }
            n=j;
            for (m=0; m<n; m++) {
                for (p=n+1; p<Nmodes3; p++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(j, m, p, k, l);
                }
            }

            p=k;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(k, m, n, l, j);
                }
            }
            p=l;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(l, m, n, k, j);
                }
            }
            p=j;
            for (m=0; m<k; m++) {
                for (n=m+1; n<p; n++) {
                    v[klj_index(m, n, p)] = m1_n1_p1_f_m1_l1_j1(j, m, n, k, l);
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
                    v[R++] = m1_n1_p1_f_m1_l1_j1(k, n, p, l, j);
                }
            }

            p=l;
            for (n=m+1; n<p; n++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(k, l, n, j);
            }
            p=j;
            for (n=m+1; n<p; n++) {
                v[klj_index(m, n, p)] = m1_n1_p1_f_m1_n1_j1(k, j, n, l);
            }

            n=l;
            R = klj_index(m, n, n+1);
            for (p=n+1; p<j; p++) {
                v[R++] = m1_n1_p1_f_m1_n1_j1(k, l, p, j);
            }

            p=j;
            v[R] = m1_n1_p1_f_m1_n1_p1(k, l, j);

            if ( R!= C) {
                cerr << "R != C\n";
                exit(EXIT_FAILURE);
            }

            M.col(C).subvec(0, C) += v.subvec(0, C);
        }
    }

public:
    void addEpot(vec &q_, double V_, vec& Vq_, mat &M)
    {
        q = q_.memptr();
        V = V_;
        Vq = Vq_.memptr();

        M(0,0) += V;

#pragma omp parallel
        {
            addEpotSingles(M);
            addEpotDoubles(M);
            addEpotTriples(M);
        }
    }
    void test_index()
    {
        cout << klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1) << " ";
        cout << 1 + Nmodes + Nmodes2*(Nmodes2 + 1)/2 + (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6;
        cout << endl;
    }

    void addHODiagonal(mat &M)
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

        if (C-1 != klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1) ) {
            cerr << "HO diagonal term ended in wrong col number." << endl;
            cerr << C-1 << " != " << klj_index(Nmodes3-3, Nmodes3-2, Nmodes3-1);
            cerr << endl;
            exit(EXIT_FAILURE);
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
