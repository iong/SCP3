/*
 *  MatrixElements.h
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 1/8/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */


#ifndef __MATRIX_ELEMENTS_H__
#define __MATRIX_ELEMENTS_H__

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

class MatrixElements {
    vec& omega;
    vec alpha, alphaq;
    double V, *Vq, *q;
    int Nmodes;
    
public:
    MatrixElements(vec &omega_) : omega(omega_)
    {
        alpha = sqrt(omega_);
        Nmodes = omega_.n_rows;
    }
    
private:
    double m_f_0(int m)
    {
        return Vq[m]/(sqrt(2.0) * alpha(m));
    }
    
    double m_f_k(int m, int k)
    {
        double res;
        res = 0.5 * (alphaq(m)*Vq[k]/alpha(k)
                     + alphaq(k)*Vq[m]/alpha(m));
        
        return res;
    }
    
    double m_f_m(int m)
    {
        return V + q[m]*Vq[m];
    }
    
    double m2_f_0(int m)
    {
        return Vq[m]*q[m]/sqrt(2.0);
    }
    
    double m2_f_k(int m, int k)
    {
        double res;
        
        res = (sqr(alphaq(m)) - 0.5) * Vq[k]/alpha(k);
        res += alpha(k) * q[m] * Vq[m];
        
        return 0.5*res;
    }
    
    double m2_f_m(int m)
    {
        return alphaq(m) * q[m] * Vq[m] + 0.5 * Vq[m]/alpha(m);
    }
    
    double m2_f_k2(int m, int k)
    {
        double res;
        
        res = (sqr(alphaq(m)) - 0.5) * q[k] * Vq[k] + (sqr(alphaq(k)) - 0.5) * q[m] * Vq[m];
        res *= 0.5;
        return res;
    }
    
    double m2_f_m2(int m)
    {
        return   (sqr(alphaq(m)) + 0.5) * q[m] * Vq[m] + V;
    }
    
    double mn_f_0(int m, int n)
    {
        return 0.5 * ( Vq[m]*alphaq(n)/alpha(m)
                      + Vq[n]*alphaq(m)/alpha(n) );
    }

    double mn_f_k(int m, int n, int k)
    {
        double res;
        
        res = sqrt(2.0) * ( alphaq(m)*alphaq(n)*Vq[k]/alpha(k) 
                           + alphaq(k)*alphaq(m)*Vq[n]/alpha(n)
                           + alphaq(n)*alphaq(k)*Vq[m]/alpha(m) );
        
        return res/3.0;
    }
    
    double mn_f_m(int m, int n)
    {
        double res;
        
        res = ( sqr(alphaq(m)) + 0.5) * Vq[n]/alpha(n) + alphaq(n)*q[m]*Vq[m];
        
        return res*sqrt(2.0)/2.0;
    }
    
    
    double mn_f_k2(int m, int n, int k)
    {
        double res;
        
        res = sqrt(2.0)*alphaq(m)*alphaq(n)*q[k]*Vq[k];
        res += (2*sqr(alphaq(k)) - 1.0) *  ( Vq[m]*alphaq(n)/alpha(m) 
                                            + Vq[n]*alphaq(m)/alpha(n) ) / sqrt(2.0);
        
        return res/3.0;
    }
    
    double mn_f_m2(int m, int n)
    {
        return alphaq(n)*Vq[m] * ( alphaq(m)*q[m] + 0.5/alpha(m) ) * sqrt(2.0);
    }
    
    double mn_f_kl(int m, int n, int k, int l)
    {
        double res;
        
        res = 0.5 * (alphaq(m)*alphaq(n)*alphaq(k) * Vq[l]/alpha(l)
                     + alphaq(l)*alphaq(m)*alphaq(n) * Vq[k]/alpha(k)
                     + alphaq(k)*alphaq(l)*alphaq(m) * Vq[n]/alpha(n)
                     + alphaq(n)*alphaq(k)*alphaq(l) * Vq[m]/alpha(m) );
        
        return res;
    }
    
    double mn_f_ml(int m, int n, int l)
    {
        double res;
        
        res = sqr(alphaq(m)) * ( alphaq(n) * Vq[l]/alpha(l)
                                + alphaq(l) * Vq[n]/alpha(n));
        
        return res;
    }
    
    double mn_f_mn(int m, int n)
    {
        double res;
        
        res = V + 0.5*(Vq[m] * q[m] + Vq[n] * q[n]);
        res += sqr(alphaq(m)) * q[n] * Vq[n];
        res += sqr(alphaq(n)) * q[m] * Vq[m];
        
        return res;
    }
    
    int double_index(int m)
    {
        return 1+Nmodes + m;
    }
    
    int single2_index(int m, int n)
    {
        if (m>n) {
            swap(m, n);
        }
        return 1 + 2*Nmodes + (Nmodes-1 + Nmodes - 1 - (m-1))*m/2 + n-m - 1;
    }

public:
    void addEpotSingles(vec &q_, double V_, vec& Vq_, mat &M)
    {
        q = q_.memptr();
        alphaq = alpha%q_;
        V = V_;
        Vq = Vq_.memptr();

        M(0,0) += V;
        
           
        for (int m=0; m<Nmodes; m++) {
            double v = m_f_0(m);
            M(0, 1+m) += v;
            M(1+m, 0) += v;
            M(1+m, 1+m) += m_f_m(m);
        }
        
        for (int m=0; m<Nmodes; m++) {
            for (int n=m+1; n<Nmodes; n++) {
                double v = m_f_k(m, n);
                M(1+m, 1+n) += v;
                M(1+n, 1+m) += v;
            }
        }
    }
    
    void addEpotDoubles(vec &q_, double V_, vec& Vq_, mat &M)
    {
        q = q_.memptr();
        alphaq = alpha%q_;
        V = V_;
        Vq = Vq_.memptr();

        M(0,0) += V;

        for (int m=0; m<Nmodes; m++) {
            M(0, 1+m) += m_f_0(m);
            M(0, double_index(m)) += m2_f_0(m);
            M(1+m, 1+m) += m_f_m(m);
            M(double_index(m), double_index(m)) += m2_f_m2(m);
            
            M(1+m, double_index(m)) += m2_f_m(m);
        }

        for (int m=0; m<Nmodes; m++) {
            for (int n=m+1; n<Nmodes; n++) {
                M(0,single2_index(m, n)) += mn_f_0(m, n);
                
                M(1+m, 1+n) += m_f_k(m, n);
                
                M(1+m, double_index(n)) += m2_f_k(n, m);
                M(1+n, double_index(m)) += m2_f_k(m, n);
                
                
                M(1+m, single2_index(m, n)) += mn_f_m(m, n);
                M(1+n, single2_index(m, n)) += mn_f_m(n, m);
                
                M(double_index(m), double_index(n)) += m2_f_k2(m, n);
                
                M(double_index(m), single2_index(m,n)) += mn_f_m2(m, n);
                M(double_index(n), single2_index(m,n)) += mn_f_m2(n, m);
                
                M(single2_index(m, n), single2_index(m, n)) += mn_f_mn(m, n);
            }
        }

        for (int m=0; m<Nmodes; m++) {
            for (int n=m+1; n<Nmodes; n++) {
                for (int k=n+1; k<Nmodes; k++) {
                    double v = mn_f_k(m, n, k);

                    M(k+1, single2_index(m, n)) += v;
                    M(n+1, single2_index(m, k)) += v;
                    M(m+1, single2_index(n, k)) += v;
                    
                    M(double_index(k), single2_index(m, n)) += mn_f_k2(m, n, k);
                    M(double_index(n), single2_index(m, k)) += mn_f_k2(k, m, n);
                    M(double_index(m), single2_index(n, k)) += mn_f_k2(n, k, m);
                    
                    M(single2_index(m, n), single2_index(m, k)) += mn_f_ml(m, n, k);
                    M(single2_index(m, n), single2_index(n, k)) += mn_f_ml(n, m, k);
                    M(single2_index(m, k), single2_index(n, k)) += mn_f_ml(k, m, n);
                }
            }
        }

#pragma omp parallel for schedule(dynamic)
        for (int m=0; m<Nmodes; m++) {
            for (int n=m+1; n<Nmodes; n++) {
                for (int k=n+1; k<Nmodes; k++) {
                    for (int l=k+1; l<Nmodes; l++) {
                        double v = mn_f_kl(m, n, k, l);

                        /**
                         * Symmetries:
                         * mn    <-> nm     : 2
                         * kl    <-> lk     : 2
                         * mn,kl <-> kl,mn  : 2
                         * unique permutations: 24 / (2*2*2) = 3
                         */
                        /*
                        if (k==(Nmodes-2)&&(l==Nmodes-1)) {
                            cout << m <<", " << n << " : " << single2_index(m, n);
                            cout << "," << single2_index(k, l) << " : " << v<< endl;
                        }
                         */
                        M(single2_index(m, n), single2_index(k, l)) += v;
                        M(single2_index(m, k), single2_index(n, l)) += v;
                        M(single2_index(m, l), single2_index(n, k)) += v;
                    }
                }
            }
        }
    }
    

    void addEkinSingles(mat& M)
    {
        double Ekin = sum(omega)/4.0;
        
        M(0,0) += Ekin;
        for (int i=0; i<Nmodes; i++) {
            M(i+1,i+1) += Ekin + 0.5*omega(i);
        }
    }
    
    void addEkinDoubles(mat& M)
    {
        double Ekin = sum(omega)/4.0;
        
        M(0,0) += Ekin;
        
        for (int i=0; i<Nmodes; i++) {
            M(i+1,i+1) += Ekin + 0.5*omega(i);
            M(0, double_index(i)) += - omega(i)*sqrt(2.0)/4.0;
            M(double_index(i), double_index(i)) += Ekin + omega(i);
            for (int j=i+1; j < Nmodes; j++) {
                M(single2_index(i, j), single2_index(i, j)) += Ekin + 0.5*(omega(i) + omega(j));
            }
        }
    }
    
    mat transitionDipole(vec& x0, vec& charges, mat& MU, mat& C)
    {
        vec mu0(3);
        
        mu0.fill(0.0);
        for (int i=0; i<x0.n_rows; i += 3) {
            mu0 += x0(span(i,i+2)) * charges[i/3];
        }
        mat Q(3,MU.n_rows);
        Q.fill(0);
        for (int i=0; i<Q.n_cols; i += 3) {
            Q(0, i  ) = charges[i/3];
            Q(1, i+1) = charges[i/3];
            Q(2, i+2) = charges[i/3];
        }
        
        mat mu_m_m2 = Q*MU;
        for (int i=0; i<3; i++) {
            mu_m_m2.row(i) /= alpha.t();
        }
        
        mat mu_0_m = mu_m_m2 / sqrt(2.0);
        
        mat mu(3,C.n_cols-1);
        for (int f=1; f<C.n_cols; f++) {
            mu.col(f-1) = mu_0_m* (C(0,0)*C(span(1,Nmodes),f) + C(span(1,Nmodes), 0)*C(0,f));
            
            if (C.n_cols == Nmodes + 1) continue;
            
            mu.col(f-1) += mu_m_m2* ( C(span(1,Nmodes),0) % C(span(Nmodes+1,2*Nmodes), f)
                                   + C(span(Nmodes+1,2*Nmodes),0) % C(span(1,Nmodes),f) );
            for (int m=0; m<Nmodes-1; m++) {
                span s(single2_index(m,m+1), single2_index(m,Nmodes-1));
                mu.col(f-1) += mu_0_m.cols(m+1,Nmodes-1) * (C(m, 0)*C(s,f) + C(s,0)*C(m,f));
            }
        }
        
        return mu;
    }
};

#endif