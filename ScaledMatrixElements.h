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


#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

static double sqr(double x)
{
    return x*x;
}

class ScaledMatrixElements {
protected:
    double V;
    int Nmodes, Nmodes2, Nmodes3;
    vec omega, Vq, q;

    int k_fence, k2_fence, kl_fence, k2l_fence, k3_fence, klj_fence;

//#include "DoubleExcitations.h"
#include "TripleExcitations.h"

    int k_index(int k);
    int k2_index(int k);

    /**
     1 + 2*Nmodes + (Nmodes - 1 + Nmodes - k)*k/2 + l - k - 1
     */
    int kl_index(int k, int l);
    int k3_index(int k);
    int k2l_index(int k,int l);

    /**
     \f[
     \sum_{n=Nmodes3-(k-1)-1}^{Nmodes3-1} n*(n-1)/2 +l*(Nmodes3-k-1) - ((l + 2) * (l + 1))/2 + j;
     */
    int klj_index(int k, int l, int j);


    void addEpotSingles(mat &M);
    void addEpotDoubles(mat &M);
    void addEpotTriples(mat &M);
    
    vec get_bra(const vec& q);
    
public:
    ScaledMatrixElements(const vec &omega_, int Nmodes2_ = 0, int Nmodes3_ = 0) :
    omega(omega_), Nmodes2(Nmodes2_), Nmodes3(Nmodes3_)
    {
        Nmodes = omega.n_rows;
        k_fence = 1;
        k2_fence = k_fence + Nmodes;
        kl_fence = k2_fence + Nmodes2;
        k3_fence = kl_fence + (Nmodes2*(Nmodes2-1)) / 2;
        k2l_fence = k3_fence + Nmodes3;
        klj_fence = k2l_fence + Nmodes3*(Nmodes3-1);
    }
    
    size_t getBasisSize();
    size_t getSubBasisSize(int n);

    void addEpot(const vec &q_, double V_, const vec& Vq_, mat &M);
    void test_index();
    void addHODiagonal(mat &M);
    mat transitionDipole(vec& charges, mat& MUa, mat& C);
};

#endif
