//
//  TestScaledMatrixElements.cpp
//  SCP_Double_Excitations
//
//  Created by Ionut Georgescu on 3/18/13.
//
//

#include "TestScaledMatrixElements.h"


void TestScaledMatrixElements::addEpot(vec &q_, double V_, vec& Vq_, vec &v)
{
    q = q_.memptr();
    V = V_;
    Vq = Vq_.memptr();
    
    if (v.n_rows != 32) {
        v.set_size(32);
        v.fill(0.0);
    }
    
    int R=0;
    
    t0 = m2_n1_f_m1_l1_j1( 0, 17, 10, 15);
    v[R++] += t0;
    v[R++] += m2_n1_f_m1_l1_j1(10, 17,  0, 15);
    v[R++] += m2_n1_f_m1_l1_j1(15, 17,  0, 10);
    
    v[R++] += m2_n1_f_m1_l1_j1( 0, 15, 10, 17);
    v[R++] += m2_n1_f_m1_l1_j1(10, 15,  0, 17);
    v[R++] += m2_n1_f_m1_l1_j1(17, 15,  0, 10);
    
    v[R++] += m2_n1_f_m1_l1_j1( 0, 10, 15, 17);
    v[R++] += m2_n1_f_m1_l1_j1(15, 10,  0, 17);
    v[R++] += m2_n1_f_m1_l1_j1(17, 10,  0, 15);
    
    v[R++] += m2_n1_f_m1_l1_j1(10, 0, 15, 17);
    v[R++] += m2_n1_f_m1_l1_j1(15, 0, 10, 17);
    v[R++] += m2_n1_f_m1_l1_j1(17, 0, 10, 15);
    
    v[R++] += m3_f_m3(0);
    v[R++] += m3_f_m3(10);
    v[R++] += m3_f_m3(15);
    v[R++] += m3_f_m3(17);
    
    v[R++] += m1_n1_f_m1_l1_j1( 0, 17, 10, 15);
    v[R++] += m1_n1_f_m1_l1_j1(10, 17,  0, 15);
    v[R++] += m1_n1_f_m1_l1_j1(15, 17,  0, 10);
    
    v[R++] += m1_n1_f_m1_l1_j1( 0, 15, 10, 17);
    v[R++] += m1_n1_f_m1_l1_j1(10, 15,  0, 17);
    v[R++] += m1_n1_f_m1_l1_j1(17, 15,  0, 10);
    
    v[R++] += m1_n1_f_m1_l1_j1( 0, 10, 15, 17);
    v[R++] += m1_n1_f_m1_l1_j1(15, 10,  0, 17);
    v[R++] += m1_n1_f_m1_l1_j1(17, 10,  0, 15);
    
    v[R++] += m1_n1_f_m1_l1_j1(10, 0, 15, 17);
    v[R++] += m1_n1_f_m1_l1_j1(15, 0, 10, 17);
    v[R++] += m1_n1_f_m1_l1_j1(17, 0, 10, 15);
    
    v[R++] += f_k1_l1(0, 10);
    v[R++] += f_k1_l1(0, 10);
    v[R++] += f_k1_l1(0, 15);
    v[R++] += f_k1_l1(0, 17);
}

