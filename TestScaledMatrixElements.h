//
//  TestScaledMatrixElements.h
//  SCP_Double_Excitations
//
//  Created by Ionut Georgescu on 3/18/13.
//
//

#ifndef __SCP_Double_Excitations__TestScaledMatrixElements__
#define __SCP_Double_Excitations__TestScaledMatrixElements__

#include <iostream>
#include "ScaledMatrixElements.h"

class TestScaledMatrixElements : public ScaledMatrixElements {
public:
    double t0;
    
    TestScaledMatrixElements(vec &omega_, int Nmodes2_ = 0, int Nmodes3_ = 0) :
        ScaledMatrixElements(omega_, Nmodes2_, Nmodes3_) {}
    
    void addEpot(vec &q_, double V_, vec& Vq_, vec &v);
};

#endif /* defined(__SCP_Double_Excitations__TestScaledMatrixElements__) */
