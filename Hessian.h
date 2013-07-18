//
//  Hessian.h
//  SCP3
//
//  Created by Ionut Georgescu on 7/15/13.
//
//

#ifndef __SCP3__Hessian__
#define __SCP3__Hessian__

#include <iostream>

#include <armadillo>

using namespace arma;

#include "PES.h"

extern void getHessian(h2o::PES& pot, const vec& r_au, double s, mat& H);
extern void massScaleHessian(const vec& mass, mat& H);
extern mat transrotBasis(const vec& q);
extern mat transrotBasis(const vec& q, const vec& mass);
extern void regtransrot(const mat& U, mat &H);

#endif /* defined(__SCP3__Hessian__) */
