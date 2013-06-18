/*
 *  HarmonicOscillatorBasis.cpp
 *  SCP1d
 *
 *  Created by Ionuț Georgescu on 4/30/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <cmath>

using namespace std;

#include "HarmonicOscillatorBasis.h"

#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Sqrt(x)         (sqrt((double)(x)))

static const double    pi_one_4th = 1.0 / sqrt(sqrt(M_PI));

static double 
ho_0(double y)
{
	return 1.0;
}


static double 
ho_1(double y)
{
	return 2.0 * y / sqrt(2.0);
}


static double 
ho_2(double y)
{
	return (-1.0 + 2.0 * y*y) / sqrt(2.0);
}


static double 
ho_3(double y)
{
	return (y * (-3.0 + 2.0 * y*y) ) / sqrt(3.0);
}


basis_function  ho_basis[] = {&ho_0, &ho_1, &ho_2, &ho_3};

int ho_basis_size = sizeof(ho_basis) / sizeof(basis_function);