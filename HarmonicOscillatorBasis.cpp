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

static double 
ho_0(double y)
{
	return 1.0;
}


static double 
ho_1(double y)
{
	return sqrt(2.0) * y;
}


static double 
ho_2(double y)
{
	return (2.0 * y*y - 1.0) / sqrt(2.0);
}


static double 
ho_3(double y)
{
	return (2.0 * y*y - 3.0) * y / sqrt(3.0);
}


basis_function  ho_basis[] = {&ho_0, &ho_1, &ho_2, &ho_3};

int ho_basis_size = sizeof(ho_basis) / sizeof(basis_function);