/*
 *  HarmonicOscillatorBasis.h
 *  SCP1d
 *
 *  Created by Ionuț Georgescu on 4/30/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */



typedef double (*basis_function)(double x);

extern basis_function ho_basis[];
extern int ho_basis_size;