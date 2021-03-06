/*
 *  DiskIO.h
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/2/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#ifndef __DISKIO_H__
#define __DISKIO_H__

#include <string>

#include <armadillo>

using namespace std;
using namespace arma;

extern void
load_from_vladimir(string &name, vec &mass, vec& x0, mat& H);

extern void
save_for_vladimir(const string &name, double F, vec& x0, mat& H);

extern uvec
OHHOHH(vec& mass, vec& r);

extern uvec
OHHOHH(vec& mass, vec& r, mat& H);

extern void
load_xyz(string &name, vec &mass, vec& x0);

extern void
save_hdf5(mat &M, char *name);

void save_hdf5(fcube &C, char *name);

extern void
load_hdf5(string& name, mat &M);

#endif
