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
load_from_vladimir(string &name, int &N, vec &mass, vec& x0, mat& H);

extern void
load_xyz(string &name, vec &mass, vec& x0);

extern void
save_hdf5(mat &M, char *name);

extern void
load_hdf5(string& name, mat &M);

#endif
