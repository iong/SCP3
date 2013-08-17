/*
 *  DiskIO.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/2/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>

#include <algorithm>

#include <H5Cpp.h>

#include "Constants.h"
#include "DiskIO.h"

using namespace H5;

static double NuclearMass(string &species)
{ 
    if (species.compare("H") == 0 ) {
        return Hmass;
    }
    else if (species.compare("D") == 0 ) {
        return Dmass;
    }
    else if (species.compare("C") == 0 ) {
        return Cmass;
    }
    else if (species.compare("O") == 0 ) {
        return Omass;
    }
    else {
        return -1.0;
    }
    
}

void
load_xyz(string &name, vec &mass, vec& x0)
{
    int N;
    

    ifstream fin(name.c_str());
    
    fin >> N;
    
    mass.resize(3*N);
    x0.resize(3*N);
    
    char skipbuf[256];
    fin.getline(skipbuf, sizeof(skipbuf));
    fin.getline(skipbuf, sizeof(skipbuf));
    
    for (int i=0; i<N; i++) {
        string species;
        
        fin >> species;
        mass( span(3*i, 3*i+2) ).fill(NuclearMass(species));
        fin >> x0[3*i] >> x0[3*i+1] >> x0[3*i+2];
    }

    fin.close();
}



void
load_from_vladimir(string &name, vec &mass, vec& x0, mat& H)
{
    int N;
    ifstream fin(name.c_str());
    
    fin >> N;
    
    mass.resize(3*N);
    x0.resize(3*N);
    
    char skipbuf[256];
    fin.getline(skipbuf, sizeof(skipbuf));
    fin.getline(skipbuf, sizeof(skipbuf));
    
    for (int i=0; i<N; i++) {
        string species;
        
        fin >> species;
        mass( span(3*i, 3*i+2) ).fill(NuclearMass(species));
        fin >> x0[3*i] >> x0[3*i+1] >> x0[3*i+2];
    }
    x0 /= bohr;

    H.set_size(3*N, 3*N);
    
    for (int i=0; i < 3*N; i++) {
        for (int j=0; j <= i; j++) {
            fin >> H(j, i);
            if (!fin.good())
            {
                H.reset();
                goto done;
            }
            H(i, j) = H(j, i);
        }
    }
done: 
    fin.close();
}

void
save_for_vladimir(const string &name, double F, vec& x0, mat& H)
{
    int N = x0.n_rows / 3;
    ofstream fout(name.c_str());
    
    fout << N << endl
        << "Free energy = " << F << endl;
    
    for (int i=0; i<x0.n_rows; i += 9) {
        fout << "O " << x0.subvec(i  , i+2).t()*bohr;
        fout << "H " << x0.subvec(i+3, i+5).t()*bohr;
        fout << "H " << x0.subvec(i+6, i+8).t()*bohr;
    }
    
    for (int i=0; i < H.n_cols; i++) {
        for (int j=0; j <= i; j++) {
            fout << H(j, i) << endl;
        }
    }
    
    fout.close();
}

uvec OHHOHH(vec& mass, vec& r)
{
    int iO = 0;
    int iH = 3;
    
    uvec p(mass.n_rows);
    
    for (int i=0; i<mass.n_rows;) {
        if (mass[i] == Omass) {
            for (int j=0; j<3; j++) {
                p[iO + j] = i + j;
            }
            
            iO += 9;
            i += 3;
        }
        else if (mass[i] == Hmass || mass[i] == Dmass) {
            for (int j=0; j<3; j++) {
                p[iH + j] = i + j;
            }
            
            iH +=  3;
            iH += 3 * (iH%9 == 0);
            i += 3;
        }
        else {
            cerr << "Uknown mass at i = "<< i << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    vec r_tmp(r);
    for (int i=0; i<p.n_rows; i++) {
        r[i] = r_tmp[p[i]];
    }

    vec mass2 = mass(p);
    mass = mass2;
    
    return p;
}


uvec OHHOHH(vec& mass, vec& r, mat& H)
{  
    uvec p = OHHOHH(mass, r);
    
    if (!H.is_empty()) {
        mat H_tmp = H.cols(p);
        H = H_tmp.rows(p);
    }
    
    return p;
}




/** Save in column major mode **/
void save_hdf5(mat &M, char *name)
{
    H5File fout(name, H5F_ACC_TRUNC);
    
    hsize_t dims[2];
    dims[0] = M.n_cols;
    dims[1] = M.n_rows;
    DataSpace dspace(2, dims);
    
    DSetCreatPropList plist;
    
    hsize_t chunk_size[2];
    chunk_size[0] = min(dims[0], (128*1024 + dims[1] - 1)/ dims[1]);
    chunk_size[1] = dims[1];
    
    plist.setChunk(2, chunk_size);
    plist.setDeflate(6);
    
    FloatType dtype(PredType::NATIVE_FLOAT);
    DataSet dset = fout.createDataSet("M", dtype, dspace, plist);
    
    FloatType mem_dtype(PredType::NATIVE_DOUBLE);
    dset.write(M.memptr(), mem_dtype);
    
    dset.close();
    fout.close();
}

/** Save in column major mode **/

void save_hdf5(fcube &C, char *name)
{
    H5File fout(name, H5F_ACC_TRUNC);
    
    hsize_t dims[3];
    dims[0] = C.n_slices;
    dims[1] = C.n_cols;
    dims[2] = C.n_rows;
    DataSpace dspace(3, dims);
    
    DSetCreatPropList plist;
    
    hsize_t chunk_size[3];
    chunk_size[0] = (128*1024 + dims[1]*dims[2] - 1)/ (dims[1]*dims[2]);
    chunk_size[1] = dims[1];
    chunk_size[2] = dims[2];
    
    chunk_size[0] = min((hsize_t)C.n_slices, max((hsize_t)1, chunk_size[0]));
    
    if (chunk_size[0] < 1) {
        chunk_size[0] = 1;
    }
    
    plist.setChunk(3, chunk_size);
    plist.setDeflate(6);
    
    FloatType dtype(PredType::NATIVE_FLOAT);
    DataSet dset = fout.createDataSet("M", dtype, dspace, plist);
    
    FloatType mem_dtype(PredType::NATIVE_FLOAT);
    dset.write(C.memptr(), mem_dtype);
    
    dset.close();
    fout.close();
}

void load_hdf5(string& name, mat &M)
{
    H5File fin(name, H5F_ACC_RDONLY);
    DataSet dset = fin.openDataSet("M");
    
    DataSpace file_space = dset.getSpace();
    
    hsize_t dims[2];
    int ndims = file_space.getSimpleExtentDims(dims);
    
    if (ndims != 2 || dims[0] != dims[1]) {
        cerr << "HDF5 dataset \"M\" has wrong shape. ndims = " << ndims
        << ", dims = "<<dims[0]<<"x"<<dims[1]<<endl;
    }
    
    M.set_size(dims[1], dims[0]);
    
    FloatType mem_dtype(PredType::NATIVE_DOUBLE);
    dset.read(M.memptr(), mem_dtype);
    
    dset.close();
    fin.close();
}


