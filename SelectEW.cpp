/*
 *  SelectEW.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/2/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */
#include <string>
#include <H5Cpp.h>

#include "DiskIO.h"
#include "F90.h"

using namespace std;
using namespace H5;

vec dsyevr(mat &A, int IU=0, mat *C=NULL)
{
    
    char JOBZ='N';
    
    char RANGE='A';
    char UPLO='U';
    int N = A.n_rows;
    int LDA = A.n_rows;
    double VU = 0.0, VL = 0.0;
    int IL = 0;
    double ABSTOL = dlamch_("Safe minimum");
    int M;
    vec W;
    double *Z = NULL;
    int LDZ = LDA;
    Col<int> ISUPPZ(2*N);
    int LWORK = 36*N;
    vec WORK(LWORK);
    int LIWORK = 10*N;
    Col<int> IWORK(LIWORK);
    int INFO;

    if (IU > 0) {
        IL = 1;
        RANGE = 'I';
        
        W.set_size(IU - IL + 1);
    }
    else {
        W.set_size(A.n_cols);
    }

    if (C != NULL) {
        JOBZ='V';
        if (IU > 0) {
            C->set_size(A.n_rows, W.n_rows);
        }
        else {
            C->set_size(A.n_rows, A.n_cols);
        }
        Z = C->memptr();
    }
    
    dsyevr_(&JOBZ, &RANGE, &UPLO,  &N,  A.memptr(),  &LDA,  &VL,  &VU,
            &IL,  &IU, &ABSTOL, &M,  W.memptr(),  Z, &LDZ, ISUPPZ.memptr(),
            WORK.memptr(), &LWORK, IWORK.memptr(), &LIWORK, &INFO );
    
    if (M != W.n_rows) {
        cerr << "#EW returned is different from requested.\n"
             << "\trequested: "<<W.n_rows<<endl
             << "\treturned: "<< M << endl;
    }
    
    return W;
}


int main(int argc, char *argv[])
{
    mat M;

    string h5name_in(argv[1]);
    int subspace_size = atoi(argv[2]);

    load_hdf5(h5name_in, M);
    
    mat EV;
    vec EW = dsyevr(M, subspace_size, &EV);
    
    string h5name_out = h5name_in.substr(0, h5name_in.length() - 3) + "_EWV.h5";
    H5File h5out(h5name_out, H5F_ACC_TRUNC);
    
    FloatType dtype(PredType::NATIVE_FLOAT);
    FloatType mem_dtype(PredType::NATIVE_DOUBLE);

    hsize_t dims[2];
    dims[0] = EW.n_rows;
    dims[1] = EV.n_rows;
    
    DataSpace dspace(1, dims);
    DataSet dset = h5out.createDataSet("EW", dtype, dspace);
    dset.write(EW.memptr(), mem_dtype);
    dset.close();
    
    dspace = DataSpace(2, dims);
    dset = h5out.createDataSet("EV", dtype, dspace);
    dset.write(EV.memptr(), mem_dtype);
    dset.close();

    h5out.close();
}