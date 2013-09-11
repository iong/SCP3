/*
 *  SelectEW.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/2/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */
#include <fstream>
#include <string>
#include <H5Cpp.h>

#include "Constants.h"
#include "DiskIO.h"
#include "F90.h"

using namespace std;
using namespace H5;

vec dsyevr(mat &A, int IU=0, mat *C=NULL)
{
    
    char JOBZ='N';
    
    char RANGE='A';

    char UPLO='U';
    if (A(0, A.n_cols-1) == 0.0) {
        UPLO='L';
    }

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

    string fin_name(argv[1]);
    int n_eigenpairs = atoi(argv[2]);

    load_hdf5(fin_name, M);
    
    mat EV;
    vec EW = dsyevr(M, n_eigenpairs, &EV);
    
    string fout_name = "eigenpairs" 
	    + fin_name.substr(2, fin_name.length() - 5)
	    + ".dat";

    ofstream fout(fout_name.c_str());
    fout.precision(15);

    fout <<"# Lowest "<< n_eigenpairs <<" eigenvalues\n";
    EW.raw_print(fout);

    fout << "# Corresponding eigenvectors, one per line." << endl;
    EV.t().raw_print(fout);

    fout.close();
}
