/*
 *  nm_vmd.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 1/30/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

static const double bohr=0.52917721092;
static const double autocm=2.194746313e5;
static const double qM = 1.1128;

unsigned int NSobol, continue_skip=0, Nmodes=0;
string input_file, continue_from;

     
 double NuclearMass(string &species)
{
    static const double Hmass=1837.15137;
    static const double Cmass=21891.6543;
    
    
    if (species.compare("H") == 0 ) {
        return Hmass;
    }
    else if (species.compare("C") == 0 ) {
        return Cmass;
    }
    else if (species.compare("O") == 0 ) {
        return Cmass*15.9994/12.0107;
    }
    else {
        return -1.0;
    }
    
}
     
     
 void
 load_from_vladimir(string &name, int &N, vec &mass, vec& x0, mat& H)
{
    ifstream fin(name.c_str());
    
    fin >> N;
    
    mass.resize(3*N);
    x0.resize(3*N);
    H.resize(3*N, 3*N);
    
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
    
    for (int i=0; i < 3*N; i++) {
        for (int j=0; j <= i; j++) {
            fin >> H(j, i);
            H(i, j) = H(j, i);
        }
    }
    
    fin.close();
}
 
 extern "C" {
     void dsyevr_(char *, char *, char *,  int *,  double *,  int *,  double *,
                  double* ,  int *IL,  int *IU, double *ABSTOL,  int *M,
                  double *W, double *Z, int *LDZ, int *ISUPPZ, double *WORK,
                  int *LWORK, int *IWORK, int *LIWORK, int *INFO);
     double dlamch_(char *);
 }
 
 
 
 vec dsyevr(mat &A, mat *C=NULL)
{
    
    char JOBZ='N';
    
    char RANGE='A';
    char UPLO='U';
    int N = A.n_rows;
    int LDA = A.n_rows;
    double VU = 0.0, VL = 0.0;
    int IL = 0, IU = 0;
    double ABSTOL = dlamch_("Safe minimum");
    int M;
    vec W(A.n_rows);
    double *Z = NULL;
    int LDZ = LDA;
    Col<int> ISUPPZ(2*N);
    int LWORK = 36*N;
    vec WORK(LWORK);
    int LIWORK = 10*N;
    Col<int> IWORK(LIWORK);
    int INFO;
    
    if (C != NULL) {
        JOBZ='V';
        C->set_size(A.n_rows, A.n_cols);
        Z = C->memptr();
    }
    dsyevr_(&JOBZ, &RANGE, &UPLO,  &N,  A.memptr(),  &LDA,  &VL,  &VU,
            &IL,  &IU, &ABSTOL, &M,  W.memptr(),  Z, &LDZ, ISUPPZ.memptr(),
            WORK.memptr(), &LWORK, IWORK.memptr(), &LIWORK, &INFO );
    
    return W;
}
 
 
 
 int main (int argc, char *  argv[]) {
     int N;
     mat H, Up, U;
     vec mass, x0, omegasq0;

     string fin(argv[1]);
     load_from_vladimir(fin, N, mass, x0, H);
     omegasq0 = dsyevr(H, &Up);
     
     U = Up.cols(6, Up.n_cols - 1);
     float nm_arrow_len = 1.0;
     float nm_arrow_radius = 0.1;
     float nm_arrow_tip_radius = 2.0 * nm_arrow_radius;
     float nm_arrow_tip_len = 1.5 * nm_arrow_tip_radius;
     
     x0 *= bohr;
     ofstream fout("nm.vmd");
     for (int i=0; i<U.n_cols; i++) {
         fout << "display projection orthographic\n"
             << "mol new \"" + fin + "\" type xyz\n"
             << "mol rename " << i << " \"mode " << i << "\"\n"
             << "mol off " << i << endl;
         
         for (int j=0; j < U.n_rows / 3; j++) {
             vec x1 = x0(span(3*j, 3*j+2));
             vec x2 = x1 + nm_arrow_len * U(span(3*j, 3*j+2), i);
             vec x3 = x2 + nm_arrow_tip_len * U(span(3*j, 3*j+2), i);
             
             fout<<"draw cylinder ";
             fout<<" {"<< x1(0) <<" "<< x1(1) <<" "<< x1(2) <<"}";
             fout<<" {"<< x2(0) <<" "<< x2(1) <<" "<< x2(2) <<"}";
             fout<<" radius " << nm_arrow_radius <<" filled yes" << endl;
             
             fout<<"draw cone ";
             fout<<" {"<< x2(0) <<" "<< x2(1) <<" "<< x2(2) <<"}";
             fout<<" {"<< x3(0) <<" "<< x3(1) <<" "<< x3(2) <<"}";
             fout<<" radius " << 2*nm_arrow_radius << endl;
         }
     }
     fout.close();
 }
