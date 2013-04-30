/*
 *  SCP3.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 2/5/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include <tr1/random>

#include <getopt.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"
#include "ScaledMatrixElements.h"

#include "sobol.hpp"

#include "qtip4pf.h"
#include "ttm3f.h"

#ifdef HAVE_BOWMAN
#include "bowman.h"
#endif



using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "doubles", required_argument, NULL, '2'},
    { "triples", required_argument, NULL, '3'},
    { "spectrum", required_argument, NULL, 's'},
    { "potential", no_argument, NULL, 'p'},
    {NULL, 0, NULL, 0}
};

static unsigned int  Nmodes2=0, Nmodes3=0;
static int64_t NSobol=1<<20, sobol_skip=1<<30;
static string input_file;
static string spectrum_file;

static string continue_from_file;
static unsigned int continue_from;
static ifstream rng_in;

string  h2o_potential("qtip4pf");

void dump_spectrum(vec &charges, mat &MU, mat &H, ScaledMatrixElements& me,
                   ostream& specout, ostream& dipoleout, ostream& E0out,
                   string Cout="")
{
    mat C;
    vec E;
    eig_sym(E, C, H);

    /*
    if (Cout.length() > 0) {
        C.save(Cout, raw_ascii);
    }
     */

    mat mu = me.transitionDipole(charges, MU, C);
    for (int ii=1; ii< E.n_rows; ii++) {
        specout << (E(ii) - E(0))*autocm << " ";
        dipoleout << norm(mu.col(ii-1), 2) << " ";
    }
    (specout  << endl).flush();
    (dipoleout << endl).flush();
}


void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "S:N:2:3:c:s:r:p:", program_options, NULL)) != -1) {
        int p2;
        switch (ch) {
            case 'N':
                p2 = atoi(optarg);
                if (p2>32) {
                    cerr << "-N <N>: <N>is too large, NSobol = 2^<N>\n";
                    exit(EXIT_FAILURE);
                }
                NSobol = 1<<p2;
                break;
            case 'S':
                p2 = atoi(optarg);
                if (p2>32) {
                    cerr << "-S <N>: <N>is too large, sobol_skip = 2^<N>\n";
                    exit(EXIT_FAILURE);
                }
                sobol_skip = 1<<p2;
                break;
            case '2':
                Nmodes2 = atoi(optarg);
                break;
            case '3':
                Nmodes3 = atoi(optarg);
                break;
            case 'r':
                rng_in.open(optarg);
                break;
            case 's':
                spectrum_file = optarg;
                break;
            case 'c':
                continue_from_file = optarg;
                continue_from = strtol(strrchr(optarg, '_')+1, NULL, 10);
                break;
            case 'p':
                h2o_potential = optarg;
                break;
            default:
                cerr << "Unknown option: " << ch << endl;
                exit(EXIT_FAILURE);
                break;
        }
    }
    
    argc -= optind;
    argv += optind;
    
    if (argc < 1) {
        cerr << "No input file supplied\n";
        exit(EXIT_FAILURE);
    }
    
    input_file = argv[0];
}

void SCP3(h2o::Potential& pot, vec& x0, vec& omega, mat& MUa)
{
    int Nmodes0 = x0.n_rows - 6;
    ScaledMatrixElements  sme(omega, Nmodes2, Nmodes3);
    
    mat M(sme.getBasisSize(), sme.getBasisSize());
    
    if (continue_from_file.empty()) {
        M.fill(0.0);
    }
    else {
        load_hdf5(continue_from_file, M);
        M *= -1.0;
        
        sme.addHODiagonal(M);
        M *= -(double)continue_from;
    }
    
    int nw = x0.n_rows / 9;
    
    vec y(Nmodes0), Vy(Nmodes0), r(9 * nw), Vr(9 * nw);
    
    double V;
    r = x0 * bohr;
    V = pot(nw, r.memptr(), Vr.memptr());
    V /= autokcalpmol;
    Vr *= bohr/autokcalpmol;
    
    cout << "V0 = "<<V<<endl;
    
    
    ofstream E0out_sd, E0out_t;
    
    if (continue_from_file.empty()) {
        E0out_sd.open("E0_singles+doubles.dat");
        E0out_t.open("E0_triples.dat");
    }
    else {
        char fname[128];
        
        sprintf(fname, "E0_singles+doubles_%07d.dat", continue_from);
        E0out_sd.open(fname);
        
        sprintf(fname, "E0_triples_%07d.dat", continue_from);
        E0out_t.open(fname);
    }
    
    fixed(E0out_sd);
    E0out_sd.precision(10);
    
    fixed(E0out_t);
    E0out_t.precision(10);
    
    
    for (int i=0; i<NSobol; i++) {
        if (rng_in.is_open()) {
            for (int j=0; j<Nmodes0; j++) {
                rng_in >> y[j];
            }
        }
        else {
            sobol::std_normal(y.n_rows, &sobol_skip, y.memptr());
        }
        
        if (i+1 <= continue_from) continue;
        
        y /=  sqrt(2.0);
        r = bohr *(MUa*y + x0);
        V = pot(nw, r.memptr(), Vr.memptr());
        V /= autokcalpmol;
        Vr *= bohr/autokcalpmol;
        
        if (isnan(V)) {
            cerr << "V=NaN at i=" << i << endl;
            i--;
            continue;
        }
        
        V -= 0.5 * dot(y, omega%y);
        Vy = MUa.t() * Vr - omega % y;
        sme.addEpot(y, V, Vy, M);
        
        if ( (i+1)%(1<<14)==0) {
            cout << i+1 << endl;
            mat Mout = M / (i+1);
            sme.addHODiagonal(Mout);
            
            vec eigvals = eig_sym(Mout.submat(0,0, Nmodes0, Nmodes0));
            double E0[3];
            
            E0[0] = eigvals[0]*autocm;
            
            int Nstates2 = sme.getSubBasisSize(2);
            eigvals = eig_sym(Mout.submat(0, 0, Nstates2-1, Nstates2-1));
            E0[1] = eigvals[0]*autocm;
            
            E0out_sd << i+1 << " " << E0[0] <<" "<< E0[1] <<" ";
            E0out_sd << E0[1] - E0[0] << endl;
            
            if ( (i+1)%(1<<17)==0) {
                char s[128];
                sprintf(s, "sM_%07d.h5", i+1);
                save_hdf5(Mout, s);
                
                eigvals = eig_sym(Mout);
                E0[2] = eigvals[0]*autocm;
                
                E0out_t << i+1 <<" "<< E0[2] <<" "<< E0[2] - E0[0] <<endl;
            }
                //dump_spectrum(TIP4P_charges, MUa, Mout, sme, specout, dipoleout, E0out);
        }
    }
    E0out_sd.close();
    E0out_t.close();
    
    M /= NSobol;
}


void OHHOHH(vec& mass, vec& r, mat& H)
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
        else if (mass[i] == Hmass) {
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
    
    mat H_tmp = H.cols(p);
    H = H_tmp.rows(p);
}


int main (int argc, char *  argv[]) {
    process_options(argc, argv);

    int N;
    mat H;
    vec mass, x0;
        // if (harmonic_approximation) {
            //load_from_vladimir(input_file, N, mass, x0);
            //    H = getHessian(x0);
            //   massScaleHessian(mass, H);
            //   }
            //  else {

    load_from_vladimir(input_file, N, mass, x0, H);
        //  }
    
    OHHOHH(mass, x0, H);
    
    vec omegasq0;
    mat U;
    eig_sym(omegasq0, U, H);

    ofstream sout("omega0.dat");
    sout << sqrt(abs(omegasq0))*autocm << endl;
    sout.close();

    if (Nmodes2 == 0) Nmodes2 = 3*N - 6;

    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);


    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }
    mat MUaT=MUa.t();

    h2o::Potential *pot;
    if (h2o_potential == "whbb") {
#ifdef HAVE_BOWMAN
        pot = new h2o::bowman();
#else
        cerr << "Support for Bowman's WHBB was not included." << endl;
        exit(EXIT_FAILURE);
#endif
    }
    else if (h2o_potential == "qtip4pf") {
        pot = new h2o::qtip4pf();
    }
    else {
        cerr << "Unknown H2O potential: " << h2o_potential << endl;
        exit(EXIT_FAILURE);
    }

        
    if ( !spectrum_file.empty() ) {
        mat M;
        load_hdf5(spectrum_file, M);

        exit(EXIT_SUCCESS);
    }

    SCP3(*pot, x0, omega, MUa);
    
    rng_in.close();
}
