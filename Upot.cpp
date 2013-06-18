/*
 *  Upot.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/9/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include <getopt.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"

#include "sobol.hpp"

#include "qtip4pf.h"
#include "ttm3f.h"

#ifdef HAVE_BOWMAN
#include "bowman.h"
#include "bowman-fortran.h"
#endif

using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "hydrogen-mass", required_argument, NULL, 'H'} ,
    { "potential", required_argument, NULL, 'p'},
    {NULL, 0, NULL, 0}
};

static int64_t NSobol=1<<20;
static long long int sobol_skip=1<<20;
static string input_file;

static ifstream rng_in;

string  h2o_potential("qtip4pf");

void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "N:S:r:H:p:", program_options, NULL)) != -1) {
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
            case 'r':
                rng_in.open(optarg);
                break;
            case 'H':
                Hmass *= atof(optarg);
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



void getHessian(h2o::Potential& pot, vec& r_au, double s, mat& H)
{
    int N = r_au.n_rows;
    int Nw = N/9;

    vec Vrp(N), Vrm(N);
    H.set_size(N,N);
    
    vec r(r_au * bohr);

    for (int i=0; i<N; i++) {
        double ri0 = r[i];
        double V;
        
        r[i] = ri0 + s * bohr;
        V = pot(Nw, r.memptr(), Vrp.memptr());
        
        r[i] = ri0 - s * bohr;
        V = pot(Nw, r.memptr(), Vrm.memptr());
        
        H.col(i) = (Vrp - Vrm) * bohr / (2.0*s*autokcalpmol);
        
        r[i] = ri0;
    }
    
    for (int i=1; i<N; i++) {
        H.col(i).subvec(0, i-1) = 0.5 * (H.col(i).subvec(0, i-1)
                                         + H.row(i).subvec(0, i-1).t() );
    }
}


void massScaleHessian(vec& mass, mat& H)
{
    vec isqrt_mass = sqrt(mass);
    for (int i=0; i<isqrt_mass.n_rows; i++) {
	    isqrt_mass[i] = 1.0/isqrt_mass[i];
    }
    
    for (int i=0; i<H.n_cols; i++) {
        H.col(i) %= isqrt_mass * isqrt_mass[i];
    }
}


void print_header(ostream& os)
{
    os <<"Hmass\tV0\tEkin=sum(omega)/4\t<V>\t<V>+Ekin\n";
}


void PrintHA(double V0, vec& omega)
{
    cout.precision(10);
    cout << "V0 = " << V0 * autocm << " cm^-1"
        << "\t\t= " << V0 * autokcalpmol << " kcal/mol\n"
        << "E0 = " << (V0 + 0.5*sum(omega)) * autocm << " cm^-1"
        << "\t\t= " << (V0 + 0.5*sum(omega)) * autokcalpmol << " kcal/mol\n"
        << "omega = " << omega.t()*autokcalpmol << endl;
}


void OHHOHH(vec& mass, vec& r)
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
    r = r_tmp(p);
}

int main (int argc, char *  argv[]) {
    vec mass, x0, omegasq0;
    
    process_options(argc, argv);
    
    // if (harmonic_approximation) {
    //load_from_vladimir(input_file, N, mass, x0);
    //    H = getHessian(x0);
    //   massScaleHessian(mass, H);
    //   }
    //  else {
    //load_from_vladimir(input_file, N, mass, x0, H);
    
    load_xyz(input_file, mass, x0);
    OHHOHH(mass, x0);
    x0 /= bohr;

    int N = x0.n_rows / 3;

    h2o::Potential *pot;
    if (h2o_potential == "whbb") {
#ifdef HAVE_BOWMAN
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
        pot = new h2o::bowman();

#else
        cerr << "Support for Bowman's WHBB was not included." << endl;
        exit(EXIT_FAILURE);
#endif
    }
    else if (h2o_potential == "qtip4pf") {
        pot = new h2o::qtip4pf();
    }
    else if (h2o_potential == "ttm3f") {
        pot = new h2o::ttm3f();
    }
    else {
        cerr << "Unknown H2O potential: " << h2o_potential << endl;
        exit(EXIT_FAILURE);
    }
   /* 
    mat H;
    getHessian(*pot, x0, 1e-3, H);
    massScaleHessian(mass, H);

    
    mat U;
    eig_sym(omegasq0, U, H);
    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));

    cout << "Equilibrium test: " << omegasq0(span(0,5)).t() << endl;
    */ 
    int nw = N/3;
    vec r = x0 * bohr;
    double V0 = (*pot)(nw, r.memptr())  / autokcalpmol;
    //PrintHA(V0, omega);
    cout << "V0 = " << V0 * autocm << " cm^-1\n";
    
    exit(EXIT_SUCCESS);
    
}
