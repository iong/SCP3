/*
 *  Upot.cpp
 *  SCP3
 *
 *  Created by Ionuț Georgescu on 4/9/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */


/*
 *  SCP_scaled_nm.cpp
 *  SCP_Double_Excitations
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

#include <getopt.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"
#include "F90.h"

using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "hydrogen-mass", required_argument, NULL, 'H'} ,
    {NULL, 0, NULL, 0}
};

static int64_t NSobol=1<<20, sobol_skip=1<<20;
static string input_file;

static ifstream rng_in;

void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "N:S:r:H:", program_options, NULL)) != -1) {
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



mat getHessian(vec& r, double s)
{
    int N = r.n_rows;
    int Nw = N/9;

    vec Vrp(N), Vrm(N);
    mat H(N,N);
    
    for (int i=0; i<N; i++) {
        double ri0 = r[i];
        double V;
        
        r[i] = ri0 + s;
        TIP4P_UF(Nw, r.memptr(), &V, Vrp.memptr());
        
        r[i] = ri0 - s;
        TIP4P_UF(Nw, r.memptr(), &V, Vrm.memptr());
        
        H.col(i) = (Vrm - Vrp) / (2.0*s);
        
        r[i] = ri0;
    }
    
    for (int i=1; i<N; i++) {
        H.col(i).subvec(0, i-1) = 0.5 * (H.col(i).subvec(0, i-1)
                                         + H.row(i).subvec(0, i-1).t() );
    }
    
    return H;
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
    x0 /= bohr;

    int N = x0.n_rows / 3;

    mat H = getHessian(x0, 0.01);
    massScaleHessian(mass, H);

    //  }
    
    mat U;
    eig_sym(omegasq0, U, H);
    
    int Nmodes0 = 3*N - 6;
    
    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);
    
    
    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }
    mat MUaT=MUa.t();
    
    double V, V0, V_avg;
    vec y(Nmodes0), Vy(Nmodes0), Vy_avg(Nmodes0), r(3*N), Vr(3*N), Vr0(3*N);
    
    V_avg = 0.0;
    Vy_avg.fill(0.0);
    
    
    int NO = N/3;
    TIP4P_UF(NO, x0.memptr(), &V0, Vr0.memptr());
    
    sobol_skip = 2*NSobol;
    for (int i=0; i<NSobol; i++) {
        if (rng_in.is_open()) {
            for (int j=0; j<Nmodes0; j++) {
                rng_in >> y[j];
            }
        }
        else {
            sobol_stdnormal_c(y.n_rows, &sobol_skip, y.memptr());
        }
        
        y /=  sqrt(2.0);
        r = MUa*y + x0;
        TIP4P_UF(NO, r.memptr(), &V, Vr.memptr());
        
        if (isnan(V)) {
            cerr << "V=NaN at i=" << i << endl;
            i--;
            continue;
        }
        
        V -= 0.5 * dot(y, omega%y);
        Vy = -MUaT * Vr - omega % y;
        
        V_avg += V;
        Vy_avg += Vy;
    }
    rng_in.close();
    
    V_avg /= NSobol;
    Vy_avg /= NSobol;
    
    /*
    cout << "V0 = "<< V0*autocm <<" cm^-1\tsum(omega)/4 ="<< 0.25*sum(omega)*autocm <<" cm^-1\n"
    <<"V0 + sum(omega)/4 =" << (V0 + 0.25*sum(omega))*autocm <<" cm^-1\n"
    <<"<V> = " << V_avg*autocm<<" cm^-1\n"
    <<"<V> - sum(omega)/4 - V0 =" << (V_avg - 0.25*sum(omega) - V0)*autocm <<" cm^-1\n"
    <<"<V> + sum(omega)/4 = " << (V_avg + 0.25*sum(omega))*autocm <<" cm^-1\n";
    */
    
    V0 *= autocm;
    omega *= autocm;
    V_avg *= autocm;
    
    cout  << Hmass/1837.15137 <<" "<< V0 <<" "<< sum(omega)/4.0 <<" "<< V_avg
        <<" "<< V_avg + sum(omega)/4.0 <<endl;
}
