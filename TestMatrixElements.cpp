/*
 *  TestMatrixElements.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 3/15/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <armadillo>

using namespace arma;

#include <getopt.h>

#include "TestScaledMatrixElements.h"


extern "C" {
    void sobol_stdnormal_c(int64_t d, int64_t *skip, void *x);
    void TIP4P_UF(int N, double *r, double *U, double *UX);
}


using namespace std;
using namespace arma;


static const double bohr=0.52917721092;
static const double autocm=2.194746313e5;
static const double qM = 1.1128;


static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    {NULL, 0, NULL, 0}
};

static unsigned int NSobol=1000000, Nmodes2=0, Nmodes3=0;
static string input_file, continue_from;
static string spectrum_file;

void process_options(int argc,  char *  argv[])
{
    int ch, i, j;
    while ( (ch = getopt_long(argc, argv, "N:", program_options, NULL)) != -1) {
        switch (ch) {
            case 'N':
                NSobol = atoi(optarg);
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

#include <cstdio>

void init_fortran_runtime(int argc, char *argv[])
{
#ifdef __GNU__
    _gfortran_set_args(argc, argv);
    _gfortran_set_args(1<<6-1);|
#endif
}

#include <xmmintrin.h>

int main(int argc, char *argv[])
{
    int N;
    mat H, U;
    vec mass, x0, omegasq0;
    
    init_fortran_runtime(argc, argv);
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    process_options(argc, argv);
    load_from_vladimir(input_file, N, mass, x0, H);
    eig_sym(omegasq0, U, H);
    
    
    int Nmodes0 = 3*N - 6;
    Nmodes2 = 3*N - 6;
    Nmodes3 = 3*N - 6;
    
    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);
    
    
    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }
    mat MUaT=MUa.t();
    
    TestScaledMatrixElements  sme(omega, Nmodes2, Nmodes3);

    int64_t sobol_skip=1<<int(ceil(log2((double)NSobol)));

    vec y(Nmodes0), Vy(Nmodes0), r(3*N), Vr(3*N);
    double V;
    
    int NO = N/3;
    TIP4P_UF(NO, x0.memptr(), &V, Vr.memptr());
    cout << "V0 = "<<V<<endl;

    vec weird_triples;
    ofstream test_out("weird_triples.dat");
    for (int i=0; i<NSobol; i++) {
        sobol_stdnormal_c(y.n_rows, &sobol_skip, y.memptr());
        y = y / (sqrt(2.0));
        r = MUa*y + x0;
        TIP4P_UF(NO, r.memptr(), &V, Vr.memptr());
        
        if (isnan(V)) {
            cerr << "V=NaN at i=" << i << endl;
            i--;
            continue;
        }
        
        V -= 0.5 * dot(y, omega%y);
        Vy = -MUaT * Vr - omega % y;
        
        sme.addEpot(y, V, Vy, weird_triples);
        
        if (i+1>130000) {
            printf("%d %g %g %g %g %g %g %g %g %g %g\n", i, y[0], y[10], y[15], y[17], V, Vy[0], Vy[10], Vy[15], Vy[17], sme.t0);
        }
        
        if ( (i+1)%10000==0) {
            cout << i+1 << endl;
            vec vavg = weird_triples / (i+1);
            
            
            test_out << i+1;
            for (int j=0; j<vavg.n_rows; j++) {
                test_out << " " << vavg[j];
            }
            test_out << endl;
                //dump_spectrum(TIP4P_charges, MUa, Mout, sme, specout, dipoleout, E0out);
        }
    }
    test_out.close();
}