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

#include <tr1/random>

#include <getopt.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"
#include "F90.h"
#include "Potentials.h"
#include "ScaledMatrixElements.h"

using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "doubles", required_argument, NULL, '2'},
    { "triples", required_argument, NULL, '3'},
    { "spectrum", required_argument, NULL, 's'},
    { "whbb", no_argument, NULL, 'W'},
    {NULL, 0, NULL, 0}
};

static unsigned int  Nmodes2=0, Nmodes3=0;
static int64_t NSobol=1<<20, sobol_skip=1<<30;
static string input_file;
static string spectrum_file;

static string continue_from_file;
static unsigned int continue_from;
static ifstream rng_in;

bool    use_whbb = false;
Potential* pot;


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

vec TIP4P_charges(size_t N)
{
    static const double qM = 1.1128;
    vec q(N);

    for (int i=0; i<N; ) {
        q[i++] = -qM;
        q[i++] = 0.5*qM;
        q[i++] = 0.5*qM;
    }
    
    return q;
}
/*
mat getHessian(vec& r)
{
    int N = r.n_cols;

    vec Vrp(N);
    mat H(N,N);

    for (int i=0; i<N; i++) {
        double ri0 = r[i];
        
        r[i] = ri0 + s;
        TIP4P_UF(N, rp.memptr(), &V, H.col(i).memptr());

        r[i] = ri0 - s;
        TIP4P_UF(N, rp.memptr(), &V, Vrp.memptr());

        H.col(i) = (H.col(i) - Vrp) / 2.0;
        
        r[i] = ri0;
    }

    for (int i=1; i<N; i++) {
        H.col(i).subvec(0, i-1) = 0.5 * (H.col(i).subvec(0, i-1) 
                                         + H.row(i).subvec(0, i-1) );
    }
    
    return H;
}


void massScaleHessian(vec& mass, mat& H)
{
    vec isqrt_mass = 1.0/sqrt(mass);
    
    for (int i=0; i<H.n_cols; i++) {
        H.col(i) *= isqrt_mass * isqrt_mass[i];
    }
}
*/


void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "S:N:2:3:c:s:r:W", program_options, NULL)) != -1) {
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
            case 'W':
                use_whbb = true;
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




int main (int argc, char *  argv[]) {
    int N;
    mat H, U;
    vec mass, x0, omegasq0;

    process_options(argc, argv);

    if (use_whbb) {
        pot = new WHBB(1e-3);
    }
    else {
        pot = new qTIP4Pf();
    }

    
        // if (harmonic_approximation) {
            //load_from_vladimir(input_file, N, mass, x0);
            //    H = getHessian(x0);
            //   massScaleHessian(mass, H);
            //   }
            //  else {
        load_from_vladimir(input_file, N, mass, x0, H);
        //  }
    
    eig_sym(omegasq0, U, H);

    ofstream sout("omega0.dat");
    sout << sqrt(abs(omegasq0))*autocm << endl;
    sout.close();

    int Nmodes0 = 3*N - 6;
    if (Nmodes2 == 0) Nmodes2 = 3*N - 6;

    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);


    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }
    mat MUaT=MUa.t();

    int Nstates2 = 1 + Nmodes0 + Nmodes2*(Nmodes2 + 1)/2;
    int Nstates = Nstates2 + (Nmodes3 * (1 + Nmodes3) * (2 + Nmodes3))/6;

    mat M(Nstates, Nstates);
    ScaledMatrixElements  sme(omega, Nmodes2, Nmodes3);

    
    if (continue_from_file.empty()) {
        M.fill(0.0);
    }
    else {
        load_hdf5(continue_from_file, M);
        M *= -1.0;
        sme.addHODiagonal(M);
        M *= -(double)continue_from;
    }
    
    vec charge = TIP4P_charges(N);

    if ( !spectrum_file.empty() ) {
        load_hdf5(spectrum_file, M);

        string tname = "freq" + spectrum_file.substr(1);
        ofstream specout(tname.c_str());
        ofstream dipoleout(("dipole" + spectrum_file.substr(1)).c_str());

        dump_spectrum(charge, MUa, M, sme, specout, dipoleout, cout, "C" + spectrum_file.substr(1));
        exit(EXIT_SUCCESS);
    }

    vec y(Nmodes0), Vy(Nmodes0), r(3*N), Vr(3*N);
    double V;


    int NO = N/3;
    (*pot)(x0, V, Vr);
    cout << "V0 = "<<V<<endl;

    /*
    ofstream specout("sfreq.dat");
    ofstream dipoleout("sdipole.dat");
    */

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
            sobol_stdnormal_c(y.n_rows, &sobol_skip, y.memptr());
        }
        
        if (i+1 <= continue_from) continue;

        y /=  sqrt(2.0);
        r = MUa*y + x0;
        (*pot)(r, V, Vr);

        if (isnan(V)) {
            cerr << "V=NaN at i=" << i << endl;
            i--;
            continue;
        }

        V -= 0.5 * dot(y, omega%y);
        Vy = -MUaT * Vr - omega % y;
        sme.addEpot(y, V, Vy, M);

        if ( (i+1)%(1<<14)==0) {
            cout << i+1 << endl;
            mat Mout = M / (i+1);
            sme.addHODiagonal(Mout);

            vec eigvals = eig_sym(Mout.submat(0,0, Nmodes0, Nmodes0));
            double E0[3];

            E0[0] = eigvals[0]*autocm;

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
    /*
    specout.close();
    dipoleout.close();
    */
    rng_in.close();

    M /= NSobol;
}
