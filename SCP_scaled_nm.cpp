/*
 *  SCP_scaled_nm.cpp
 *  SCP_Double_Excitations
 *
 *  Created by Ionuț Georgescu on 2/5/13.
 *  Copyright 2013 UCI Chemistry. All rights reserved.
 *
 */

#include <armadillo>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <string>

#include <tr1/random>

#include <getopt.h>

#include "ScaledMatrixElements.h"


extern "C" {
    void sobol_stdnormal_c(int64_t d, int64_t *skip, void *x);
    void TIP4P_UF(int N, double *r, double *U, double *UX);
}


using namespace std;
using namespace arma;


static const double bohr=0.52917721092;
static const double autocm=2.194746313e5;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "doubles", required_argument, NULL, '2'},
    { "triples", required_argument, NULL, '3'},
    { "spectrum", required_argument, NULL, 's'},
    {NULL, 0, NULL, 0}
};

static unsigned int  Nmodes2=0, Nmodes3=0;
static int64_t NSobol=1<<20, sobol_skip=1<<30;
static string input_file;
static string spectrum_file;

void process_options(int argc,  char *  argv[])
{
    int ch, i, j;
    while ( (ch = getopt_long(argc, argv, "N:2:3:c:s:", program_options, NULL)) != -1) {
        switch (ch) {
            case 'N':
                NSobol = atoi(optarg);
                break;
            case 'S':
                sobol_skip = 1<<atoi(optarg);
                break;
            case '2':
                Nmodes2 = atoi(optarg);
                break;
            case '3':
                Nmodes3 = atoi(optarg);
                break;
            case 's':
                spectrum_file = optarg;
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




void dump_spectrum(vec &charges, mat &MU, mat &H, ScaledMatrixElements& me,
                   ostream& specout, ostream& dipoleout, ostream& E0out,
                   string Cout="")
{
    char s[256];
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


int main (int argc, char *  argv[]) {
    int N;
    mat H, U;
    vec mass, x0, omegasq0;

    process_options(argc, argv);
    load_from_vladimir(input_file, N, mass, x0, H);
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
    M.fill(0.0);
    ScaledMatrixElements  sme(omega, Nmodes2, Nmodes3);

    vec charge = TIP4P_charges(N);

    if ( !spectrum_file.empty() ) {
        M.load(spectrum_file, raw_ascii);

        string tname = "freq" + spectrum_file.substr(1);
        ofstream specout(tname.c_str());
        ofstream dipoleout(("dipole" + spectrum_file.substr(1)).c_str());

        dump_spectrum(charge, MUa, M, sme, specout, dipoleout, cout, "C" + spectrum_file.substr(1));
        exit(EXIT_SUCCESS);
    }

    vec y(Nmodes0), Vy(Nmodes0), r(3*N), Vr(3*N);
    double V;


    int NO = N/3;
    TIP4P_UF(NO, x0.memptr(), &V, Vr.memptr());
    cout << "V0 = "<<V<<endl;

    ofstream specout("sfreq.dat");
    ofstream dipoleout("sdipole.dat");

    ofstream E0out("E0.dat");
    fixed(E0out);
    E0out.precision(10);
    
    mat sobol_sequence;
    sobol_sequence.load("MatousekAffineOwen30.dat", raw_ascii);
    
    /*
    for (int i=0; i<NSobol; i++) {
        sobol_stdnormal_c(sobol_sequence.n_rows, &sobol_skip,
                          sobol_sequence.colptr(i));
    }
    
    tr1::mt19937 rng;
    for (uint64_t i=0; i<NSobol-1; i++) {
        uint64_t j = rng();
        j = (j * (NSobol - i)) >> 32;
        cout << i << " " << i+j << endl;
        sobol_sequence.swap_cols(i, i + j);
    }
     */
    
    for (int i=0; i<NSobol; i++) {
        y = sobol_sequence.row(i) / sqrt(2.0);
        r = MUa*y + x0;
        TIP4P_UF(NO, r.memptr(), &V, Vr.memptr());

        if (isnan(V)) {
            cerr << "V=NaN at i=" << i << endl;
            i--;
            continue;
        }

        V -= 0.5 * dot(y, omega%y);
        Vy = -MUaT * Vr - omega % y;
            //cout << y(Nmodes-2) << " " << y(Nmodes-1) << " " << Vy(Nmodes-2) << " " << Vy(Nmodes-1) << endl;
        sme.addEpot(y, V, Vy, M);

        if ( (i+1)%(1<<15)==0) {
            cout << i+1 << endl;
            mat Mout = M / (i+1);
            sme.addHODiagonal(Mout);

      //      if ( (i+1)%1<<17==0) {
                char s[100];
                sprintf(s, "sM_%07d.dat", i+1);
                Mout.save(s, raw_ascii);
      //      }

            vec eigvals = eig_sym(Mout.submat(0,0, Nmodes0, Nmodes0));
            double E0[3];

            E0[0] = eigvals[0]*autocm;

            eigvals = eig_sym(Mout.submat(0, 0, Nstates2-1, Nstates2-1));
            E0[1] = eigvals[0]*autocm;

            eigvals = eig_sym(Mout);
            E0[2] = eigvals[0]*autocm;

            E0out << i+1 << " " << E0[0] <<" "<< E0[1] <<" "<< E0[2] <<" ";
            E0out << E0[1] - E0[0] <<" "<< E0[2] - E0[0] << endl;
                //dump_spectrum(TIP4P_charges, MUa, Mout, sme, specout, dipoleout, E0out);
        }
    }
    specout.close();
    dipoleout.close();
    M /= NSobol;
}
