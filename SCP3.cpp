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
#include <iterator>
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
    { "potential", required_argument, NULL, 'p'},
    { "select-modes", required_argument, NULL, 'm'},
    { "VMD", no_argument, NULL, 'V'},
    { "freeze-deselected", no_argument, NULL, 'f'},
    {NULL, 0, NULL, 0}
};

static unsigned int  Nmodes2=0, Nmodes3=0;
static int64_t NSobol=1<<20;
long long int sobol_skip=1<<30;

static string input_file;
static string spectrum_file;

static string continue_from_file;
static unsigned int continue_from;
static ifstream rng_in;

static uvec selected_modes;

static bool vmd_normal_modes = false;
static bool freeze_deselected = false;

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


void parse_mode_description(const string& desc, uvec& selected_modes)
{
    istringstream iss(desc);
    vector<int> t_modes;
    
    int start, stop;
    while (iss.good()) {
        if (iss.peek() == ',') {
            iss.get();
        }
        
        iss >> start;
        stop = start;
        
        if (iss.peek() == '-') {
            iss.get();
            iss >> stop;
        }
        
        for (;start <= stop; start++) {
            t_modes.push_back(start);
        }
    }
    
    selected_modes.set_size ( t_modes.size() );
    copy(t_modes.begin(), t_modes.end(), selected_modes.begin());
}


void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "S:N:2:3:c:s:r:p:m:Vf", program_options, NULL)) != -1) {
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
            case 'm':
                parse_mode_description(optarg, selected_modes);
                break;
            case 'V':
                vmd_normal_modes = true;
                break;
            case 'f':
                freeze_deselected = true;
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

void potentialMap(h2o::Potential& pot, vec& x0, mat& MUa, int nb_dim_pts, double grid_size)
{
    int nw = x0.n_rows / 9;
    
    int Nmodes = MUa.n_cols;

    fcube V(nb_dim_pts, nb_dim_pts, Nmodes*(Nmodes-1));

    for (int ni=0; ni < Nmodes; ni++) {
        for (int nj=0; nj < Nmodes; nj++) {
            if (ni == nj) {
                continue;
            }
            
            vec y(Nmodes);
            vec r(MUa.n_rows);
            
            y.fill(0.0);
            
            for (int i=0; i<nb_dim_pts; i++) {
                y[ni] = (-0.5  + (double)(i)/(double)(nb_dim_pts - 1)) * grid_size;
                for (int j=0; j<nb_dim_pts; j++) {
                    y[nj] = (-0.5  + (double)(j)/(double)(nb_dim_pts - 1)) * grid_size;
                    r = bohr * (MUa * y + x0);
                    V(j, i, ni*(Nmodes-1) + nj) = pot(nw, r.memptr()) / autokcalpmol;
                }
            }
        }
    }

    save_hdf5(V, "V.h5");
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

                /*
                eigvals = eig_sym(Mout);
                E0[2] = eigvals[0]*autocm;
                
                E0out_t << i+1 <<" "<< E0[2] <<" "<< E0[2] - E0[0] <<endl;
                */
            }
                //dump_spectrum(TIP4P_charges, MUa, Mout, sme, specout, dipoleout, E0out);
        }
    }
    E0out_sd.close();
    E0out_t.close();
    
    M /= NSobol;
}

void SCP3_a(h2o::Potential& pot, const vec& x0, const vec& omega, const mat& MUa, const uvec& modes)
{
    if (Nmodes2 > modes.n_rows || Nmodes3 > modes.n_rows) {
        cerr << "Number of double or triple excitations exceed number of modes."
        << endl;
        exit(EXIT_FAILURE);
    }
    
    int Nmodes0 = omega.n_rows;
    
    cout << MUa.n_cols <<", "<< modes.t() << endl;

    ScaledMatrixElements  sme(omega(modes), Nmodes2, Nmodes3);
    int Nstates = sme.getBasisSize();
    
    mat M(Nstates, Nstates);
    M.fill(0.0);

    int nw = x0.n_rows / 9;
    
    vec y(Nmodes0), Vy(Nmodes0), r(9 * nw), Vr(9 * nw);
    
    double V;
    r = x0 * bohr;
    V = pot(nw, r.memptr(), Vr.memptr());
    V /= autokcalpmol;
    Vr *= bohr/autokcalpmol;


    ofstream E0out_sd("E0_sd.dat"), E0out_t("E0_t.dat");

    fixed(E0out_sd);
    E0out_sd.precision(10);
    
    fixed(E0out_t);
    E0out_t.precision(10);
    
    double E0[3];
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
        sme.addEpot(y(modes), V, Vy(modes), M);

        
        if ( (i+1)%(1<<14)==0) {
            mat Mout = M / (i+1);
            sme.addHODiagonal(Mout);
            Mout.diag() += 0.5 * (sum(omega) - sum(omega(modes)));
            

            
            int Nstates1 = sme.getSubBasisSize(1);
            vec eigvals = eig_sym(Mout.submat(0,0, Nstates1-1, Nstates1-1));
            E0[0] = eigvals[0]*autocm;
            
            int Nstates2 = sme.getSubBasisSize(2);
            eigvals = eig_sym(Mout.submat(0,0, Nstates2-1, Nstates2-1));
            E0[1] = eigvals[0]*autocm;
            
            E0out_sd << i+1 << " " << E0[0] <<" "<< E0[1] <<" "
                << E0[1] - E0[0] << endl;
            
            if ( (i+1)%(1<<17)==0) {
                char s[128];
                sprintf(s, "sM_%07d.h5", i+1);
                save_hdf5(Mout, s);

                eigvals = eig_sym(Mout);
                E0[2] = eigvals[0]*autocm;
                
                E0out_t << i+1 <<" "<< E0[2] <<" "<< E0[2] - E0[0] <<endl;
            }
        }
    }
    E0out_sd.close();
    E0out_t.close();

    M /= NSobol;
    sme.addHODiagonal(M);
    M.diag() += 0.5 * (sum(omega) - sum(omega(modes)));
    
        //cout << E0[2] - M(0,0)*autocm << endl;
        //cout << endl << M*autocm << endl;
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


void display_normal_modes(string& fin, vec& x0, mat U)
{
    float nm_arrow_len = 1.0;
    float nm_arrow_radius = 0.1;
    float nm_arrow_tip_radius = 2.0 * nm_arrow_radius;
    float nm_arrow_tip_len = 1.5 * nm_arrow_tip_radius;
    
        //x0 *= bohr;
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

    if (vmd_normal_modes) {
        display_normal_modes(input_file, x0, U.cols(6, 3*N-6) );
    }
        //potentialMap(*pot, x0, MUa, 65, 8.0);
    if (selected_modes.is_empty() ) {
        selected_modes = linspace<uvec>(0, MUa.n_cols - 1, MUa.n_cols);
    }
     
    if (freeze_deselected) {
        SCP3_a(*pot, x0, omega(selected_modes), MUa.cols(selected_modes),
               linspace<uvec>(0, selected_modes.n_rows - 1, selected_modes.n_rows));
    }
    else {
        SCP3_a(*pot, x0, omega, MUa, selected_modes);
    }
    
    rng_in.close();
}
