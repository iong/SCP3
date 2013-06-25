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

#include <getopt.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"
#include "ScaledMatrixElements.h"

#include "sobol.hpp"

#include "qtip4pf.h"
#include "ttm3f.h"

#ifdef HAVE_BOWMAN
#include "bowman-fortran.h"
#include "bowman.h"
#include "ps.h"
#endif



using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "doubles", required_argument, NULL, '2'},
    { "triples", required_argument, NULL, '3'},
    { "potential", required_argument, NULL, 'p'},
    { "select-modes", required_argument, NULL, 'm'},
    { "distributed", required_argument, NULL, 'd'},
    { "no-gradient", no_argument, NULL, 'G'},
    { "block-width", required_argument, NULL, 'B'},
    {NULL, 0, NULL, 0}
};

static unsigned int  Nmodes2=0, Nmodes3=0;
static int64_t NSobol=1<<20;
long long int sobol_skip=1<<30;

static string input_file;

static ifstream rng_in;

static uvec selected_modes;

static bool use_gradient = true;

static int proc_id = 0;
static int nb_proc = 1;

static int block_width = 1;

string  h2o_potential("qtip4pf");


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
    while ( (ch = getopt_long(argc, argv, "S:N:2:3:r:p:m:d:GB:", program_options, NULL)) != -1) {
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
            case 'p':
                h2o_potential = optarg;
                break;
            case 'm':
                parse_mode_description(optarg, selected_modes);
                break;
            case 'd':
                proc_id = atoi(optarg);
                nb_proc = atoi(strrchr(optarg, '/') + 1);
                break;
            case 'G':
                use_gradient = false;
                break;
            case 'B':
                block_width = atoi(optarg);
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


h2o::Potential *getPotential(const string& name)
{
    h2o::Potential *pot;
    if (name == "whbb") {
#ifdef HAVE_BOWMAN
        pot = new h2o::bowman();
#else
        cerr << "Support for Bowman's WHBB was not included." << endl;
        exit(EXIT_FAILURE);
#endif
    }
    else if (name == "qtip4pf") {
        pot = new h2o::qtip4pf();
    }
    else {
        cerr << "Unknown H2O potential: " << h2o_potential << endl;
        exit(EXIT_FAILURE);
    }

    return pot;
}



void SCP3_a(const string& h2o_potential, const vec& x0, const vec& omega, const mat& MUa, const uvec& modes)
{
    if (Nmodes2 > modes.n_rows || Nmodes3 > modes.n_rows) {
        cerr << "Number of double or triple excitations exceed number of modes."
        << endl;
        exit(EXIT_FAILURE);
    }
    
    int Nmodes0 = omega.n_rows;

    ScaledMatrixElements  sme(omega(modes), Nmodes2, Nmodes3);
    int Nstates = sme.getBasisSize();
    
    mat M(Nstates, Nstates);
    M.fill(0.0);
    
    ofstream E0out_sd, E0out_t;
    if (nb_proc > 1) {
        char s[32];
        sprintf(s, "E0_sd-%03d.dat", proc_id);
        E0out_sd.open(s);
    }
    else {
        E0out_sd.open("E0_sd.dat");
        E0out_t .open("E0_t.dat" );
    }


    fixed(E0out_sd);
    E0out_sd.precision(10);
    
    fixed(E0out_t);
    E0out_t.precision(10);
    
    double E0[3];
 
    int seq_len = NSobol / nb_proc;
    int seq_start = seq_len * proc_id;
    int seq_stop = seq_start + seq_len;

    int nw = x0.n_rows / 9;
    
    vec bV(block_width);
    mat by(Nmodes0, block_width), bVy(Nmodes0, block_width); 

    if (seq_start % block_width != 0) {
        cerr << "seq_start must me a multiple of block_width!\n";
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < seq_stop; i += block_width) {
        if (rng_in.is_open()) {
            for (int j=0; j < block_width; j++) {
                for (int k=0; k<Nmodes0; k++) {
                    rng_in >> by(k, j);
                }
            }
        }
        else {
            for (int j=0; j < block_width; j++) {
                sobol::std_normal(by.n_rows, &sobol_skip, by.colptr(j));
            }
        }

        if (i < seq_start ) continue;

        by /=  sqrt(2.0);

#pragma omp parallel
        {
            h2o::Potential *pot = getPotential(h2o_potential);

#pragma omp for schedule(static)
            for (int j = 0; j < block_width; j++) {
                double V;

                vec r = bohr *(x0 + MUa*by.col(j));
                
                if (use_gradient) {
                    vec Vr(r.n_rows);

                    V = (*pot)(nw, r.memptr(), Vr.memptr()) / autokcalpmol;
                    Vr *= bohr/autokcalpmol;
                    
                    bVy.col(j) = MUa.t() * Vr - omega % by.col(j);
                }
                else {
                    V = (*pot)(nw, r.memptr()) / autokcalpmol;
                }

                bV[j] = V - 0.5 * dot(by.col(j), omega % by.col(j));

            }

            delete pot;
        }
            
        if (use_gradient) {
            sme.addEpot(by.rows(modes), bV, bVy.rows(modes), M);
        }
        else {
            sme.addEpot(by.rows(modes), bV, M);
        }

            if ( (i+1)%(1<<14)==0) {
                mat Mout = M / (i+1 - seq_start);
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
                
                if ( (i+1)%(1<<17)==0 && nb_proc == 1) {
                    char s[128];
                    sprintf(s, "sM_%07d.h5", i+1);
                    save_hdf5(Mout, s);
                   /* 
                    eigvals = eig_sym(Mout);
                    E0[2] = eigvals[0]*autocm;
                    
                    E0out_t << i+1 <<" "<< E0[2] <<" "<< E0[2] - E0[0] <<endl;
                    */
                }
            }
    }
    E0out_sd.close();
    E0out_t.close();
    
    M /= seq_len;
    sme.addHODiagonal(M);
    M.diag() += 0.5 * (sum(omega) - sum(omega(modes)));
    
    if (nb_proc > 1) {
        char s[128];
        sprintf(s, "sM_%07d.h5", seq_stop);
        save_hdf5(M, s);
    }
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

    vec mass2 = mass(p);
    mass = mass2;
}


int main (int argc, char *  argv[]) {
    process_options(argc, argv);

    int N;
    mat H;
    vec mass, x0;

    load_from_vladimir(input_file, N, mass, x0, H);
    OHHOHH(mass, x0, H);
    
    vec omegasq0;
    mat U;
    eig_sym(omegasq0, U, H);

    if (proc_id == 0) {
        ofstream sout("omega0.dat");
        sout << sqrt(abs(omegasq0))*autocm << endl;
        sout.close();
    }

    if (Nmodes2 == 0) Nmodes2 = 3*N - 6;

    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);


    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }

#ifdef HAVE_BOWMAN
    if (h2o_potential == "whbb") {
        ps::pot_nasa_init();
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
    }
#endif

    if (selected_modes.is_empty() ) {
        selected_modes = linspace<uvec>(0, MUa.n_cols - 1, MUa.n_cols);
    }

    SCP3_a(h2o_potential, x0, omega, MUa, selected_modes);
    
    rng_in.close();
}
