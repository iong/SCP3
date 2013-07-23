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

#include <mpi.h>

#include <armadillo>

#include "Constants.h"
#include "DiskIO.h"
#include "Hessian.h"
#include "ScaledMatrixElements.h"

#include "sobol.hpp"

#include "h2o.h"

using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "doubles", required_argument, NULL, '2'},
    { "triples", required_argument, NULL, '3'},
    { "PES", required_argument, NULL, 'p'},
    { "select-modes", required_argument, NULL, 'm'},
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

static int myrank = 0;
static int nprocs = 1;

static int block_width = 1024;

string  h2o_pes("qtip4pf");


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
    while ( (ch = getopt_long(argc, argv, "S:N:2:3:r:p:m:GB:", program_options, NULL)) != -1) {
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
                h2o_pes = optarg;
                break;
            case 'm':
                parse_mode_description(optarg, selected_modes);
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


void SCP3_a(const string& h2o_pes, const vec& x0, const vec& omega, const mat& MUa, const uvec& modes)
{
    if (Nmodes2 > modes.n_rows || Nmodes3 > modes.n_rows) {
        cerr << "Number of double or triple excitations exceed number of modes."
        << endl;
        exit(EXIT_FAILURE);
    }
    
    int Nmodes0 = omega.n_rows;

    ScaledMatrixElements  sme(omega(modes), Nmodes2, Nmodes3);
    int Nstates = sme.getBasisSize();
    
    mat M;
    ofstream E0out_sd, E0out_t;
    if (myrank ==0) {
        M.zeros(Nstates, Nstates);
        E0out_sd.open("E0_sd.dat");
        E0out_t .open("E0_t.dat" );

        fixed(E0out_sd);
        E0out_sd.precision(10);
        
        fixed(E0out_t);
        E0out_t.precision(10);
    }
    
 
    int nw = x0.n_rows / 9;
    
    int local_work = block_width / nprocs;
    vec bV(local_work), global_bV;
    mat by(Nmodes0, local_work), global_by,
        bVy(Nmodes0, local_work), global_bVy; 

    for (int i = 0; i < NSobol; i += block_width) {
        MPI::COMM_WORLD.Barrier();

        int local_start = i + local_work * myrank;
        int local_stop = local_start + local_work;

        if (rng_in.is_open()) {
            for (int j=i; j < local_start; j++) {
                double tmp;
                for (int k=0; k<Nmodes0; k++) rng_in >> tmp; 
            }
            for (int j=0; j < local_work; j++) {
                for (int k=0; k<Nmodes0; k++) {
                    rng_in >> by(k, j);
                }
            }
        }
        else {
            for (int j=0; j < local_work; j++) {
                long long int current_skip = sobol_skip + local_start + j;
                sobol::std_normal(by.n_rows, &current_skip, by.colptr(j));
            }
        }

        by /=  sqrt(2.0);

        h2o::PES *pot = h2o::PESFromString(h2o_pes);

        for (int j = 0; j < local_work; j++) {
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

        if (myrank == 0) {
            global_bV.zeros(block_width);
            global_by.zeros(Nmodes0, block_width);
            global_bVy.zeros(Nmodes0, block_width);
        }

        MPI::COMM_WORLD.Gather(bV.memptr(), local_work, MPI::DOUBLE,
                global_bV.memptr(), local_work, MPI::DOUBLE, 0);

        MPI::COMM_WORLD.Gather(
                by.memptr(), by.n_rows * local_work, MPI::DOUBLE,
                global_by.memptr(), by.n_rows * local_work, MPI::DOUBLE, 0);

        MPI::COMM_WORLD.Gather(
                bVy.memptr(), bVy.n_rows * local_work, MPI::DOUBLE,
                global_bVy.memptr(), bVy.n_rows * local_work, MPI::DOUBLE, 0);


        if (myrank != 0) continue;

        cout << "finished pot eval\n";
            
        if (use_gradient) {
            sme.addEpot(global_by.rows(modes), global_bV,
                    global_bVy.rows(modes), M);
        }
        else {
            sme.addEpot(global_by.rows(modes), global_bV, M);
        }

        if ( (i+block_width)%(1<<14)==0) {
            double E0[3];

            mat Mout = M / (i + block_width);
            sme.addHODiagonal(Mout);
            Mout.diag() += 0.5 * (sum(omega) - sum(omega(modes)));

            int Nstates1 = sme.getSubBasisSize(1);
            vec eigvals = eig_sym(Mout.submat(0,0, Nstates1-1, Nstates1-1));
            E0[0] = eigvals[0]*autocm;
            
            int Nstates2 = sme.getSubBasisSize(2);
            eigvals = eig_sym(Mout.submat(0,0, Nstates2-1, Nstates2-1));
            E0[1] = eigvals[0]*autocm;
            
            E0out_sd << i+block_width << " " << E0[0] <<" "<< E0[1] <<" "
            << E0[1] - E0[0] << endl;
            
            if ( (i+block_width)%(1<<17)==0) {
                char s[128];
                sprintf(s, "sM_%07d.h5", i+block_width);
                save_hdf5(Mout, s);
               /* 
                eigvals = eig_sym(Mout);
                E0[2] = eigvals[0]*autocm;
                
                E0out_t << i+1 <<" "<< E0[2] <<" "<< E0[2] - E0[0] <<endl;
                */
            }
        }
    }

    if (myrank != 0) return;

    E0out_sd.close();
    E0out_t.close();
    
    M /= NSobol;
    sme.addHODiagonal(M);
    M.diag() += 0.5 * (sum(omega) - sum(omega(modes)));
}

int main (int argc, char *  argv[])
{
    int provided = MPI::Init_thread(argc, argv, MPI::THREAD_FUNNELED);

    if (provided < MPI::THREAD_FUNNELED) {
        cerr << "Fucked!\n" << provided <<" "<< MPI::THREAD_FUNNELED <<endl;
        exit(EXIT_FAILURE);
    }

    myrank   = MPI::COMM_WORLD.Get_rank();
    nprocs = MPI::COMM_WORLD.Get_size();

    process_options(argc, argv);

#ifdef HAVE_BOWMAN
    if (h2o_pes == "whbb") {
        ps::pot_nasa_init();
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
    }
    else if (h2o_pes == "hbb2-pol") {
        ps::pot_nasa_init();
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
        h2o::x3b_bits::load("x3b.nc");
    }
#endif

    mat H;
    vec mass, x0;

    load_from_vladimir(input_file, mass, x0, H);
    OHHOHH(mass, x0, H);
    x0 /= bohr;

        
    vec omegasq0;
    mat U;
    eig_sym(omegasq0, U, H);

    if (myrank == 0) {
        ofstream sout("omega0.dat");
        sout << sqrt(abs(omegasq0))*autocm << endl;
        sout.close();
    }

    int N = x0.n_rows / 3;
    if (Nmodes2 == 0) Nmodes2 = 3*N - 6;

    vec omega = sqrt(omegasq0.rows(6, 3*N - 1));
    vec alpha = sqrt(omega);
    mat MUa = U.cols(6, 3*N-1);


    for (int i=0; i<MUa.n_cols; i++) {
        MUa.col(i) /= sqrt(mass)*alpha(i);
    }

    if (selected_modes.is_empty() ) {
        selected_modes = linspace<uvec>(0, MUa.n_cols - 1, MUa.n_cols);
    }
    
    SCP3_a(h2o_pes, x0, omega, MUa, selected_modes);
    
    rng_in.close();

    MPI::Finalize();
}
