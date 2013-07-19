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

#include "SCP1.h"


using namespace std;
using namespace arma;

static struct option program_options[] = {
    { "NSobol", required_argument, NULL, 'N'} ,
    { "skip", required_argument, NULL, 'S'} ,
    { "rng-file", required_argument, NULL, 'r'} ,
    { "PES", required_argument, NULL, 'p'},
    { "max-itertions", required_argument, NULL, 'L'},
    {NULL, 0, NULL, 0}
};

static int64_t NSobol=1<<20;
long long int sobol_skip=1<<22;

static string input_file;

static ifstream rng_in;

static int max_iterations=200;

string  h2o_pes("qtip4pf");


void process_options(int argc,  char *  argv[])
{
    int ch;
    while ( (ch = getopt_long(argc, argv, "S:N:r:p:L:", program_options, NULL)) != -1) {
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
            case 'p':
                h2o_pes = optarg;
                break;
            case 'L':
                max_iterations = atoi(optarg);
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

int main (int argc, char *  argv[])
{
    MPI::Init(argc, argv);

    process_options(argc, argv);

#ifdef HAVE_BOWMAN
    if (h2o_pes == "whbb" || h2o_pes == "hbb2-pol") {
        ps::pot_nasa_init();
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
        h2o::x3b_bits::load("x3b.nc");
    }
#endif

    mat H;
    vec mass, x0;

    load_from_vladimir(input_file, mass, x0, H);
    cout << H.n_elem << endl;
    OHHOHH(mass, x0, H);
        
    SCP1 scp1(mass, h2o_pes, NSobol, !H.is_empty());
    double F0 = scp1(x0, 0.0, H, max_iterations);
    
    cout << "F0 = " << F0 << endl;

    MPI::Finalize();
}
