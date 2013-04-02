#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <armadillo>
#include <string>
#include <boost/program_options.hpp>

#include "DiskIO.h"
#include "F90.h"
#include "MatrixElements.h"


using namespace std;
using namespace arma;
namespace po = boost::program_options;

static const double bohr=0.52917721092;
static const double autocm=2.194746313e5;
static const double qM = 1.1128;

unsigned int NSobol, continue_skip=0, Nmodes=0;
string input_file, continue_from;

po::variables_map process_options(int argc,  char *  argv[])
{
    
    po::options_description cmd_line_opts("Start Parameters");
    cmd_line_opts.add_options()
    ("help", "Help message")
    ("NSobol,N", po::value<unsigned int>(&NSobol)->default_value(100000), "Length of the Sobol sequence")
    ("continue,c", po::value<string>(), "Name (in format M_%06d.dat) of the file to continue from.")
    ("input-file,f",po::value<string>(), "Input file name from Vladimir")
    ("modes,m", po::value<unsigned int>(&Nmodes)->default_value(0), "Number of modes to include from highest to lowest, inclusively.")
    ("spectrum,s", po::value<string>(), "Spectrum for Hamiltonian saved in M_%06d.dat");
    
    po::positional_options_description pos_opts;
    pos_opts.add("input-file", 1);
    
    po::variables_map vm;
    
    try {
        po::store(po::command_line_parser(argc, argv)
                  .options(cmd_line_opts).positional(pos_opts)
                  .run(), vm);
        po::notify(vm);
    }
    catch (po::error_with_option_name &e) {
        cerr << e.what() << endl;
        cerr << "Finished.\n";
        exit(EXIT_FAILURE);
    }
    
    input_file = vm["input-file"].as<string>();
    if (vm.count("continue")) {
        continue_from = vm["continue"].as<string>();
        int i = continue_from.find_last_of('_')+1;
        int j = continue_from.find('.', i);
        
        istringstream ss(continue_from.substr(i, j-i + 1));
        ss >> continue_skip;
    }
    return vm;
}


vec dsyevr(mat &A, mat *C=NULL)
{
    
    char JOBZ='N';
    
    char RANGE='A';
    char UPLO='U';
    int N = A.n_rows;
    int LDA = A.n_rows;
    double VU = 0.0, VL = 0.0;
    int IL = 0, IU = 0;
    double ABSTOL = dlamch_("Safe minimum");
    int M;
    vec W(A.n_rows);
    double *Z = NULL;
    int LDZ = LDA;
    Col<int> ISUPPZ(2*N);
    int LWORK = 36*N;
    vec WORK(LWORK);
    int LIWORK = 10*N;
    Col<int> IWORK(LIWORK);
    int INFO;
    
    if (C != NULL) {
        JOBZ='V';
        C->set_size(A.n_rows, A.n_cols);
        Z = C->memptr();
    }
    dsyevr_(&JOBZ, &RANGE, &UPLO,  &N,  A.memptr(),  &LDA,  &VL,  &VU,
            &IL,  &IU, &ABSTOL, &M,  W.memptr(),  Z, &LDZ, ISUPPZ.memptr(),
            WORK.memptr(), &LWORK, IWORK.memptr(), &LIWORK, &INFO );
    
    return W;
}



void dump_spectrum(vec &charges, mat &MU, mat &H, MatrixElements& me, 
                   ostream& specout, ostream& dipoleout, string Cout="")
{
    char s[256];
    mat C;
    vec evals = dsyevr(H, &C);

    if (Cout.length() > 0) {
        C.save(Cout, raw_ascii);
    }
    
    mat mu = me.transitionDipole(charges, MU, C);
    for (int ii=1; ii< evals.n_rows; ii++) {
       specout << (evals(ii) - evals(0))*autocm << " ";
       dipoleout << norm(mu.col(ii-1), 2) << " ";
    }
    (specout  << endl).flush();
    (dipoleout << endl).flush();
}


int main (int argc, char *  argv[]) {
    int N;
    mat H, U;
    vec mass, x0, omegasq0;
    po::variables_map vm;
    
    vm = process_options(argc, argv);
    load_from_vladimir(input_file, N, mass, x0, H);
    omegasq0 = dsyevr(H, &U);
    
    ofstream sout("omega0.dat");
    sout << sqrt(abs(omegasq0))*autocm << endl;
    sout.close();
    
    if (Nmodes == 0) Nmodes = 3*N - 6;
    
    vec omega = sqrt(omegasq0.rows(3*N-Nmodes,3*N-1));
    vec alpha = sqrt(omega);
    mat MU = U.cols(3*N-Nmodes,3*N-1);
    

    for (int i=0; i<MU.n_cols; i++) {
        MU.col(i) /= sqrt(mass);
    }
    mat MUT=MU.t();

    int Nstates = 1 + Nmodes*(Nmodes + 3)/2;
    
    mat M(Nstates, Nstates);
    MatrixElements  me(omega);
    
    vec TIP4P_charges(N);
    TIP4P_charges.fill(0.5*qM);
    for (int i=0; i<N; i+= 3) TIP4P_charges[i] = -qM;
    

    if (vm.count("spectrum") >0) {
        string iname = vm["spectrum"].as<string>();
        
        M.load(iname, raw_ascii);
        
        string tname = "freq" + iname.substr(1);
        ofstream specout(tname.c_str());
        ofstream dipoleout(("dipole" + iname.substr(1)).c_str());
        
        dump_spectrum(TIP4P_charges, MU, M, me, specout, dipoleout, "C" + iname.substr(1)); 
        exit(EXIT_SUCCESS);
    }

    int64_t sobol_skip=1<<int(ceil(log2((double)NSobol)));
    if (continue_skip > 0) {
        sobol_skip += continue_skip;
        M.load(continue_from, raw_ascii);
        M *= -1.0;
        me.addEkinDoubles(M);
        M *= -(double)continue_skip;
    }
    
    vec q(Nmodes), Vq(Nmodes), r(3*N), Vr(3*N);
    double V;

    
    int NO = N/3;
    TIP4Pc(NO, x0.memptr(), &V, Vr.memptr());
    cout << "V0 = "<<V<<endl;

    mat Mout;
    ofstream specout("freq.dat");
    ofstream dipoleout("dipole.dat");
    
    
    for (int i=continue_skip; i<continue_skip+NSobol; i++) {
        
        sobol_stdnormal_c(q.n_rows, &sobol_skip, q.memptr());
        q = q / (sqrt(2.0)*sqrt(omega));
        r = MU*q + x0;
        TIP4Pc(NO, r.memptr(), &V, Vr.memptr());
        Vq = -MUT * Vr;
        me.addEpotDoubles(q, V, Vq, M);
        
        if ( (i+1)%10000==0) {
            cout << i+1 << endl;
            Mout = M / (i+1);
            me.addEkinDoubles(Mout);
            
            char s[100];
            sprintf(s, "M_%07d.dat", i+1);
            Mout.save(s, raw_ascii);
            
            dump_spectrum(TIP4P_charges, MU, Mout, me, specout, dipoleout);
        }
    }
    specout.close();
    dipoleout.close();
    M /= NSobol;
}
