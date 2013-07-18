//
//  SCP1.h
//  SCP3
//
//  Created by Ionut Georgescu on 7/15/13.
//
//

#ifndef __SCP3__SCP1__
#define __SCP3__SCP1__

#include <iostream>
#include <fstream>
#include <armadillo>


using namespace arma;
using namespace std;

#include "Constants.h"
#include "h2o.h"

class SCP1 {
    double omega_cutoff;
    double epsrel;
    vec mass;
    string PES_name;
    long long sobol_skip;
    long long NSobol;
    double qeps;
    
    ofstream omega_out;
    ofstream free_energy_out;

    vec disp_disp_corr_diag(const vec& omega, double kT);
    double average(const vec& q, const vec& d, const mat& U, vec& UX, mat& UX_q);
    h2o::PES * getPES()
    {
        return h2o::PESFromString(PES_name);
    }
public:
    SCP1(const vec& mass_, const string& PES_name_, long long NSobol_) :
        mass(mass_),
        PES_name(PES_name_),
        NSobol(NSobol_)
    {
        omega_cutoff = 50.0 / autocm;
        epsrel = 1e-4;
        qeps = 2;
        
        sobol_skip = 1 << ((int)floor(log2((float)NSobol)) + 1);
        //sobol_skip = 1000;
        
        omega_out.open("omega_scp1.dat");
        free_energy_out.open("free_energy.dat");
    }
    
    SCP1& setSobolSkip(long long s)
    {
        sobol_skip = s;
        
        return *this;
    }
    
    SCP1& setEPSREL(double e)
    {
        epsrel = e;
        
        return *this;
    }
    
    SCP1& setQEps(double e)
    {
        qeps = e;
        
        return *this;
    }
    
    SCP1& setCutoff(double w)
    {
        omega_cutoff = w;
        
        return *this;
    }
    
    
    
    double operator()(vec& q, double kT, mat& Ks);
};

#endif /* defined(__SCP3__SCP1__) */
