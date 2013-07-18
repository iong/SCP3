#ifndef __H2O_H__
#define __H2O_H__

#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

#include "qtip4pf.h"

#ifdef HAVE_BOWMAN
#include "bowman.h"
#include "ttm4-hbb2-x3b.h"
#endif

namespace h2o {
    static PES *PESFromString(const string& name)
    {
        PES *pot;
#ifdef HAVE_BOWMAN
        if (name == "whbb") {
            pot = new h2o::bowman();
        }
        else if (name == "hbb2-pol") {
            h2o::ttm4_hbb2_x3b *p = new h2o::ttm4_hbb2_x3b();
            p->load("x3b.nc");
            pot = p;
        }
#else
        if (name == "whbb" || name == "hbb2-pol") {
            cerr << "Support for Bowman's WHBB was not included." << endl;
            exit(EXIT_FAILURE);
        }
#endif

        else if (name == "qtip4pf") {
            pot = new qtip4pf();
        }
        else {
            cerr << "Unknown H2O potential: " << name << endl;
            exit(EXIT_FAILURE);
        }
        
        return pot;
    }
    
}

#endif // __H2O_H__
