#ifndef __H2O_H__
#define __H2O_H__

namespace h2o {
    class Potential {
    public:
        Potential() {}
        
        virtual const char* name() const = 0;
        virtual double operator()(size_t nw, const double*) = 0; // O H H O H H
        virtual double operator()(size_t nw, const double*, double*) = 0; // O H H O H H
    };
}

/*
namespace h2o_i {
#include <armadillo>
    
#include "F90.h"
    
    using namespace arma;
class Potential {
    public:
        Potential() {}
        virtual void operator()(vec& x, double& V, vec& Vx) = 0;
};

class qTIP4Pf : public Potential {
    public:
        qTIP4Pf() {}

        virtual void operator()(vec& x, double& V, vec& Vx)
        {
            int nb_h2o = x.n_rows / 9;
            TIP4P_UF(nb_h2o, x.memptr(), &V, Vx.memptr());
        }
};

class WHBB : public Potential {
    bool    initialized;
    double  eps;
    public:
        WHBB(double eps) :initialized(false) {
            this->eps = eps;
        }

        virtual void operator()(vec& x, double& V, vec& Vx)
        {
            int nb_h2o = x.n_rows / 9;

            if (!initialized) {
                whbb_pes_init(nb_h2o);
                initialized = true;
            }

            int nb_atoms = x.n_rows / 3;
            whbb_fgrad(nb_atoms, x.memptr(), eps, &V, Vx.memptr());
        }
};
}
 */

#endif // __H2O_H__