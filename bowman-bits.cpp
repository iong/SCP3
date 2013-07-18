#include <cstddef>

#include "macros.h"

#include "bowman-bits.h"
#include "bowman-fortran.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

const double auang = 0.5291772083;
const double aukcal = 627.51;

const double EPS = 1.0e-5;

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace bowman_bits {

//----------------------------------------------------------------------------//

double pes2b(const double* RESTRICT wa,
             const double* RESTRICT wb)
{
    double x2[18];

    for (size_t k = 0; k < 3; ++k) {
        x2[0 + 6*k] = wa[0 + k]/auang; // Oa
        x2[3 + 6*k] = wb[0 + k]/auang; // Ob

        x2[1 + 6*k] = wa[3 + k]/auang; // Ha1
        x2[2 + 6*k] = wa[6 + k]/auang; // Ha2

        x2[4 + 6*k] = wb[3 + k]/auang; // Hb1
        x2[5 + 6*k] = wb[6 + k]/auang; // Hb2
    }

    return aukcal*h2o::fortran::pes2b(x2);
}

//----------------------------------------------------------------------------//

double pes3b(const double* RESTRICT wa,
             const double* RESTRICT wb,
             const double* RESTRICT wc)
{
    double x3[27];

    for (size_t k = 0; k < 3; ++k) {
        x3[0  + k] = wa[3 + k]/auang; // Ha1
        x3[3  + k] = wa[6 + k]/auang; // Ha2
        x3[6  + k] = wb[3 + k]/auang; // Hb1
        x3[9  + k] = wb[6 + k]/auang; // Hb2
        x3[12 + k] = wc[3 + k]/auang; // Hc1
        x3[15 + k] = wc[6 + k]/auang; // Hc2
        x3[18 + k] = wa[0 + k]/auang; // Oa
        x3[21 + k] = wb[0 + k]/auang; // Ob
        x3[24 + k] = wc[0 + k]/auang; // Oc
    }

    return aukcal*h2o::fortran::pes3b(x3);
}

//----------------------------------------------------------------------------//

double pes2b(const double* RESTRICT wa,
             const double* RESTRICT wb,
             double* RESTRICT dwa,
             double* RESTRICT dwb)
{
    double x2[18];

    for (size_t k = 0; k < 3; ++k) {
        x2[0 + 6*k] = wa[0 + k]/auang; // Oa
        x2[3 + 6*k] = wb[0 + k]/auang; // Ob

        x2[1 + 6*k] = wa[3 + k]/auang; // Ha1
        x2[2 + 6*k] = wa[6 + k]/auang; // Ha2

        x2[4 + 6*k] = wb[3 + k]/auang; // Hb1
        x2[5 + 6*k] = wb[6 + k]/auang; // Hb2
    }

    double E = aukcal*h2o::fortran::pes2b(x2);

    // finite differences

    double factor = aukcal/auang/EPS/2;

    double dX[18];
    for (size_t n = 0; n < 18; ++n) {
        double x0 = x2[n];

        x2[n] = x0 + EPS;
        double Ep = h2o::fortran::pes2b(x2);

        x2[n] = x0 - EPS;
        double Em = h2o::fortran::pes2b(x2);

        dX[n] = factor*(Ep - Em);
        x2[n] = x0;
    }

    for (size_t k = 0; k < 3; ++k) {
        dwa[0 + k] = dX[0 + 6*k]; // Oa
        dwb[0 + k] = dX[3 + 6*k]; // Ob

        dwa[3 + k] = dX[1 + 6*k]; // Ha1
        dwa[6 + k] = dX[2 + 6*k]; // Ha2

        dwb[3 + k] = dX[4 + 6*k]; // Hb1
        dwb[6 + k] = dX[5 + 6*k]; // Hb2
    }

    return E;
}

//----------------------------------------------------------------------------//

double pes3b(const double* RESTRICT wa,
             const double* RESTRICT wb,
             const double* RESTRICT wc,
             double* RESTRICT dwa,
             double* RESTRICT dwb,
             double* RESTRICT dwc)
{
    double x3[27];

    for (size_t k = 0; k < 3; ++k) {
        x3[0  + k] = wa[3 + k]/auang; // Ha1
        x3[3  + k] = wa[6 + k]/auang; // Ha2
        x3[6  + k] = wb[3 + k]/auang; // Hb1
        x3[9  + k] = wb[6 + k]/auang; // Hb2
        x3[12 + k] = wc[3 + k]/auang; // Hc1
        x3[15 + k] = wc[6 + k]/auang; // Hc2
        x3[18 + k] = wa[0 + k]/auang; // Oa
        x3[21 + k] = wb[0 + k]/auang; // Ob
        x3[24 + k] = wc[0 + k]/auang; // Oc
    }

    double E = aukcal*h2o::fortran::pes3b(x3);

    // finite differences

    double factor = aukcal/auang/EPS/2;

    double dX[27];
    for (size_t n = 0; n < 27; ++n) {
        double x0 = x3[n];

        x3[n] = x0 + EPS;
        double Ep = h2o::fortran::pes3b(x3);

        x3[n] = x0 - EPS;
        double Em = h2o::fortran::pes3b(x3);

        dX[n] = factor*(Ep - Em);
        x3[n] = x0;
    }

    for (size_t k = 0; k < 3; ++k) {
        dwa[3 + k] = dX[0  + k]; // Ha1
        dwa[6 + k] = dX[3  + k]; // Ha2
        dwb[3 + k] = dX[6  + k]; // Hb1
        dwb[6 + k] = dX[9  + k]; // Hb2
        dwc[3 + k] = dX[12 + k]; // Hc1
        dwc[6 + k] = dX[15 + k]; // Hc2
        dwa[0 + k] = dX[18 + k]; // Oa
        dwb[0 + k] = dX[21 + k]; // Ob
        dwc[0 + k] = dX[24 + k]; // Oc
    }

    return E;
}

//----------------------------------------------------------------------------//

} // namespace bowman_bits

////////////////////////////////////////////////////////////////////////////////
