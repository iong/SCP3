#include <cmath>
#include <cstdlib>
#include <iostream>

#include "ttm3f.h"
#include "qtip4pf.h"
#include "ttm4-hbb2-x3b.h"

#include "bowman-fortran.h"
#include "bowman.h"

namespace {

const double prism_xyz[] = { // O H H O H H ... O H H
   -1.501133590e+00,  -1.973467100e-01,   1.441262950e+00,
   -2.000985830e+00,  -3.880940000e-01,   2.232286820e+00,
   -6.130182100e-01,  -6.045016400e-01,   1.565741520e+00,
   -1.738921070e+00,  -3.566716500e-01,  -1.311647220e+00,
   -2.485655290e+00,  -7.566635700e-01,  -1.757003500e+00,
   -1.883593590e+00,  -4.958262000e-01,  -3.594679700e-01,
   -5.900491100e-01,   2.045369740e+00,  -1.139971100e-01,
   -1.009913650e+00,   1.627267800e+00,  -8.732886700e-01,
   -9.614419100e-01,   1.549440890e+00,   6.291669800e-01,
    9.999999100e-01,  -1.179514770e+00,   1.454043080e+00,
    1.550433810e+00,  -3.959168300e-01,   1.334139300e+00,
    1.004009200e+00,  -1.558530670e+00,   5.592624600e-01,
    9.829462100e-01,  -1.411538550e+00,  -1.350690300e+00,
    1.478885410e+00,  -5.856945200e-01,  -1.278498130e+00,
    7.624997000e-02,  -1.127406420e+00,  -1.530339410e+00,
    1.967206930e+00,   1.048112740e+00,  -1.378723900e-01,
    2.659992820e+00,   1.698771500e+00,  -2.392144000e-01,
    1.124990980e+00,   1.547180090e+00,  -1.466355100e-01
};

template <typename POT_TYPE>
void call_potential(POT_TYPE& potential, size_t nw, const double* xyz)
{
    double grd[9*nw];

    const double E = potential(nw, xyz, grd);

    double grdsq(0);
    for (size_t i = 0; i < 9*nw; ++i)
        grdsq += grd[i]*grd[i];

    std::cout << "==============================================\n\n"
              << "potential: " << potential.name() << "\n\n"
              << "E = " << E << " kcal/mol\n"
              << "|g| = " << std::sqrt(grdsq) << " kcal/mol/A\n\n";
}

} // namespace

int main(int argc, char** argv)
{
    // initialize WHBB

    if (getenv("BOWMAN_DATADIR") == 0)
        putenv("BOWMAN_DATADIR=./bowman/data");

    try {
        h2o::fortran::pes2b_init();
        h2o::fortran::pes3b_init();
    } catch (const std::exception& e) {
        std::cerr << "\n ** Error ** : " << e.what()
                  << "\n ** Error ** : (try seting BOWMAN_DATADIR variable)\n"
                  << std::endl;
        return EXIT_FAILURE;
    }

    h2o::ttm3f   pot1;
    h2o::qtip4pf pot2;
    h2o::bowman  pot3;

    h2o::ttm4_hbb2_x3b pot4;

    // and 3B-part of the HBB2-pol

    pot4.load("x3b.nc");

    const size_t num_water = (sizeof(prism_xyz)/sizeof(prism_xyz[0]))/9;

    call_potential(pot1, num_water, prism_xyz);
    call_potential(pot2, num_water, prism_xyz);
    call_potential(pot3, num_water, prism_xyz);
    call_potential(pot4, num_water, prism_xyz);
}
