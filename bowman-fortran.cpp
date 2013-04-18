#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cassert>
#include <cstdlib>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <string>
#include <stdexcept>

#include "bowman-fortran.h"

////////////////////////////////////////////////////////////////////////////////

extern "C" {

// 2-body
void prepot_(const char*, size_t);
void calcpot_(double*, const double*, const size_t*);

// 3-body
void pes0_mp_pes0_init_(const char*, const char*, size_t, size_t);
void pes1c_mp_pes1_init_(const char*, size_t);

double pes_x6y3_mp_pes_x6y3_pot_(const size_t*);

} // extern "C"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

static bool pes2b_initialized = false;
static bool pes3b_initialized = false;

#ifdef HAVE_BOWMAN
const char pes2b_file[] = "fit0000.out.pes4X1-T1.090915";
const char pes3b_dir5[] = "WHBB_MP2_3b5";
const char pes3b_dir6[] = "WHBB_MP2_3b6";
#else
void bail_out()
{
    throw std::runtime_error("Bowman PES is not supported in this build");
}
#endif // HAVE_BOWMAN

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o { namespace fortran {

//----------------------------------------------------------------------------//

void pes2b_init()
{
#ifdef HAVE_BOWMAN
    const char* dir = std::getenv("BOWMAN_DATADIR");
    std::string pes2b_full;

    pes2b_full += (dir ? dir : ".");
    pes2b_full += '/';
    pes2b_full += pes2b_file;

    if (pes2b_full.length() > 80)
        throw std::runtime_error("file name for Bowman pes2b data is too long");

    struct stat buf;
    if (stat(pes2b_full.c_str(), &buf))
        throw std::runtime_error("could not stat Bowman pes2b data file");

    prepot_(pes2b_full.c_str(), pes2b_full.length());
    pes2b_initialized = true;
#else
    bail_out();
#endif // HAVE_BOWMAN
}

//----------------------------------------------------------------------------//

bool pes2b_ready()
{
    return pes2b_initialized;
}

//----------------------------------------------------------------------------//

double pes2b(const double* xyz)
{
#ifdef HAVE_BOWMAN
    assert(pes2b_initialized);

    size_t N = 2;
    double E;
    calcpot_(&E, xyz, &N);

    return E;
#else
    bail_out();
    return 0.0; // makes g++ happier
#endif // HAVE_BOWMAN
}

//----------------------------------------------------------------------------//

void pes3b_init()
{
#ifdef HAVE_BOWMAN
    const char* dir = std::getenv("BOWMAN_DATADIR");
    std::string pes3b_full;

    pes3b_full += (dir ? dir : ".");
    pes3b_full += '/';
    pes3b_full += (std::getenv("WHBB_DO_6TH") ? pes3b_dir6 : pes3b_dir5);

    if (pes3b_full.length() > 80)
        throw std::runtime_error("file name for Bowman pes3b data is too long");

    struct stat buf;
    if (stat(pes3b_full.c_str(), &buf))
        throw std::runtime_error("could not stat Bowman pes3b datadir");

    const char name[] = "x6y3";

    pes0_mp_pes0_init_(pes3b_full.c_str(), 0, pes3b_full.length(), 0);
    pes1c_mp_pes1_init_(name, 4);

    pes3b_initialized = true;
#else
    bail_out();
#endif // HAVE_BOWMAN
}

//----------------------------------------------------------------------------//

bool pes3b_ready()
{
    return pes3b_initialized;
}

//----------------------------------------------------------------------------//

double pes3b(const double* xyz)
{
#ifdef HAVE_BOWMAN
    assert(pes3b_initialized);
    assert(sizeof(size_t) == sizeof(void*));

    // this is not portable (ifort-specific)

    size_t ad[6 + 3*2]; // "array descriptor" ... what a wonderful technology indeed

    ad[0] = (size_t) xyz;
    ad[1] = sizeof(double); // element size
    ad[2] = 0; // offset
    ad[3] = 0x01 | 0x02 | 0x04 ; // flags
    ad[4] = 2; // rank
    // ad[5] reserved

    ad[6] = 3; // extent 1
    ad[7] = sizeof(double); // e-e distance
    ad[8] = 0; // start 1

    ad[9]  = 9; // extent 2
    ad[10] = 3*sizeof(double); // e-e distance
    ad[11] = 0; // start 2

    return pes_x6y3_mp_pes_x6y3_pot_(ad);
#else
    bail_out();
    return 0.0; // makes g++ complain less
#endif // HAVE_BOWMAN
}

//----------------------------------------------------------------------------//

}} // namespace h2o::fortran

////////////////////////////////////////////////////////////////////////////////
