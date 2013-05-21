#ifndef BOWMAN_BITS_H
#define BOWMAN_BITS_H

namespace bowman_bits {

//----------------------------------------------------------------------------//

double pes2b(const double* OHH1, const double* OHH2);
double pes3b(const double* OHH1, const double* OHH2, const double* OHH3);

//----------------------------------------------------------------------------//

double pes2b(const double* OHH1, const double* OHH2,
             double* dOHH1, double* dOHH2);

double pes3b(const double* OHH1, const double* OHH2, const double* OHH3,
             double* dOHH1, double* dOHH2, double* dOHH3);

//----------------------------------------------------------------------------//

const double r2i = 6.5; // A
const double r2f = 7.5; // A
const double r3i = 7.0; // A
const double r3f = 8.0; // A

//----------------------------------------------------------------------------//

inline double f_switch(const double& r, const double& ri, const double& rf)
{
    const double ra = (r - ri)/(rf - ri);
    const double ra3 = ra*ra*ra;

    return ra3*(10.0 + ra*(6.0*ra - 15.0));
}

//----------------------------------------------------------------------------//

inline double f_switch
    (const double& r, const double& ri, const double& rf, double& g)
{
    const double t1 = 1.0/(rf - ri);
    const double ra = t1*(r - ri);
    const double ra2 = ra*ra;
    const double ra3 = ra2*ra;
    const double t2 = 1.0 - ra;

    g = 30*ra2*t2*t2*t1;

    return ra3*(10.0 + ra*(6.0*ra - 15.0));
}

//----------------------------------------------------------------------------//

} // namespace bowman_bits

#endif // BOWMAN_BITS_H
