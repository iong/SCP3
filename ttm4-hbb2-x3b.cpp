#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <algorithm>

#include "ps.h"
#include "macros.h"
#include "bowman-bits.h"
#include "bowman-fortran.h"
#include "ttm4-hbb2-x3b.h"

#define dont_DISABLE_HBB2 yes
#define dont_DISABLE_X3B yes

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

inline double dist(const double* RESTRICT x1,
                   const double* RESTRICT x2)
{
    const double d1 = x1[0] - x2[0];
    const double d2 = x1[1] - x2[1];
    const double d3 = x1[2] - x2[2];

    return std::sqrt(d1*d1 + d2*d2 + d3*d3);
}

//----------------------------------------------------------------------------//

const double r2i = 5.5;
const double r2f = 7.5;

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o {

//----------------------------------------------------------------------------//

ttm4_hbb2_x3b::ttm4_hbb2_x3b()
{
}
    
ttm4_hbb2_x3b::~ttm4_hbb2_x3b()
{
}

//----------------------------------------------------------------------------//

double ttm4_hbb2_x3b::operator()
    (size_t nw, const double* RESTRICT crd)
{
    assert(crd || nw == 0);

    double E1b(0);
    for (size_t i = 0; i < nw; ++i)
        E1b += ps::pot_nasa(crd + 9*i, 0);

    // 2-body = (1-s)*[E_HBB2 - E_ind] + s*E_elec ( + s*E_vdw, but it is zero)

    double E2b(0), U2b_ind(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;
            const size_t j9 = 3*j3;

            // O-O distance
            const double dRij = dist(crd + i9, crd + j9);

            double dimer[18];
            std::copy(crd + i9, crd + i9 + 9, dimer);
            std::copy(crd + j9, crd + j9 + 9, dimer + 9);

            double Ettm4_elec, Ettm4_ind;
            m_ttm4(2, dimer, Ettm4_elec, 0, Ettm4_ind, 0);

            if (dRij < r2f) {
#               ifdef DISABLE_HBB2
                const double ET1 = 0.0;
#               else
                const double ET1 = bowman_bits::pes2b(crd + i9, crd + j9);
#               endif // DISABLE_HBB2

                if (dRij < r2i) {
                    E2b += ET1;
                    U2b_ind += Ettm4_ind;
                } else {
                    const double s = bowman_bits::f_switch(dRij, r2i, r2f);
                    E2b += (1.0 - s)*ET1 + s*Ettm4_elec;
                    U2b_ind += (1.0 - s)*Ettm4_ind;
                }
            } else {
                E2b += Ettm4_elec;
            }
        }
    }

    // 3-body terms

    double E3b(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i9 = 9*i;
        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j9 = 9*j;
            for (size_t k = j + 1; k < nw; ++k) {
                const size_t k9 = 9*k;
                const double E3 = m_x3b(crd + i9, crd + j9, crd + k9, 0, 0, 0);
                E3b += E3;
            }
        }
    }

    // N-body terms

    double ENb, Ettm4_elec;
    m_ttm4(nw, crd, Ettm4_elec, 0, ENb, 0);

    ENb -= U2b_ind;

    return E1b + E2b + E3b + ENb;
}

//----------------------------------------------------------------------------//

double ttm4_hbb2_x3b::operator()
    (size_t nw, const double* RESTRICT crd,
                      double* RESTRICT grad)
{
    assert((crd && grad) || nw == 0);

    // N-body terms

    double E1b, ENb(0);
    m_ttm4(nw, crd, E1b, 0, ENb, grad);

    // 1-body terms

    E1b = 0.0;
    for (size_t i = 0; i < nw; ++i) {
        const size_t i9 = 9*i;
        double ps_grad[9];
        E1b += ps::pot_nasa(crd + i9, ps_grad);
        for (size_t k = 0; k < 9; ++k)
            grad[i9 + k] += ps_grad[k];
    }

    // 2-body = (1-s)*[E_HBB2 - E_ind] + s*E_elec

    double E2b(0), U2b_ind(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;
            const size_t j9 = 3*j3;

            // O-O distance
            double Rij[3], dRij(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[i9 + k] - crd[j9 + k];
                dRij += Rij[k]*Rij[k];
            }

            dRij = std::sqrt(dRij);

            double dimer[18];
            std::copy(crd + i9, crd + i9 + 9, dimer);
            std::copy(crd + j9, crd + j9 + 9, dimer + 9);

            double Ettm4_elec, gEttm4_elec[18];
            double Ettm4_ind, gEttm4_ind[18];
            m_ttm4(2, dimer, Ettm4_elec, gEttm4_elec, Ettm4_ind, gEttm4_ind);

            if (dRij < r2f) {
                double dpesi[9], dpesj[9];
#               ifdef DISABLE_HBB2
                const double ET1 = 0.0;
                for (size_t k = 0; k < 9; ++k)
                    dpesi[k] = dpesj[k] = 0.0;
#               else
                const double ET1 =
                    bowman_bits::pes2b(crd + i9, crd + j9, dpesi, dpesj);
#               endif // DISABLE_HBB2

                for (size_t k = 0; k < 9; ++k) {
                    dpesi[k] -= gEttm4_ind[k];
                    dpesj[k] -= gEttm4_ind[k + 9];
                }

                if (dRij < r2i) {
                    E2b += ET1;
                    U2b_ind += Ettm4_ind;

                    for (size_t k = 0; k < 9; ++k) {
                        grad[i9 + k] += dpesi[k];
                        grad[j9 + k] += dpesj[k];
                    }
                } else {
                    double g;
                    const double s = bowman_bits::f_switch(dRij, r2i, r2f, g);
                    E2b += (1.0 - s)*ET1 + s*Ettm4_elec;
                    U2b_ind += (1.0 - s)*Ettm4_ind;

                    for (size_t k = 0; k < 9; ++k) {
                        grad[i9 + k] +=
                            (1.0 - s)*dpesi[k] + s*gEttm4_elec[k];
                        grad[j9 + k] +=
                            (1.0 - s)*dpesj[k] + s*gEttm4_elec[k + 9];
                    }

                    g *= (Ettm4_elec + Ettm4_ind - ET1)/dRij;

                    for (size_t k = 0; k < 3; ++k) {
                        grad[i9 + k] += g*Rij[k];
                        grad[j9 + k] -= g*Rij[k];
                    }
                }
            } else {
                E2b += Ettm4_elec;
                for (size_t k = 0; k < 9; ++k) {
                    grad[i9 + k] += gEttm4_elec[k];
                    grad[j9 + k] += gEttm4_elec[k + 9];
                }
            }
        }
    }

    // 3-body terms

    double E3b(0);
    for (size_t i = 0; i < nw; ++i) {
        const size_t i9 = 9*i;
        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j9 = 9*j;
            for (size_t k = j + 1; k < nw; ++k) {
                const size_t k9 = 9*k;

                double d3bi[9], d3bj[9], d3bk[9];

#               ifdef DISABLE_X3B
                const double E3 = 0.0;
                for (size_t a = 0; a < 9; ++a)
                    d3bi[a] = d3bj[a] = d3bk[a] = 0.0;
#               else
                const double E3 = m_x3b(crd + i9, crd + j9, crd + k9,
                                        d3bi, d3bj, d3bk);
#               endif // DISABLE_X3B

                for (size_t a = 0; a < 9; ++a) {
                     grad[i9 + a] += d3bi[a];
                     grad[j9 + a] += d3bj[a];
                     grad[k9 + a] += d3bk[a];
                 }

                E3b += E3;
            }
        }
    }

    ENb -= U2b_ind;

    return E1b + E2b + E3b + ENb;
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
