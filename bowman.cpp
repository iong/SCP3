#include <cmath>
#include <cassert>

#include <algorithm>

#include "macros.h"

#include "ps.h"

#include "bowman.h"
#include "bowman-bits.h"
#include "bowman-fortran.h"

#include "ttm3f-bits.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

inline void compute_M_site_crd
    (const double* RESTRICT O,
     const double* RESTRICT H1,
     const double* RESTRICT H2,
     double* RESTRICT M)
{
    using ttm3f_bits::gamma1;
    using ttm3f_bits::gamma2;

    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o {

//----------------------------------------------------------------------------//

bowman::bowman()
{
    m_nw = 0;
    m_mem = 0;
}

//----------------------------------------------------------------------------//

bowman::~bowman()
{
    if (m_nw > 0)
        delete[] m_mem;
}

//----------------------------------------------------------------------------//

void bowman::allocate(size_t nw)
{
    if (m_nw > nw)
        return;

    delete[] m_mem;

    // operator() below
    m_mem = new double[
        3*nw  // m-site
      + 3*nw  // charge
      + 27*nw // grdq
      + 9*nw*(2 + nw) // scratch for the TTM3F bits
    ];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

double bowman::operator()(size_t nw, const double* RESTRICT crd)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers
    double* msite = m_mem;
    double* charge = msite + 3*nw;
    double* scratch = charge + 3*nw;

    // compute M-sites, charges, 1-body terms
    double E1b(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        compute_M_site_crd(crd + i9, crd + i9 + 3, crd + i9 + 6, msite + i3);

        E1b += ps::pot_nasa(crd + i9, 0);

        using ttm3f_bits::gammaM;

        using ttm3f_bits::dms_param1;
        using ttm3f_bits::dms_param2;
        using ttm3f_bits::dms_param3;

        double q3[3];
        ps::dms_nasa(dms_param1, dms_param2, dms_param3,
                     crd + i9, q3, 0, true);

        double tmp = 0.5*gammaM/(1.0 - gammaM);

        using ttm3f_bits::CHARGECON;

        charge[i3 + 0] = CHARGECON*q3[0]/(1.0 - gammaM);          // M
        charge[i3 + 1] = CHARGECON*(q3[1] + tmp*(q3[1] + q3[2])); // H1
        charge[i3 + 2] = CHARGECON*(q3[2] + tmp*(q3[1] + q3[2])); // H2
    }

    // 2-body terms
    double E2b(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        for (size_t j = i + 1; j < nw; ++j) {
            size_t j3 = 3*j;
            size_t j9 = 3*j3;

            // O-O distance
            double dRijsq = 0.0, Rij[3];
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[i9 + k] - crd[j9 + k];
                dRijsq += Rij[k]*Rij[k];
            }

            double dRij = std::sqrt(dRijsq);

            double Ettm3_vdw(0), Ettm3_CC(0), Ettm3_CD_DD(0);

            using bowman_bits::r2i;
            using bowman_bits::r2f;

            ttm3f_bits::U_coul_CC_CD_DD_2w
                (msite + i3, crd + i9 + 3, charge + i3,
                 msite + j3, crd + j9 + 3, charge + j3,
                 Ettm3_CC, Ettm3_CD_DD);

            E2b += (nw - 3.0)*Ettm3_CD_DD;

            if (dRij >= r2i) {
                using ttm3f_bits::vdwC;
                using ttm3f_bits::vdwD;
                using ttm3f_bits::vdwE;

                double dR6 = dRijsq*dRijsq*dRijsq;
                double expon = vdwD*std::exp(-vdwE*dRij);

                Ettm3_vdw = vdwC/dR6 + expon;
            } // dRij >= r2i

            double Ettm3 = Ettm3_vdw + Ettm3_CC + Ettm3_CD_DD;

            if (dRij < r2f) {
                double ET1 = bowman_bits::pes2b(crd + i9, crd + j9);

                if (dRij < r2i) {
                    E2b += ET1;
                } else {
                    double s = bowman_bits::f_switch(dRij, r2i, r2f);
                    E2b += (1.0 - s)*ET1 + s*Ettm3;
                }
            } else {
                E2b += Ettm3;
            }
        }
    }

    // 3-body contributions

    double E3b(0);
    for (size_t a = 0; a < nw; ++a) {
        size_t a9 = 9*a;

        for (size_t b = a + 1; b < nw; ++b) {
            size_t b9 = 9*b;

            double Rab[3], dRab(0);
            for (size_t k = 0; k < 3; ++k) {
                Rab[k] = crd[a9 + k] - crd[b9 + k];
                dRab += Rab[k]*Rab[k];
            }

            dRab = std::sqrt(dRab);

            using bowman_bits::r3i;
            using bowman_bits::r3f;

            if (dRab > r3f)
                continue;

            for (size_t c = b + 1; c < nw; ++c) {
                size_t c9 = 9*c;

                double Rac[3], dRac(0);
                double Rbc[3], dRbc(0);
                for (size_t k = 0; k < 3; ++k) {
                    Rac[k] = crd[a9 + k] - crd[c9 + k];
                    dRac += Rac[k]*Rac[k];
                    Rbc[k] = crd[b9 + k] - crd[c9 + k];
                    dRbc += Rbc[k]*Rbc[k];
                }

                dRac = std::sqrt(dRac);
                dRbc = std::sqrt(dRbc);

                double rmax = std::max(dRab, std::max(dRac, dRbc));

                if (rmax < r3f) {
                    double E3 =
                        bowman_bits::pes3b(crd + a9, crd + b9, crd + c9);

                    if (rmax < r3i) {
                        E3b += E3;
                    } else {
                        double s = bowman_bits::f_switch(rmax, r3i, r3f);
                        E3b += E3*(1.0 - s);
                    }
                }
            }
        }
    }

    // N-body terms (except 2-body U_coul_CD_DD, to avoid recomputing them)

    double ENb(0);

    ttm3f_bits::U_coul_CD_DD(nw, crd, msite, charge, scratch, ENb);

    for (size_t a = 0; a < nw; ++a) {
        size_t a3 = 3*a;

        for (size_t b = a + 1; b < nw; ++b) {
            size_t b3 = 3*b;

            for (size_t c = b + 1; c < nw; ++c) {
                size_t c3 = 3*c;

                double Eind;
                ttm3f_bits::U_coul_CD_DD_3w
                    (msite + a3, crd + 3*a3 + 3, charge + a3,
                     msite + b3, crd + 3*b3 + 3, charge + b3,
                     msite + c3, crd + 3*c3 + 3, charge + c3,
                     Eind);

                ENb -= Eind;
            }
        }
    }

    return E1b + E2b + E3b + ENb;
}

//----------------------------------------------------------------------------//

double bowman::operator()
    (size_t nw, const double* RESTRICT crd, double* RESTRICT grad)
{
    assert((crd && grad) || nw == 0);

    allocate(nw);

    // setup pointers
    double* msite = m_mem;
    double* charge = msite + 3*nw;
    double* grdq = charge + 3*nw;
    double* scratch = grdq + 27*nw;

    // M-site, charges, 1-body terms
    double E1b(0);
    for (size_t n = 0; n < nw; ++n) {
        size_t n3 = 3*n;
        size_t n9 = 3*n3;

        compute_M_site_crd(crd + n9, crd + n9 + 3, crd + n9 + 6, msite + n3);

        E1b += ps::pot_nasa(crd + n9, grad + n9);

        using ttm3f_bits::gamma1;
        using ttm3f_bits::gamma2;

        using ttm3f_bits::dms_param1;
        using ttm3f_bits::dms_param2;
        using ttm3f_bits::dms_param3;

        double q3[3], dq3[27];
        ps::dms_nasa(dms_param1, dms_param2, dms_param3,
                     crd + n9, q3, dq3, true);

        double tmp = gamma2/gamma1;

        charge[n3 + 0] = q3[0]/gamma1;                // M
        charge[n3 + 1] = q3[1] + tmp*(q3[1] + q3[2]); // H1
        charge[n3 + 2] = q3[2] + tmp*(q3[1] + q3[2]); // H2

#define DQ3(i,j,k) dq3[k + 3*(j + 3*i)]
#define GRDQ(i,j,k) grdq[k + 3*(j + 3*(i + 3*n))]

        for (size_t k = 0; k < 3; ++k) {
            GRDQ(0, 0, k) = DQ3(0, 0, k) + tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 0, k) = DQ3(1, 0, k) + tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 0, k) = DQ3(2, 0, k) + tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));

            GRDQ(0, 1, k) = DQ3(0, 1, k) + tmp*(DQ3(0, 1, k) + DQ3(0, 0, k));
            GRDQ(1, 1, k) = DQ3(1, 1, k) + tmp*(DQ3(1, 1, k) + DQ3(1, 0, k));
            GRDQ(2, 1, k) = DQ3(2, 1, k) + tmp*(DQ3(2, 1, k) + DQ3(2, 0, k));

            GRDQ(0, 2, k) = DQ3(0, 2, k) - 2*tmp*(DQ3(0, 0, k) + DQ3(0, 1, k));
            GRDQ(1, 2, k) = DQ3(1, 2, k) - 2*tmp*(DQ3(1, 0, k) + DQ3(1, 1, k));
            GRDQ(2, 2, k) = DQ3(2, 2, k) - 2*tmp*(DQ3(2, 0, k) + DQ3(2, 1, k));
        }

        using ttm3f_bits::CHARGECON;

        for (size_t k = 0; k < 3; ++k)
            charge[n3 + k] *= CHARGECON;

        for (size_t k = 0; k < 27; ++k)
            grdq[27*n + k] *= CHARGECON;
    }

    // 2-body terms
    double E2b(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        for (size_t j = i + 1; j < nw; ++j) {
            size_t j3 = 3*j;
            size_t j9 = 3*j3;

            // O-O distance
            double dRijsq = 0.0, Rij[3];
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[i9 + k] - crd[j9 + k];
                dRijsq += Rij[k]*Rij[k];
            }

            double dRij = std::sqrt(dRijsq);

            double Ettm3_vdw(0), Ettm3_CC(0), Ettm3_CD_DD(0);
            double dttm3i[9], dttm3j[9], dtmpi[9], dtmpj[9];

            ttm3f_bits::U_coul_CC_CD_DD_2w
                (msite + i3, crd + i9 + 3, charge + i3,
                 msite + j3, crd + j9 + 3, charge + j3,
                 Ettm3_CC, Ettm3_CD_DD,
                 grdq + 27*i, grdq + 27*j,
                 dttm3i, dtmpi, dttm3j, dtmpj);

            double w_CD_DD = nw - 3.0;

            E2b += w_CD_DD*Ettm3_CD_DD;
            for (size_t k = 0; k < 9; ++k) {
                grad[i9 + k] += w_CD_DD*dtmpi[k];
                grad[j9 + k] += w_CD_DD*dtmpj[k];
            }

            using bowman_bits::r2i;
            using bowman_bits::r2f;

            if (dRij >= r2i) {
                using ttm3f_bits::vdwC;
                using ttm3f_bits::vdwD;
                using ttm3f_bits::vdwE;

                double dR6 = dRijsq*dRijsq*dRijsq;
                double expon = vdwD*std::exp(-vdwE*dRij);

                Ettm3_vdw = vdwC/dR6 + expon;

                double tmp = -(6.0*vdwC/dR6)/dRijsq - vdwE*expon/dRij;

                // assemble Ettm3 gradients in dttm3?
                for (size_t k = 0; k < 3; ++k) {
                    dttm3i[k] += tmp*Rij[k];
                    dttm3j[k] -= tmp*Rij[k];
                }
                for (size_t k = 0; k < 9; ++k) {
                    dttm3i[k] += dtmpi[k];
                    dttm3j[k] += dtmpj[k];
                }
            } // dRij >= r2i

            double Ettm3 = Ettm3_vdw + Ettm3_CC + Ettm3_CD_DD;

            if (dRij < r2f) {
                double dpesi[9], dpesj[9];
                double ET1 =
                    bowman_bits::pes2b(crd + i9, crd + j9, dpesi, dpesj);

                if (dRij < r2i) {
                    E2b += ET1;

                    for (size_t k = 0; k < 9; ++k) {
                        grad[i9 + k] += dpesi[k];
                        grad[j9 + k] += dpesj[k];
                    }
                } else {
                    double g; // ds/dr
                    double s = bowman_bits::f_switch(dRij, r2i, r2f, g);

                    g *= (Ettm3 - ET1)/dRij;

                    E2b += (1.0 - s)*ET1 + s*Ettm3;
                    for (size_t k = 0; k < 3; ++k) {
                        grad[i9 + k] += (1.0 - s)*dpesi[k] + s*dttm3i[k]
                                      + g*Rij[k];
                        grad[j9 + k] += (1.0 - s)*dpesj[k] + s*dttm3j[k]
                                      - g*Rij[k];
                    }
                    for (size_t k = 3; k < 9; ++k) {
                        grad[i9 + k] += (1.0 - s)*dpesi[k] + s*dttm3i[k];
                        grad[j9 + k] += (1.0 - s)*dpesj[k] + s*dttm3j[k];
                    }
                }
            } else {
                E2b += Ettm3;
                for (size_t k = 0; k < 9; ++k) {
                    grad[i9 + k] += dttm3i[k];
                    grad[j9 + k] += dttm3j[k];
                }
            }
        }
    }

    // 3-body contributions

    double E3b(0);
    for (size_t a = 0; a < nw; ++a) {
        size_t a9 = 9*a;

        for (size_t b = a + 1; b < nw; ++b) {
            size_t b9 = 9*b;

            double Rab[3], dRab(0);
            for (size_t k = 0; k < 3; ++k) {
                Rab[k] = crd[a9 + k] - crd[b9 + k];
                dRab += Rab[k]*Rab[k];
            }

            dRab = std::sqrt(dRab);

            using bowman_bits::r3i;
            using bowman_bits::r3f;

            if (dRab > r3f)
                continue;

            for (size_t c = b + 1; c < nw; ++c) {
                size_t c9 = 9*c;

                double Rac[3], dRac(0);
                double Rbc[3], dRbc(0);
                for (size_t k = 0; k < 3; ++k) {
                    Rac[k] = crd[a9 + k] - crd[c9 + k];
                    dRac += Rac[k]*Rac[k];
                    Rbc[k] = crd[b9 + k] - crd[c9 + k];
                    dRbc += Rbc[k]*Rbc[k];
                }

                dRac = std::sqrt(dRac);
                dRbc = std::sqrt(dRbc);

                double rmax = std::max(dRab, std::max(dRac, dRbc));

                if (rmax < r3f) {
                    double dpesa[9], dpesb[9], dpesc[9];
                    double E3 =
                        bowman_bits::pes3b(crd + a9, crd + b9, crd + c9,
                                           dpesa, dpesb, dpesc);

                    if (rmax < r3i) {
                        E3b += E3;
                        for (size_t k = 0; k < 9; ++k) {
                            grad[a9 + k] += dpesa[k];
                            grad[b9 + k] += dpesb[k];
                            grad[c9 + k] += dpesc[k];
                        }
                    } else {
                        double g;
                        double s =
                            bowman_bits::f_switch(rmax, r3i, r3f, g);

                        E3b += E3*(1.0 - s);
                        for (size_t k = 0; k < 9; ++k) {
                            grad[a9 + k] += (1.0 - s)*dpesa[k];
                            grad[b9 + k] += (1.0 - s)*dpesb[k];
                            grad[c9 + k] += (1.0 - s)*dpesc[k];
                        }

                        g *= E3;

                        if (dRab > dRac && dRab > dRbc) {
                            g /= dRab;
                            for (size_t k = 0; k < 3; ++k) {
                                grad[a9 + k] -= g*Rab[k];
                                grad[b9 + k] += g*Rab[k];
                            }
                        }

                        if (dRac > dRab && dRac > dRbc) {
                            g /= dRac;
                            for (size_t k = 0; k < 3; ++k) {
                                grad[a9 + k] -= g*Rac[k];
                                grad[c9 + k] += g*Rac[k];
                            }
                        }

                        if (dRbc > dRab && dRbc > dRac) {
                            g /= dRbc;
                            for (size_t k = 0; k < 3; ++k) {
                                grad[b9 + k] -= g*Rbc[k];
                                grad[c9 + k] += g*Rbc[k];
                            }
                        }
                    }
                }
            }
        }
    }

    // N-body terms (except 2-body U_coul_CD_DD, to avoid recomputing them)

    double ENb(0);

    ttm3f_bits::U_coul_CD_DD(nw, crd, msite, charge, scratch, ENb, grdq, grad);

    for (size_t a = 0; a < nw; ++a) {
        size_t a3 = 3*a;
        size_t a9 = 3*a3;

        for (size_t b = a + 1; b < nw; ++b) {
            size_t b3 = 3*b;
            size_t b9 = 3*b3;

            for (size_t c = b + 1; c < nw; ++c) {
                size_t c3 = 3*c;
                size_t c9 = 3*c3;

                double crd3[27], msite3[9], charge3[9];
                double Eind, grdq3[3*27], grad3[27];

                for (size_t k = 0; k < 3; ++k) {
                    msite3[0 + k] = msite[a3 + k];
                    msite3[3 + k] = msite[b3 + k];
                    msite3[6 + k] = msite[c3 + k];

                    charge3[0 + k] = charge[a3 + k];
                    charge3[3 + k] = charge[b3 + k];
                    charge3[6 + k] = charge[c3 + k];
                }

                for (size_t k = 0; k < 9; ++k) {
                    crd3[0  + k] = crd[a9 + k];
                    crd3[9  + k] = crd[b9 + k];
                    crd3[18 + k] = crd[c9 + k];
                }

                for (size_t k = 0; k < 27; ++k) {
                    grdq3[0*27 + k] = grdq[3*a9 + k];
                    grdq3[1*27 + k] = grdq[3*b9 + k];
                    grdq3[2*27 + k] = grdq[3*c9 + k];
                    grad3[k] = 0.0;
                }

                ttm3f_bits::U_coul_CD_DD(3, crd3, msite3, charge3,
                                         scratch, Eind, grdq3, grad3);

                ENb -= Eind;

                for (size_t k = 0; k < 9; ++k) {
                    grad[a9 + k] -= grad3[0  + k];
                    grad[b9 + k] -= grad3[9  + k];
                    grad[c9 + k] -= grad3[18 + k];
                }
            }
        }
    }

    return E1b + E2b + E3b + ENb;
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
