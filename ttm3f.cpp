#include <cmath>
#include <cassert>

#include <algorithm>

#include "ps.h"

#include "ttm3f.h"
#include "ttm3f-bits.h"
#include "ttm3f-bits-smear.h"

#include "macros.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

const double polarM = 1.444;

const double polfacH = 0.496;
const double polfacM = 0.837;

const double AA_HH = std::pow(polfacH*polfacH, 1.0/6.0);
const double AA_HM = std::pow(polfacH*polfacM, 1.0/6.0);
const double AA_MM = std::pow(polfacM*polfacM, 1.0/6.0);

//----------------------------------------------------------------------------//

const size_t dip_maxiter = 500;
const double dip_tolerance = 1.0e-15;

const double DEBYE = 4.8033324;

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

ttm3f::ttm3f()
{
    m_nw = 0;
    m_mem = 0;
}

//----------------------------------------------------------------------------//

ttm3f::~ttm3f()
{
    if (m_nw > 0)
        delete[] m_mem;
}

//----------------------------------------------------------------------------//

void ttm3f::allocate(size_t nw)
{
    if (m_nw > nw)
        return;

    delete[] m_mem;

    // see operator() below
    m_mem = new double[51*nw + 9*nw*nw];

    m_nw = nw;
}

//----------------------------------------------------------------------------//

double ttm3f::operator()(size_t nw, const double* RESTRICT crd)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers (16*nw + 9*nw*nw)
    double* msite = m_mem;
    double* charge = msite + 3*nw;
    double* dip = charge + 3*nw;
    double* dip_prev = dip + 3*nw;
    double* Efd = dip_prev + 3*nw;
    double* Efq = Efd + 3*nw;
    double* phi = Efq + 3*nw;
    double* ddt = phi + 3*nw;

    // zero out Efd/Efq/phi
    std::fill(Efq, Efq + 9*nw, 0.0);

    // compute M-sites, charges, 1-body terms
    double Eint(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        compute_M_site_crd(crd + i9, crd + i9 + 3, crd + i9 + 6, msite + i3);

        Eint += ps::pot_nasa(crd + i9, 0);

        using ttm3f_bits::gamma1;
        using ttm3f_bits::gamma2;

        using ttm3f_bits::dms_param1;
        using ttm3f_bits::dms_param2;
        using ttm3f_bits::dms_param3;

        double q3[3];
        ps::dms_nasa(dms_param1, dms_param2, dms_param3,
                     crd + i9, q3, 0, true);

        double tmp = gamma2/gamma1;

        using ttm3f_bits::CHARGECON;

        charge[i3 + 0] = CHARGECON*q3[0]/gamma1;                  // M
        charge[i3 + 1] = CHARGECON*(q3[1] + tmp*(q3[1] + q3[2])); // H1
        charge[i3 + 2] = CHARGECON*(q3[2] + tmp*(q3[1] + q3[2])); // H2
    }

    // 2-body loop: vdw terms and C-C electrostatics

    double Evdw(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        // Hydrogens - M-sites
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            size_t j3 = 3*j;

            for (size_t l = 1; l < 3; ++l) {
                size_t j3l = j3 + l;
                size_t jh = 3*j3l;

                double Rij[3], Rsq(0);
                for (size_t k = 0; k < 3; ++k) {
                    Rij[k] = msite[i3 + k] - crd[jh + k];
                    Rsq += Rij[k]*Rij[k];
                }

                double ts0, ts1;
                ttm3f_bits::smear01(std::sqrt(Rsq), AA_HM, ts0, ts1);

                phi[i3] += ts0*charge[j3l];
                phi[j3l] += ts0*charge[i3];

                for (size_t k = 0; k < 3; ++k)
                    Efq[i3 + k] += charge[j3l]*ts1*Rij[k];
            }
        }

        // M-sites, Oxygens and Hydrogen - Hydrogen

        for (size_t j = i + 1; j < nw; ++j) {
            size_t j3 = 3*j;
            size_t j9 = 3*j3;

            // H-H
            for (size_t a = 1; a < 3; ++a) {
                size_t i3a = i3 + a;

                for (size_t b = 1; b < 3; ++b) {
                    size_t j3b = j3 + b;

                    double Rsq(0);
                    for (size_t k = 0; k < 3; ++k) {
                        double dx = crd[3*i3a + k] - crd[3*j3b + k];
                        Rsq += dx*dx;
                    }

                    double ts0;
                    ttm3f_bits::smear0(std::sqrt(Rsq), AA_HH, ts0);

                    phi[i3a] += charge[j3b]*ts0;
                    phi[j3b] += charge[i3a]*ts0;
                }
            }

            // O-O distance
            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[i9 + k] - crd[j9 + k];
                Rsq += Rij[k]*Rij[k];
            }

            double dRij = std::sqrt(Rsq);

            using ttm3f_bits::vdwC;
            using ttm3f_bits::vdwD;
            using ttm3f_bits::vdwE;

            double dR6 = Rsq*Rsq*Rsq;
            double expon = vdwD*std::exp(-vdwE*dRij);

            Evdw += vdwC/dR6 + expon;

            // M-M distance
            Rsq = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = msite[i3 + k] - msite[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            double ts0, ts1, ts2;
            ttm3f_bits::smear012(std::sqrt(Rsq), AA_MM, ts0, ts1, ts2);

            phi[i3] += ts0*charge[j3];
            phi[j3] += ts0*charge[i3];

            for (size_t k = 0; k < 3; ++k) {
                Efq[i3 + k] += charge[j3]*ts1*Rij[k];
                Efq[j3 + k] -= charge[i3]*ts1*Rij[k];
            }

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    ddt[3*nw*(i3 + k) + j3 + l] = dd3[k][l];
                    ddt[3*nw*(j3 + k) + i3 + l] = dd3[k][l];
                }
            }
        }

        // diagonal part of ddt
        for (size_t k = 0; k < 3; ++k)
            for (size_t l = 0; l < 3; ++l)
                ddt[3*nw*(i3 + k) + i3 + l] = 0.0;

    }

    // calculate induced dipoles iteratively

    for (size_t i = 0; i < 3*nw; ++i) {
        dip[i] = polarM*Efq[i];
        dip_prev[i] = dip[i];
    }

    double dmix = 0.7;
    double stath = DEBYE/ttm3f_bits::CHARGECON/(std::sqrt(double(nw)));

#ifndef NDEBUG
    bool converged = false;
#endif // NDEBUG

    for (size_t iter = 0; iter < dip_maxiter; ++iter) {

        for (size_t k = 0; k < 3*nw; ++k)
            for (size_t l = 0; l < 3*nw; ++l)
                Efd[k] += ddt[3*nw*k + l]*dip[l];

        double deltadip = 0.0;

        for (size_t i = 0; i < 3*nw; ++i) {
            dip[i] = polarM*(Efq[i] + Efd[i]);
            dip[i] = dmix*dip[i] + (1.0 - dmix)*dip_prev[i];

            double delta = dip[i] - dip_prev[i];
            deltadip += delta*delta;
        }

        deltadip = std::sqrt(deltadip)*stath;

        if (deltadip < dip_tolerance) {
#ifndef NDEBUG
            converged = true;
#endif // NDEBUG
            break; // converged!
        } else {
            std::copy(dip, dip + 3*nw, dip_prev);
            std::fill(Efd, Efd + 3*nw, 0.0);
        }
    }

    assert(converged);

    double Eelec(0), Eind(0);
    for (size_t i = 0; i < 3*nw; ++i) {
        Eelec += phi[i]*charge[i];
        Eind -= dip[i]*Efq[i];
    }

    Eelec /= 2;
    Eind /= 2;

    return Eint + Evdw + Eelec + Eind;
}

//----------------------------------------------------------------------------//

double ttm3f::operator()
    (size_t nw, const double* RESTRICT crd, double* RESTRICT grad)
{
    assert(crd || nw == 0);

    allocate(nw);

    // setup pointers (51*nw + 9*nw*nw)
    double* msite = m_mem;             //  3*nw => 3*nw
    double* charge = msite + 3*nw;     //  3*nw => 6*nw
    double* dip = charge + 3*nw;       //  3*nw => 9*nw
    double* dip_prev = dip + 3*nw;     //  3*nw => 12*nw
    double* Efd = dip_prev + 3*nw;     //  3*nw => 15*nw
    double* Efq = Efd + 3*nw;          //  3*nw => 18*nw
    double* phi = Efq + 3*nw;          //  3*nw => 21*nw
    double* dM = phi + 3*nw;           //  3*nw => 24*nw
    double* ddt = dM + 3*nw;           // (3*nw)*(3*nw) => 24*nw + 9*nw*nw
    double* grdq = ddt + 9*nw*nw;      // 27*nw => 51*nw + 9*nw*nw

    // zero out Efd/Efq/phi/dM
    std::fill(Efd, Efd + 12*nw, 0.0);

    // compute M-sites, charges, 1-body terms
    double Eint(0);
    for (size_t n = 0; n < nw; ++n) {
        size_t n3 = 3*n;
        size_t n9 = 3*n3;

        compute_M_site_crd(crd + n9, crd + n9 + 3, crd + n9 + 6, msite + n3);

        Eint += ps::pot_nasa(crd + n9,  grad + n9);

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

    // 2-body loop: vdw terms and C-C electrostatics

    double Evdw(0);
    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        // Hydrogens - M-sites
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            size_t j3 = 3*j;

            for (size_t l = 1; l < 3; ++l) {
                size_t j3l = j3 + l;
                size_t jh = 3*j3l;

                double Rij[3], Rsq(0);
                for (size_t k = 0; k < 3; ++k) {
                    Rij[k] = msite[i3 + k] - crd[jh + k];
                    Rsq += Rij[k]*Rij[k];
                }

                double ts0, ts1;
                ttm3f_bits::smear01(std::sqrt(Rsq), AA_HM, ts0, ts1);

                phi[i3] += ts0*charge[j3l];
                phi[j3l] += ts0*charge[i3];

                double tmp1 = charge[j3l]*ts1;
                double tmp2 = charge[i3]*tmp1;
                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += tmp1*Rij[k];
                    grad[jh + k] += tmp2*Rij[k];
                }
            }
        }

        // M-sites, Oxygens and Hydrogen - Hydrogen

        for (size_t j = i + 1; j < nw; ++j) {
            size_t j3 = 3*j;
            size_t j9 = 3*j3;

            // H-H
            for (size_t a = 1; a < 3; ++a) {
                size_t i3a = i3 + a;

                for (size_t b = 1; b < 3; ++b) {
                    size_t j3b = j3 + b;

                    double Rij[3], Rsq(0);
                    for (size_t k = 0; k < 3; ++k) {
                        Rij[k] = crd[3*i3a + k] - crd[3*j3b + k];
                        Rsq += Rij[k]*Rij[k];
                    }

                    double ts0, ts1;
                    ttm3f_bits::smear01(std::sqrt(Rsq), AA_HH, ts0, ts1);

                    phi[i3a] += charge[j3b]*ts0;
                    phi[j3b] += charge[i3a]*ts0;

                    // we do not store Efq on hydrogens, therefore
                    double qq = charge[i3a]*charge[j3b]*ts1;
                    for (size_t k = 0; k < 3; ++k) {
                        grad[3*i3a + k] -= qq*Rij[k];
                        grad[3*j3b + k] += qq*Rij[k];
                    }
                }
            }

            // O-O distance
            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = crd[i9 + k] - crd[j9 + k];
                Rsq += Rij[k]*Rij[k];
            }

            double dRij = std::sqrt(Rsq);

            using ttm3f_bits::vdwC;
            using ttm3f_bits::vdwD;
            using ttm3f_bits::vdwE;

            double dR6 = Rsq*Rsq*Rsq;
            double expon = vdwD*std::exp(-vdwE*dRij);

            Evdw += vdwC/dR6 + expon;

            double tmp = -(6.0*vdwC/dR6)/Rsq - vdwE*expon/dRij;
            for (size_t k = 0; k < 3; ++k) {
                grad[i9 + k] += tmp*Rij[k];
                grad[j9 + k] -= tmp*Rij[k];
            }

            // M-M distance
            Rsq = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = msite[i3 + k] - msite[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            double ts0, ts1, ts2;
            ttm3f_bits::smear012(std::sqrt(Rsq), AA_MM, ts0, ts1, ts2);

            phi[i3] += ts0*charge[j3];
            phi[j3] += ts0*charge[i3];

            for (size_t k = 0; k < 3; ++k) {
                Efq[i3 + k] += charge[j3]*ts1*Rij[k];
                Efq[j3 + k] -= charge[i3]*ts1*Rij[k];
            }

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    ddt[3*nw*(i3 + k) + j3 + l] = dd3[k][l];
                    ddt[3*nw*(j3 + k) + i3 + l] = dd3[k][l];
                }
            }
        }

        // diagonal part of ddt
        for (size_t k = 0; k < 3; ++k)
            for (size_t l = 0; l < 3; ++l)
                ddt[3*nw*(i3 + k) + i3 + l] = 0.0;
    }

    // calculate induced dipoles iteratively

    for (size_t i = 0; i < 3*nw; ++i) {
        dip[i] = polarM*Efq[i];
        dip_prev[i] = dip[i];
    }

    double dmix = 0.7;
    double stath = DEBYE/ttm3f_bits::CHARGECON/(std::sqrt(double(nw)));

#ifndef NDEBUG
    bool converged = false;
#endif // NDEBUG

    for (size_t iter = 0; iter < dip_maxiter; ++iter) {

        for (size_t k = 0; k < 3*nw; ++k)
            for (size_t l = 0; l < 3*nw; ++l)
                Efd[k] += ddt[3*nw*k + l]*dip[l];

        double deltadip = 0.0;

        for (size_t i = 0; i < 3*nw; ++i) {
            dip[i] = polarM*(Efq[i] + Efd[i]);
            dip[i] = dmix*dip[i] + (1.0 - dmix)*dip_prev[i];

            double delta = dip[i] - dip_prev[i];
            deltadip += delta*delta;
        }

        deltadip = std::sqrt(deltadip)*stath;

        if (deltadip < dip_tolerance) {
#ifndef NDEBUG
            converged = true;
#endif // NDEBUG
            break; // converged!
        } else {
            std::copy(dip, dip + 3*nw, dip_prev);
            std::fill(Efd, Efd + 3*nw, 0.0);
        }
    }

    assert(converged);

    double Eelec(0), Eind(0);
    for (size_t i = 0; i < 3*nw; ++i) {
        Eelec += phi[i]*charge[i];
        Eind -= dip[i]*Efq[i];
    }

    Eelec /= 2;
    Eind /= 2;

    // charge-dipole and dipole-dipole

    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;

        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            size_t j3 = 3*j;

            for (size_t l = 1; l < 3; ++l) {
                size_t j3l = j3 + l;
                size_t jh = 3*j3l;

                double Rij[3], Rsq(0), diR(0);
                for (size_t k = 0; k < 3; ++k) {
                    Rij[k] = msite[i3 + k] - crd[jh + k];
                    Rsq += Rij[k]*Rij[k];
                    diR += dip[i3 + k]*Rij[k];
                }

                double ts1, ts2;
                ttm3f_bits::smear12(std::sqrt(Rsq), AA_HM, ts1, ts2);

                for (size_t k = 0; k < 3; ++k) {
                    double derij =
                        charge[j3l]*(3*ts2*diR*Rij[k] - ts1*dip[i3 + k]);

                    dM[i3 + k] += derij;
                    grad[jh + k] -= derij;
                }

                phi[j3l] -= ts1*diR;
            }

            // charge-dipole: M-sites
            double Rij[3], Rsq(0), diR(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = msite[i3 + k] - msite[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dip[i3 + k]*Rij[k];
            }

            double ts1, ts2, ts3;
            ttm3f_bits::smear123(std::sqrt(Rsq), AA_MM, ts1, ts2, ts3);

            for (size_t k = 0; k < 3; ++k) {
                double derij =
                    charge[j3]*(3*ts2*diR*Rij[k] - ts1*dip[i3 + k]);

                dM[i3 + k] += derij;
                dM[j3 + k] -= derij;
            }

            phi[j3] -= ts1*diR;

            // dipole-dipole : M-sites
            if (j > i) {

                double djR(0), didj(0);
                for (size_t k = 0; k < 3; ++k) {
                    djR += dip[j3 + k]*Rij[k];
                    didj += dip[i3 + k]*dip[j3 + k];
                }

                for (size_t k = 0; k < 3; ++k) {
                    double derij =
                       - 3*ts2*(didj*Rij[k] + djR*dip[i3 + k]
                                            + diR*dip[j3 + k])
                       + 15*ts3*diR*djR*Rij[k];

                    dM[i3 + k] += derij;
                    dM[j3 + k] -= derij;
                }
            }
        }
    }

    // distribute M-site gradients

    using ttm3f_bits::gamma1;
    using ttm3f_bits::gamma2;

    for (size_t i = 0; i < nw; ++i) {
        size_t i3 = 3*i;
        size_t i9 = 3*i3;

        for (size_t k = 0; k < 3; ++k) {
            double tmp = dM[i3 + k] - charge[i3]*Efq[i3 + k];

            grad[i9 + 0 + k] += gamma1*tmp; // O
            grad[i9 + 3 + k] += gamma2*tmp; // H
            grad[i9 + 6 + k] += gamma2*tmp; // H
        }
    }

    // derivatives from the adjustable charges of the NASA PES

    for (size_t n = 0; n < nw; ++n) {
        size_t n3 = 3*n;
        size_t io  = 9*n + 0;
        size_t ih1 = 9*n + 3;
        size_t ih2 = 9*n + 6;

        for (size_t k = 0; k < 3; ++k) {
            grad[ih1 + k] += GRDQ(0, 0, k)*phi[n3 + 1]  // phi(h1)
                           + GRDQ(0, 1, k)*phi[n3 + 2]  // phi(h2)
                           + GRDQ(0, 2, k)*phi[n3 + 0]; // phi(M)

            grad[ih2 + k] += GRDQ(1, 0, k)*phi[n3 + 1]  // phi(h1)
                           + GRDQ(1, 1, k)*phi[n3 + 2]  // phi(h2)
                           + GRDQ(1, 2, k)*phi[n3 + 0]; // phi(M)

            grad[io + k] += GRDQ(2, 0, k)*phi[n3 + 1]  // phi(h1)
                          + GRDQ(2, 1, k)*phi[n3 + 2]  // phi(h2)
                          + GRDQ(2, 2, k)*phi[n3 + 0]; // phi(M)
        }
    }

    return Eint + Evdw + Eelec + Eind;
}

//----------------------------------------------------------------------------//

} // namespace h2o

////////////////////////////////////////////////////////////////////////////////
