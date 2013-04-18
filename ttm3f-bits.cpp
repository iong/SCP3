#include <cmath>
#include <cassert>

#include <limits>
#include <algorithm>

#include "macros.h"
#include "ttm3f-bits.h"
#include "ttm3f-bits-smear.h"

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

double U_CC_HH
    (const double* RESTRICT H1, const double& RESTRICT q1,
     const double* RESTRICT H2, const double& RESTRICT q2)
{
    double r2(0);
    for (size_t k = 0; k < 3; ++k) {
        const double dx = H1[k] - H2[k];
        r2 += dx*dx;
    }

    double ts0;
    ttm3f_bits::smear0(std::sqrt(r2), AA_HH, ts0);

    return q1*q2*ts0;
}

//----------------------------------------------------------------------------//

double U_CC_HH
    (const double* RESTRICT H1, const double& RESTRICT q1,
     const double* RESTRICT H2, const double& RESTRICT q2,
     double* RESTRICT Efq1, double& RESTRICT phi1,
     double* RESTRICT Efq2, double& RESTRICT phi2)
{
    double r2(0), dr[3];
    for (size_t k = 0; k < 3; ++k) {
        dr[k] = H1[k] - H2[k];
        r2 += dr[k]*dr[k];
    }

    double ts0, ts1;
    ttm3f_bits::smear01(std::sqrt(r2), AA_HH, ts0, ts1);

    for (size_t k = 0; k < 3; ++k) {
        Efq1[k] += q2*ts1*dr[k];
        Efq2[k] -= q1*ts1*dr[k];
    }

    phi1 += q2*ts0;
    phi2 += q1*ts0;

    return q1*q2*ts0;
}

//----------------------------------------------------------------------------//

double U_CC_HM
    (const double* RESTRICT H, const double& RESTRICT qH,
     const double* RESTRICT M, const double& RESTRICT qM,
     double* RESTRICT Efq)
{
    double rMH[3], r2(0);
    for (size_t k = 0; k < 3; ++k) {
        rMH[k] = M[k] - H[k];
        r2 += rMH[k]*rMH[k];
    }

    double ts0, ts1;
    ttm3f_bits::smear01(std::sqrt(r2), AA_HM, ts0, ts1);

    for (size_t k = 0; k < 3; ++k)
        Efq[k] += qH*ts1*rMH[k];

    return qH*qM*ts0;
}

//----------------------------------------------------------------------------//

double U_CC_HM
    (const double* RESTRICT H, const double& RESTRICT qH,
     const double* RESTRICT M, const double& RESTRICT qM,
     double* RESTRICT EfqH, double& RESTRICT phiH,
     double* RESTRICT EfqM, double& RESTRICT phiM)
{
    double rMH[3], r2(0);
    for (size_t k = 0; k < 3; ++k) {
        rMH[k] = M[k] - H[k];
        r2 += rMH[k]*rMH[k];
    }

    double ts0, ts1;
    ttm3f_bits::smear01(std::sqrt(r2), AA_HM, ts0, ts1);

    phiH += qM*ts0;
    phiM += qH*ts0;

    for (size_t k = 0; k < 3; ++k) {
        EfqM[k] += qH*ts1*rMH[k];
        EfqH[k] -= qM*ts1*rMH[k];
    }

    return qH*qM*ts0;
}

//----------------------------------------------------------------------------//

//
// electric field created by hydrogen at M-site
//

void Efq_HM(const double* RESTRICT H, const double& RESTRICT qH,
            const double* RESTRICT M, double* RESTRICT Efq)
{
    double dr[3], r2(0);
    for (size_t k = 0; k < 3; ++k) {
        dr[k] = M[k] - H[k];
        r2 += dr[k]*dr[k];
    }

    double ts1;
    ttm3f_bits::smear1(std::sqrt(r2), AA_HM, ts1);

    for (size_t k = 0; k < 3; ++k)
        Efq[k] += qH*ts1*dr[k];
}

//----------------------------------------------------------------------------//

void grad_CD // charge-dipole gradients helper
    (const double* RESTRICT rC, const double& RESTRICT Q,
     const double* RESTRICT rD, const double* RESTRICT dip,
     const double& RESTRICT AA,
     double& RESTRICT phiC,
     double* RESTRICT dC, double* RESTRICT dD)
{
    double Rij[3], Rsq(0), diR(0);
    for (size_t k = 0; k < 3; ++k) {
        Rij[k] = rD[k] - rC[k];
        Rsq += Rij[k]*Rij[k];
        diR += dip[k]*Rij[k];
    }

    double ts1, ts2;
    ttm3f_bits::smear12(std::sqrt(Rsq), AA, ts1, ts2);

    for (size_t k = 0; k < 3; ++k) {
        const double derij = Q*(3*ts2*diR*Rij[k] - ts1*dip[k]);

        dD[k] += derij;
        dC[k] -= derij;
    }

    phiC -= ts1*diR;
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ttm3f_bits {

//----------------------------------------------------------------------------//

void U_coul_CC_CD_DD_2w
   (const double* RESTRICT Ma, // (a)    M-site : x y z
    const double* RESTRICT Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qa, // (a)   charges : M q1 q2
    const double* RESTRICT Mb, // (b)    M-site : x y z
    const double* RESTRICT Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qb, // (b)   charges : M q1 q2
    double& Eelec, double& Eind)
{
    // H-H energies
    Eelec = U_CC_HH(Ha + 0, Qa[1], Hb + 0, Qb[1])
          + U_CC_HH(Ha + 0, Qa[1], Hb + 3, Qb[2])
          + U_CC_HH(Ha + 3, Qa[2], Hb + 0, Qb[1])
          + U_CC_HH(Ha + 3, Qa[2], Hb + 3, Qb[2]);

    // electric field at the M-sites
    double Efq_a[3] = {0.0, 0.0, 0.0};
    double Efq_b[3] = {0.0, 0.0, 0.0};

    // H-M energies (and fields at M-sites)
    Eelec += U_CC_HM(Ha + 0, Qa[1], Mb, Qb[0], Efq_b)
           + U_CC_HM(Ha + 3, Qa[2], Mb, Qb[0], Efq_b)
           + U_CC_HM(Hb + 0, Qb[1], Ma, Qa[0], Efq_a)
           + U_CC_HM(Hb + 3, Qb[2], Ma, Qa[0], Efq_a);

    // M-M energy & fields
    double x[3], r2(0);
    for (size_t k = 0; k < 3; ++k) {
        x[k] = Ma[k] - Mb[k];
        r2 += x[k]*x[k];
    }

    double ts0, ts1, ts2;
    smear012(std::sqrt(r2), AA_MM, ts0, ts1, ts2);

    Eelec += Qa[0]*Qb[0]*ts0;

    double xEa(0), xEb(0);
    for (size_t k = 0; k < 3; ++k) {
        Efq_a[k] += Qb[0]*ts1*x[k];
        Efq_b[k] -= Qa[0]*ts1*x[k];
        xEa += Efq_a[k]*x[k];
        xEb += Efq_b[k]*x[k];
    }

    // induced dipoles

    double TEa[3], TEb[3], xTa(0), xTb(0);
    for (size_t k = 0; k < 3; ++k) {
        TEa[k] = polarM*(3*ts2*x[k]*xEa - ts1*Efq_a[k]);
        TEb[k] = polarM*(3*ts2*x[k]*xEb - ts1*Efq_b[k]);
        xTa += TEa[k]*x[k];
        xTb += TEb[k]*x[k];
    }

    const double tmp0 = polarM*ts1;
    const double tmp1 = 1.0 - tmp0*tmp0;
    const double tmp2 = polarM/tmp1;
    const double tmp3 = 3*polarM*polarM*ts2*(2*ts1 - 3*ts2*r2);
    const double tmp4 = tmp3/(tmp1 + r2*tmp3);

    Eind = 0.0;
    for (size_t k = 0; k < 3; ++k) {
        const double dip_a = tmp2*(Efq_a[k] + TEb[k] - x[k]*tmp4*(xEa + xTb));
        const double dip_b = tmp2*(Efq_b[k] + TEa[k] - x[k]*tmp4*(xEb + xTa));

        Eind -= dip_a*Efq_a[k] + dip_b*Efq_b[k];
    }

    Eind /= 2;
}

//----------------------------------------------------------------------------//

void U_coul_CD_DD
   (size_t nw,
    const double* RESTRICT OHH, // coordinates : 9*nw : O H H O H H
    const double* RESTRICT M,  //     M-sites : 3*nw : x y z x y z
    const double* RESTRICT Q, //     charges : 3*nw : M H H M H H
    double* RESTRICT mem,    // size >= 12*nw + 9*nw*nw = 3*nw*(4 + 3*nw)
    double& RESTRICT Eind)
{
    // setup pointers
    double* Efq = mem;
    double* Efd = Efq + 3*nw;
    double* ddt = Efd + 3*nw;
    double* dip = ddt + (3*nw)*(3*nw);
    double* dip_prev = dip + 3*nw;

    // zero out Efq/Efd
    std::fill(mem, mem + 6*nw, 0.0);

    // compute electric field at the M-sites
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;

        // hydrogens
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            const size_t j3 = 3*j;

            Efq_HM(OHH + 3*(j3 + 1), Q[j3 + 1], M + i3, Efq + i3);
            Efq_HM(OHH + 3*(j3 + 2), Q[j3 + 2], M + i3, Efq + i3);
        }

        // M-sites
        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;

            double dr[3], r2(0);
            for (size_t k = 0; k < 3; ++k) {
                dr[k] = M[i3 + k] - M[j3 + k];
                r2 += dr[k]*dr[k];
            }

            double ts1, ts2;
            smear12(std::sqrt(r2), AA_MM, ts1, ts2);

            for (size_t k = 0; k < 3; ++k) {
                Efq[i3 + k] += Q[j3]*ts1*dr[k];
                Efq[j3 + k] -= Q[i3]*ts1*dr[k];
            }

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*dr[0]*dr[0] - ts1;
            dd3[1][1] = 3.0*ts2*dr[1]*dr[1] - ts1;
            dd3[2][2] = 3.0*ts2*dr[2]*dr[2] - ts1;
            dd3[0][1] = 3.0*ts2*dr[0]*dr[1];
            dd3[0][2] = 3.0*ts2*dr[0]*dr[2];
            dd3[1][2] = 3.0*ts2*dr[1]*dr[2];
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

    const double dmix = 0.7;
    const double stath = DEBYE/CHARGECON/(std::sqrt(double(nw)));

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

            const double delta = dip[i] - dip_prev[i];
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

    Eind = 0.0;
    for (size_t i = 0; i < 3*nw; ++i)
        Eind -= dip[i]*Efq[i];

    Eind /= 2;
}

//----------------------------------------------------------------------------//

void U_coul_CD_DD_3w
   (const double* RESTRICT Ma, // (a)    M-site : x y z
    const double* RESTRICT Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qa, // (a)   charges : M q1 q2
    const double* RESTRICT Mb, // (b)    M-site : x y z
    const double* RESTRICT Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qb, // (b)   charges : M q1 q2
    const double* RESTRICT Mc, // (c)    M-site : x y z
    const double* RESTRICT Hc, // (c) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qc, // (c)   charges : M q1 q2
    double& RESTRICT Eind)
{
    const size_t nw = 3;

    double OHH[9*nw];
    double   M[3*nw];
    double   Q[3*nw];

    for (size_t k = 0; k < 3; ++k) {

        Q[3*0 + k] = Qa[k];
        Q[3*1 + k] = Qb[k];
        Q[3*2 + k] = Qc[k];

        for (size_t l = 0; l < 2; ++l) {
            OHH[9*0 + 3*(l + 1) + k] = Ha[3*l + k];
            OHH[9*1 + 3*(l + 1) + k] = Hb[3*l + k];
            OHH[9*2 + 3*(l + 1) + k] = Hc[3*l + k];
        }

        M[3*0 + k] = Ma[k];
        M[3*1 + k] = Mb[k];
        M[3*2 + k] = Mc[k];
    }

    double mem[3*nw*(4 + 3*nw)];

    U_coul_CD_DD(nw, OHH, M, Q, mem, Eind);
}

//----------------------------------------------------------------------------//

void U_coul_CC_CD_DD_2w
   (const double* RESTRICT Ma, // (a)    M-site : x y z
    const double* RESTRICT Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qa, // (a)   charges : M q1 q2
    const double* RESTRICT Mb, // (b)    M-site : x y z
    const double* RESTRICT Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* RESTRICT Qb, // (b)   charges : M q1 q2
    double& RESTRICT Eelec, double& RESTRICT Eind,
    const double* RESTRICT grdQa, // (a) charge gradient
    const double* RESTRICT grdQb, // (b) charge gradient
    double* RESTRICT dEelec_a, double* RESTRICT dEind_a,
    double* RESTRICT dEelec_b, double* RESTRICT dEind_b)
{
    // potential due to charges: Ma Ha1 Ha2 Mb Hb1 Hb2
    double phiC[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // electric field due to charges (same layout)
    double Efq[18];
    std::fill(Efq, Efq + 18, 0.0);

    // H-H energies/fields
    Eelec =
  U_CC_HH(Ha + 0, Qa[1], Hb + 0, Qb[1], Efq + 3*1, phiC[1], Efq + 3*4, phiC[4])
+ U_CC_HH(Ha + 0, Qa[1], Hb + 3, Qb[2], Efq + 3*1, phiC[1], Efq + 3*5, phiC[5])
+ U_CC_HH(Ha + 3, Qa[2], Hb + 0, Qb[1], Efq + 3*2, phiC[2], Efq + 3*4, phiC[4])
+ U_CC_HH(Ha + 3, Qa[2], Hb + 3, Qb[2], Efq + 3*2, phiC[2], Efq + 3*5, phiC[5]);

    // H-M energies (and fields at M-sites)
    Eelec +=
      U_CC_HM(Ha + 0, Qa[1], Mb, Qb[0], Efq + 3*1, phiC[1], Efq + 3*3, phiC[3])
    + U_CC_HM(Ha + 3, Qa[2], Mb, Qb[0], Efq + 3*2, phiC[2], Efq + 3*3, phiC[3])
    + U_CC_HM(Hb + 0, Qb[1], Ma, Qa[0], Efq + 3*4, phiC[4], Efq + 3*0, phiC[0])
    + U_CC_HM(Hb + 3, Qb[2], Ma, Qa[0], Efq + 3*5, phiC[5], Efq + 3*0, phiC[0]);

    // M-M energy & fields
    double x[3], r2(0);
    for (size_t k = 0; k < 3; ++k) {
        x[k] = Ma[k] - Mb[k];
        r2 += x[k]*x[k];
    }

    double ts0, ts1, ts2, ts3;
    smear0123(std::sqrt(r2), AA_MM, ts0, ts1, ts2, ts3);

    Eelec += Qa[0]*Qb[0]*ts0;

    double xEa(0), xEb(0);
    for (size_t k = 0; k < 3; ++k) {
        Efq[3*0 + k] += Qb[0]*ts1*x[k];
        Efq[3*3 + k] -= Qa[0]*ts1*x[k];
        xEa += Efq[3*0 + k]*x[k];
        xEb += Efq[3*3 + k]*x[k];
    }

    phiC[0] += Qb[0]*ts0;
    phiC[3] += Qa[0]*ts0;

    // induced dipoles

    double TEa[3], TEb[3], xTa(0), xTb(0);
    for (size_t k = 0; k < 3; ++k) {
        TEa[k] = polarM*(3*ts2*x[k]*xEa - ts1*Efq[0 + k]);
        TEb[k] = polarM*(3*ts2*x[k]*xEb - ts1*Efq[9 + k]);
        xTa += TEa[k]*x[k];
        xTb += TEb[k]*x[k];
    }

    const double tmp0 = polarM*ts1;
    const double tmp1 = 1.0 - tmp0*tmp0;
    const double tmp2 = polarM/tmp1;
    const double tmp3 = 3*polarM*polarM*ts2*(2*ts1 - 3*ts2*r2);
    const double tmp4 = tmp3/(tmp1 + r2*tmp3);

    Eind = 0.0;

    double dip_a[3], dip_b[3];
    for (size_t k = 0; k < 3; ++k) {
        dip_a[k] = tmp2*(Efq[0 + k] + TEb[k] - x[k]*tmp4*(xEa + xTb));
        dip_b[k] = tmp2*(Efq[9 + k] + TEa[k] - x[k]*tmp4*(xEb + xTa));

        Eind -= dip_a[k]*Efq[0 + k] + dip_b[k]*Efq[9 + k];
    }

    Eind /= 2;

    // Eelec gradients

    using ttm3f_bits::gamma1;
    using ttm3f_bits::gamma2;

    for (size_t k = 0; k < 3; ++k) {
        const double dMa = Qa[0]*Efq[3*0 + k];
        dEelec_a[0 + k] = - gamma1*dMa; // Oa
        dEelec_a[3 + k] = - (Qa[1]*Efq[3*1 + k] + gamma2*dMa); // Ha1
        dEelec_a[6 + k] = - (Qa[2]*Efq[3*2 + k] + gamma2*dMa); // Ha2

        const double dMb = Qb[0]*Efq[3*3 + k];
        dEelec_b[0 + k] = - gamma1*dMb; // Ob
        dEelec_b[3 + k] = - (Qb[1]*Efq[3*4 + k] + gamma2*dMb); // Hb1
        dEelec_b[6 + k] = - (Qb[2]*Efq[3*5 + k] + gamma2*dMb); // Hb2
    }

#define GRDQa(i,j,k) grdQa[k + 3*(j + 3*i)]
#define GRDQb(i,j,k) grdQb[k + 3*(j + 3*i)]

    for (size_t k = 0; k < 3; ++k) {
        dEelec_a[3 + k] += GRDQa(0, 0, k)*phiC[1]  // phi(Ha1)
                         + GRDQa(0, 1, k)*phiC[2]  // phi(Ha2)
                         + GRDQa(0, 2, k)*phiC[0]; // phi(Ma)

        dEelec_a[6 + k] += GRDQa(1, 0, k)*phiC[1]  // phi(Ha1)
                         + GRDQa(1, 1, k)*phiC[2]  // phi(Ha2)
                         + GRDQa(1, 2, k)*phiC[0]; // phi(Ma)

        dEelec_a[0 + k] += GRDQa(2, 0, k)*phiC[1]  // phi(Ha1)
                         + GRDQa(2, 1, k)*phiC[2]  // phi(Ha2)
                         + GRDQa(2, 2, k)*phiC[0]; // phi(Ma)

        dEelec_b[3 + k] += GRDQb(0, 0, k)*phiC[4]  // phi(Hb1)
                         + GRDQb(0, 1, k)*phiC[5]  // phi(Hb2)
                         + GRDQb(0, 2, k)*phiC[3]; // phi(Mb)

        dEelec_b[6 + k] += GRDQb(1, 0, k)*phiC[4]  // phi(Hb1)
                         + GRDQb(1, 1, k)*phiC[5]  // phi(Hb2)
                         + GRDQb(1, 2, k)*phiC[3]; // phi(Mb)

        dEelec_b[0 + k] += GRDQb(2, 0, k)*phiC[4]  // phi(Hb1)
                         + GRDQb(2, 1, k)*phiC[5]  // phi(Hb2)
                         + GRDQb(2, 2, k)*phiC[3]; // phi(Mb)
    }

    // Eind gradients

    for (size_t k = 0; k < 9; ++k)
        dEind_a[k] = dEind_b[k] = 0.0;

    double phiD[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // M H H M H H

    // charge-dipole: hydrogens

    grad_CD(Ha + 0, Qa[1], Mb, dip_b, AA_HM, phiD[1], dEind_a + 3, dEind_b + 0);
    grad_CD(Ha + 3, Qa[2], Mb, dip_b, AA_HM, phiD[2], dEind_a + 6, dEind_b + 0);
    grad_CD(Hb + 0, Qb[1], Ma, dip_a, AA_HM, phiD[4], dEind_b + 3, dEind_a + 0);
    grad_CD(Hb + 3, Qb[2], Ma, dip_a, AA_HM, phiD[5], dEind_b + 6, dEind_a + 0);

    // charge-dipole, dipole-dipole: M-sites

    {
        // x[3] = Ma - Mb, r2 = x^2

        double daR(0), dbR(0), dadb(0);
        for (size_t k = 0; k < 3; ++k) {
            daR += dip_a[k]*x[k];
            dbR += dip_b[k]*x[k];
            dadb += dip_a[k]*dip_b[k];
        }

        for (size_t k = 0; k < 3; ++k) {
            const double der_a = Qb[0]*(3*ts2*daR*x[k] - ts1*dip_a[k]);

            dEind_a[k] += der_a;
            dEind_b[k] -= der_a;

            const double der_b = Qa[0]*(3*ts2*dbR*x[k] - ts1*dip_b[k]);

            dEind_b[k] += der_b;
            dEind_a[k] -= der_b;
        }

        phiD[0] += ts1*dbR;
        phiD[3] -= ts1*daR;

        // dipole-dipole
        for (size_t k = 0; k < 3; ++k) {
            const double der =
                   - 3*ts2*(dadb*x[k] + dbR*dip_a[k] + daR*dip_b[k])
                   + 15*ts3*daR*dbR*x[k];

            dEind_a[k] += der;
            dEind_b[k] -= der;
        }
    }

    // re-position M-site gradient

    for (size_t k = 0; k < 3; ++k) {
        const double dMa = dEind_a[0 + k];
        dEind_a[0 + k] = gamma1*dMa; // Oa
        dEind_a[3 + k] += gamma2*dMa; // Ha1
        dEind_a[6 + k] += gamma2*dMa; // Ha2

        const double dMb = dEind_b[0 + k];
        dEind_b[0 + k] = gamma1*dMb; // Ob
        dEind_b[3 + k] += gamma2*dMb; // Hb1
        dEind_b[6 + k] += gamma2*dMb; // Hb2
    }

    // NASA DMS effects

    for (size_t k = 0; k < 3; ++k) {
        dEind_a[3 + k] += GRDQa(0, 0, k)*phiD[1]  // phi(Ha1)
                        + GRDQa(0, 1, k)*phiD[2]  // phi(Ha2)
                        + GRDQa(0, 2, k)*phiD[0]; // phi(Ma)

        dEind_a[6 + k] += GRDQa(1, 0, k)*phiD[1]  // phi(Ha1)
                        + GRDQa(1, 1, k)*phiD[2]  // phi(Ha2)
                        + GRDQa(1, 2, k)*phiD[0]; // phi(Ma)

        dEind_a[0 + k] += GRDQa(2, 0, k)*phiD[1]  // phi(Ha1)
                        + GRDQa(2, 1, k)*phiD[2]  // phi(Ha2)
                        + GRDQa(2, 2, k)*phiD[0]; // phi(Ma)

        dEind_b[3 + k] += GRDQb(0, 0, k)*phiD[4]  // phi(Hb1)
                        + GRDQb(0, 1, k)*phiD[5]  // phi(Hb2)
                        + GRDQb(0, 2, k)*phiD[3]; // phi(Mb)

        dEind_b[6 + k] += GRDQb(1, 0, k)*phiD[4]  // phi(Hb1)
                        + GRDQb(1, 1, k)*phiD[5]  // phi(Hb2)
                        + GRDQb(1, 2, k)*phiD[3]; // phi(Mb)

        dEind_b[0 + k] += GRDQb(2, 0, k)*phiD[4]  // phi(Hb1)
                        + GRDQb(2, 1, k)*phiD[5]  // phi(Hb2)
                        + GRDQb(2, 2, k)*phiD[3]; // phi(Mb)
    }
}

//----------------------------------------------------------------------------//

void U_coul_CD_DD
   (size_t nw,
    const double* RESTRICT OHH, // coordinates : 9*nw : O H H O H H
    const double* RESTRICT M,  //     M-sites : 3*nw : x y z x y z
    const double* RESTRICT Q, //     charges : 3*nw : M H H M H H
    double* RESTRICT mem,    // size >= 18*nw + 9*nw*nw = 9*nw*(2 + nw)
    double& RESTRICT Eind,
    const double* RESTRICT grdq, // charge gradients
    double* RESTRICT dOHH) // gradients : 9*nw O H H O H H
{
    // setup pointers
    double* Efq = mem;
    double* Efd = Efq + 3*nw;
    double* phi = Efd + 3*nw;
    double* dM  = phi + 3*nw;
    double* ddt = dM + 3*nw;
    double* dip = ddt + (3*nw)*(3*nw);
    double* dip_prev = dip + 3*nw;

    // zero out Efq/Efd/phi/dM
    std::fill(mem, mem + 12*nw, 0.0);

    // compute electric field at the M-sites
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;

        // hydrogens
        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            const size_t j3 = 3*j;

            Efq_HM(OHH + 3*(j3 + 1), Q[j3 + 1], M + i3, Efq + i3);
            Efq_HM(OHH + 3*(j3 + 2), Q[j3 + 2], M + i3, Efq + i3);
        }

        // M-sites
        for (size_t j = i + 1; j < nw; ++j) {
            const size_t j3 = 3*j;

            double dr[3], r2(0);
            for (size_t k = 0; k < 3; ++k) {
                dr[k] = M[i3 + k] - M[j3 + k];
                r2 += dr[k]*dr[k];
            }

            double ts1, ts2;
            smear12(std::sqrt(r2), AA_MM, ts1, ts2);

            for (size_t k = 0; k < 3; ++k) {
                Efq[i3 + k] += Q[j3]*ts1*dr[k];
                Efq[j3 + k] -= Q[i3]*ts1*dr[k];
            }

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*dr[0]*dr[0] - ts1;
            dd3[1][1] = 3.0*ts2*dr[1]*dr[1] - ts1;
            dd3[2][2] = 3.0*ts2*dr[2]*dr[2] - ts1;
            dd3[0][1] = 3.0*ts2*dr[0]*dr[1];
            dd3[0][2] = 3.0*ts2*dr[0]*dr[2];
            dd3[1][2] = 3.0*ts2*dr[1]*dr[2];
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

    const double dmix = 0.7;
    const double stath = DEBYE/CHARGECON/(std::sqrt(double(nw)));

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

            const double delta = dip[i] - dip_prev[i];
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

    Eind = 0.0;
    for (size_t i = 0; i < 3*nw; ++i)
        Eind -= dip[i]*Efq[i];

    Eind /= 2;

    // === gradients ===

    // external loop over M-sites
    for (size_t i = 0; i < nw; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = 0; j < nw; ++j) {
            if (i == j)
                continue;

            const size_t j3 = 3*j;

            // charge-dipole: hydrogens
            for (size_t l = 1; l < 3; ++l) {
                const size_t j3l = j3 + l;
                const size_t jh = 3*j3l;

                grad_CD(OHH + jh, Q[j3l],
                        M + i3, dip + i3,
                        AA_HM, phi[j3l],
                        dOHH + jh, dM + i3);
            }

            // charge-dipole: M-sites
            double Rij[3], Rsq(0), diR(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = M[i3 + k] - M[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dip[i3 + k]*Rij[k];
            }

            double ts1, ts2, ts3;
            smear123(std::sqrt(Rsq), AA_MM, ts1, ts2, ts3);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                    Q[j3]*(3*ts2*diR*Rij[k] - ts1*dip[i3 + k]);

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
                    const double derij =
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
        const size_t i3 = 3*i;
        const size_t i9 = 3*i3;

        for (size_t k = 0; k < 3; ++k) {
            dOHH[i9 + 0 + k] += gamma1*dM[i3 + k]; // O
            dOHH[i9 + 3 + k] += gamma2*dM[i3 + k]; // H
            dOHH[i9 + 6 + k] += gamma2*dM[i3 + k]; // H
        }
    }

    // derivatives from the adjustable charges of the NASA PES

#define GRDQ(i,j,k) grdq[k + 3*(j + 3*(i + 3*n))]

    for (size_t n = 0; n < nw; ++n) {
        const size_t n3 = 3*n;
        const size_t io  = 9*n + 0;
        const size_t ih1 = 9*n + 3;
        const size_t ih2 = 9*n + 6;

        for (size_t k = 0; k < 3; ++k) {
            dOHH[ih1 + k] += GRDQ(0, 0, k)*phi[n3 + 1]  // phi(h1)
                           + GRDQ(0, 1, k)*phi[n3 + 2]  // phi(h2)
                           + GRDQ(0, 2, k)*phi[n3 + 0]; // phi(M)

            dOHH[ih2 + k] += GRDQ(1, 0, k)*phi[n3 + 1]  // phi(h1)
                           + GRDQ(1, 1, k)*phi[n3 + 2]  // phi(h2)
                           + GRDQ(1, 2, k)*phi[n3 + 0]; // phi(M)

            dOHH[io + k] += GRDQ(2, 0, k)*phi[n3 + 1]  // phi(h1)
                          + GRDQ(2, 1, k)*phi[n3 + 2]  // phi(h2)
                          + GRDQ(2, 2, k)*phi[n3 + 0]; // phi(M)
        }
    }
}

//----------------------------------------------------------------------------//

} // namespace ttm3f_bits

////////////////////////////////////////////////////////////////////////////////

