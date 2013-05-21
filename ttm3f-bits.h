#ifndef TTM3F_BITS_H
#define TTM3F_BITS_H

//
// pieces of TTM3-F needed for Bowman's potential
//

namespace ttm3f_bits {

// vdw
const double vdwC = -0.72298855E+03;
const double vdwD =  0.10211829E+06;
const double vdwE =  0.37170376E+01;

// dms
const double dms_param1 = 0.5;
const double dms_param2 = 0.9578;
const double dms_param3 = 0.012;

// M-site positioning
const double gammaM = 0.46;
const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

// constants
const double CHARGECON = 18.22234397655801030455;

// handy bits
void U_coul_CC_CD_DD_2w
   (const double* Ma, // (a)    M-site : x y z
    const double* Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qa, // (a)   charges : M q1 q2
    const double* Mb, // (b)    M-site : x y z
    const double* Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qb, // (b)   charges : M q1 q2
    double& Eelec, double& Eind);

void U_coul_CD_DD
   (size_t nw,
    const double* OHH, // coordinates : 9*nw : O H H O H H
    const double* M,  //     M-sites : 3*nw : x y z x y z
    const double* Q, //     charges : 3*nw : M H H M H H
    double* mem,    // size >= 12*nw + 9*nw*nw = 3*nw*(4 + 3*nw)
    double& Eind);

void U_coul_CD_DD_3w
   (const double* Ma, // (a)    M-site : x y z
    const double* Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qa, // (a)   charges : M q1 q2
    const double* Mb, // (b)    M-site : x y z
    const double* Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qb, // (b)   charges : M q1 q2
    const double* Mc, // (c)    M-site : x y z
    const double* Hc, // (c) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qc, // (c)   charges : M q1 q2
    double& Eind);

//----------------------------- with gradients -------------------------------//

void U_coul_CC_CD_DD_2w
   (const double* Ma, // (a)    M-site : x y z
    const double* Ha, // (a) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qa, // (a)   charges : M q1 q2
    const double* Mb, // (b)    M-site : x y z
    const double* Hb, // (b) Hydrogens : x1 y1 z1 x2 y2 z2
    const double* Qb, // (b)   charges : M q1 q2
    double& Eelec, double& Eind,
    const double* grdQa, // (a) charge gradient
    const double* grdQb, // (b) charge gradient
    double* dEelec_a, double* dEind_a, // (a) gradients
    double* dEelec_b, double* dEind_b); // (b) gradients

void U_coul_CD_DD
   (size_t nw,
    const double* OHH, // coordinates : 9*nw : O H H O H H
    const double* M,  //     M-sites : 3*nw : x y z x y z
    const double* Q, //     charges : 3*nw : M H H M H H
    double* mem,    // size >= 18*nw + 9*nw*nw = 9*nw*(2 + nw)
    double& Eind,
    const double* grdq, // charge gradients
    double* dOHH);

} // namespace ttm3f_bits

#endif // TTM3F_BITS_H
