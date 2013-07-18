#include <cmath>
#include <limits>

#include "macros.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

const double EPS = std::numeric_limits<double>::epsilon();
const double FPMIN = std::numeric_limits<double>::min()/EPS;

//----------------------------------------------------------------------------//

double gser23(const double& x)
{
    double a = 2.0/3.0;

    double ap = a;
    double sum = 1.0/a;
    double del = sum;

    for (;;) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (std::fabs(del) < std::fabs(sum)*EPS) {
            return sum*std::exp(- x + a*std::log(x));
        }
    }
}

//----------------------------------------------------------------------------//

double gcf23(const double& x)
{
    double a = 2.0/3.0;

    double b = x + 1.0 - a;
    double c = 1.0/FPMIN;
    double d = 1.0/b;
    double h = d;

    for (int i = 1;; ++i) {
        double an = -i*(i - a);

        b += 2.0;
        d = an*d + b;
        if (std::fabs(d) < FPMIN)
            d = FPMIN;

        c = b + an/c;

        if (std::fabs(c) < FPMIN)
            c = FPMIN;

        d = 1.0/d;
        double del = d*c;
        h *= del;

        if (std::fabs(del - 1.0) <= EPS)
            break;
    }

    return std::exp( - x + a*std::log(x))*h;
}

//----------------------------------------------------------------------------//

//
// gammq(2/3, x)*gammln(2/3)
//

double gammq23(const double& x)
{
    double exp_gammln_23 = 1.354117939426400e+00;
    return (x < 5.0/3.0 ? exp_gammln_23 - gser23(x) : gcf23(x));
}

//----------------------------------------------------------------------------//

const double a_XX = 0.175;
const double a_XX_13 = std::pow(a_XX, 1.0/3.0);

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace ttm3f_bits {

//----------------------------------------------------------------------------//

void smear0(const double& RESTRICT r12, const double& RESTRICT AA,
            double& RESTRICT ts0)
{
    double dri = 1.0/r12;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);

    ts0 = (1.0 - exp1 + a_XX_13*rA*gammq23(a_XX*rA3))*dri;
}

//----------------------------------------------------------------------------//

void smear1(const double& RESTRICT r12, const double& RESTRICT AA,
            double& RESTRICT ts1)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);

    ts1 = (1.0 - exp1)*dri*drsqi;
}

//----------------------------------------------------------------------------//

void smear01(const double& RESTRICT r12, const double& RESTRICT AA,
             double& RESTRICT ts0, double& RESTRICT ts1)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);

    ts0 = (1.0 - exp1 + a_XX_13*rA*gammq23(a_XX*rA3))*dri;
    ts1 = (1.0 - exp1)*dri*drsqi;
}

//----------------------------------------------------------------------------//

void smear12(const double& RESTRICT r12, const double& RESTRICT AA,
             double& RESTRICT ts1, double& RESTRICT ts2)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);

    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - exp1*a_XX/(AA*AA*AA))*drsqi;
}

//----------------------------------------------------------------------------//

void smear012(const double& RESTRICT r12, const double& RESTRICT AA,
              double& RESTRICT ts0, double& RESTRICT ts1, double& RESTRICT ts2)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);

    ts0 = (1.0 - exp1 + a_XX_13*rA*gammq23(a_XX*rA3))*dri;
    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - exp1*a_XX/(AA*AA*AA))*drsqi;
}

//----------------------------------------------------------------------------//

void smear123(const double& RESTRICT r12, const double& RESTRICT AA,
              double& RESTRICT ts1, double& RESTRICT ts2, double& RESTRICT ts3)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);
    double AA3 = AA*AA*AA;
    double AA6 = AA3*AA3;

    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - exp1*a_XX/AA3)*drsqi;
    ts3 = (ts2 - 0.6*exp1*r12*a_XX*a_XX/AA6)*drsqi;
}

//----------------------------------------------------------------------------//

void smear0123(const double& RESTRICT r12, const double& RESTRICT AA,
               double& RESTRICT ts0, double& RESTRICT ts1,
               double& RESTRICT ts2, double& RESTRICT ts3)
{
    double dri = 1.0/r12;
    double drsqi = dri*dri;

    double rA = r12/AA;
    double rA3 = rA*rA*rA;
    double exp1 = std::exp(-a_XX*rA3);
    double AA3 = AA*AA*AA;
    double AA6 = AA3*AA3;

    ts0 = (1.0 - exp1 + a_XX_13*rA*gammq23(a_XX*rA3))*dri;
    ts1 = (1.0 - exp1)*dri*drsqi;
    ts2 = (ts1 - exp1*a_XX/AA3)*drsqi;
    ts3 = (ts2 - 0.6*exp1*r12*a_XX*a_XX/AA6)*drsqi;
}

//----------------------------------------------------------------------------//

} // namespace ttm3f_bits

////////////////////////////////////////////////////////////////////////////////
