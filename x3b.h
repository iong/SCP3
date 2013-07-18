#ifndef X3B_H
#define X3B_H

#include <cstddef>

namespace h2o {

namespace x3b_bits {
    static const size_t ncoeffs = 131;
    extern double kOO;
    extern double kOH;
    extern double kHH;
    extern double r0_OO;
    extern double r0_OH;
    extern double r0_HH;
    extern double r3i;
    extern double r3f;
    extern double coeffs[];
    extern bool initialized;
    extern void load(const char*);
}


struct x3b {

    x3b();

    double operator()(const double* w1,
                      const double* w2,
                      const double* w3,
                      double* g1,
                      double* g2,
                      double* g3) const;
private:
    double f_switch3(const double&, double&) const;

};
} // namespace h2o

#endif // X3B_H
