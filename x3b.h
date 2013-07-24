#ifndef X3B_H
#define X3B_H

#include <cstddef>

namespace h2o {

struct x3b {

    x3b();

    void load(const char*);

    static const size_t ncoeffs = 131;

    inline const double& r3i() const;
    inline const double& r3f() const;

    double operator()(const double* w1,
                      const double* w2,
                      const double* w3,
                      double* g1,
                      double* g2,
                      double* g3) const;

private:
    double m_kOO;
    double m_kOH;
    double m_kHH;

private:
    double m_r0_OO;
    double m_r0_OH;
    double m_r0_HH;

private:
    double m_r3i;
    double m_r3f;

    double f_switch3(const double&, double&) const;

private:
    double m_coeffs[ncoeffs];
};

inline const double& x3b::r3i() const
{
    return m_r3i;
}

inline const double& x3b::r3f() const
{
    return m_r3f;
}

} // namespace h2o

#endif // X3B_H
