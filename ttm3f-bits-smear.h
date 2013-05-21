#ifndef TTM3F_BITS_SMEAR_H
#define TTM3F_BITS_SMEAR_H

//
// smeared interactions for TTM3-F
//

namespace ttm3f_bits {

void smear0(const double& r12, const double& AA, double& ts0);

void smear1(const double& r12, const double& AA, double& ts1);

void smear01(const double& r12, const double& AA,
             double& ts0, double& ts1);

void smear12(const double& r12, const double& AA,
             double& ts1, double& ts2);

void smear012(const double& r12, const double& AA,
              double& ts0, double& ts1, double& ts2);

void smear123(const double& r12, const double& AA,
              double& ts1, double& ts2, double& ts3);

void smear0123(const double& r12, const double& AA,
               double& ts0, double& ts1, double& ts2, double& ts3);

} // namespace ttm3f_bits

#endif // TTM3F_BITS_SMEAR_H
