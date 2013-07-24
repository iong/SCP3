#ifndef PS_H
#define PS_H

namespace ps {

double pot_nasa(const double*, double*);
void dms_nasa(const double&, const double&, const double&,
              const double*, double* q3, double* dq3, bool ttm3);

} // namespace ps

#endif // PS_H
