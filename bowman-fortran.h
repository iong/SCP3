#ifndef BOWMAN_FORTRAN_H
#define BOWMAN_FORTRAN_H

namespace h2o { namespace fortran {

void pes2b_init();
bool pes2b_ready();

double pes2b(const double*);

void pes3b_init();
bool pes3b_ready();

double pes3b(const double*);

}} // namespace h2o::fortran

#endif // BOWMAN_FORTRAN_H
