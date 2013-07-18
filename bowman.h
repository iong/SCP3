#ifndef BOWMAN_H
#define BOWMAN_H

#include <cstddef>

#include "PES.h"

namespace h2o {

class bowman : public PES  {
public:
    bowman();
    ~bowman();

    const char* name() const
    {
        return "WHBB";
    }

    double operator()(size_t nw, const double*); // O H H O H H
    double operator()(size_t nw, const double*, double*); // O H H O H H

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace h2o

#endif // BOWMAN_H
