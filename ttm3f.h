#ifndef TTM3F_H
#define TTM3F_H

#include <cstddef>

#include "PES.h"

namespace h2o {

class ttm3f : public PES {
public:
    ttm3f();
    ~ttm3f();

    const char* name() const
    {
        return "TTM3-F";
    }

    double operator()(size_t nw, const double*); // O H H O H H
    double operator()(size_t nw, const double*, double*); // O H H O H H

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);
};

} // namespace h2o

#endif // TTM3F_H
