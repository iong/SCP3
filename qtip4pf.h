#ifndef QTIP4PF_H
#define QTIP4PF_H

#include <cstddef>

#include "h2o.h"

namespace h2o {

class qtip4pf : public Potential  {
public:
    qtip4pf();
    ~qtip4pf();

    const char* name() const
    {
        return "q-TIP4P/F";
    }

    double operator()(size_t nw, const double*); // O H H O H H
    double operator()(size_t nw, const double*, double*); // O H H O H H

    // molecular dipoles (size = 3*nw); master only
    const double* dip(size_t, const double*);

private:
    size_t  m_nw;
    double* m_mem;

    void allocate(size_t);

private:
    bool m_dip_ready;
};

} // namespace h2o

#endif // QTIP4PF_H
