#ifndef TTM4_HBB2_X3B_H
#define TTM4_HBB2_X3B_H

#include <string>
#include <cstddef>

#include "x3b.h"
#include "ttm4-es.h"

#include "h2o.h"

namespace h2o {
    
class ttm4_hbb2_x3b : public Potential  {
public:
    ttm4_hbb2_x3b();

    void load(const char*);

    const char* name() const
    {
        return m_name.c_str();
    }

    double operator()(size_t nw, const double*); // O H H O H H
    double operator()(size_t nw, const double*, double*); // O H H O H H

private:
    std::string m_name;

private:
    ttm::ttm4_es m_ttm4;

private:
    x3b m_x3b;
};

} // namespace h2o

#endif // TTM4_HBB2_X3B_H
