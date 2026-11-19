#ifndef SEAPODYM_PARAMETER_DIAGNOSTICS_H
#define SEAPODYM_PARAMETER_DIAGNOSTICS_H

#include <admodel.h>
#include <cmath>
#include <limits>

namespace seapodym_diagnostics {

inline double parameter_jacobian(double internal, double lower, double upper)
{
    double penalty = 0.0;
    return ndfboundp(internal, lower, upper, penalty);
}

inline double physical_cv(double estimate, double internal_se, double jacobian)
{
    if (!std::isfinite(estimate) || estimate == 0.0 ||
        !std::isfinite(internal_se) || internal_se < 0.0 || !std::isfinite(jacobian))
        return std::numeric_limits<double>::quiet_NaN();
    return std::fabs(jacobian) * internal_se / std::fabs(estimate);
}

} // namespace seapodym_diagnostics
#endif
