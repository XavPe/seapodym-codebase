#ifndef SEAPODYM_REGULARIZATION_H
#define SEAPODYM_REGULARIZATION_H
#include <cmath>
namespace seapodym_regularization {
// Keep the original arithmetic when disabled, including sqrt(g*g).
inline double habitat_abs(double g, double d, double eta) {
    if (eta == 0.0) return std::sqrt(g*g);
    const double q = d*g;
    return (q * (q/(std::hypot(q,eta)+eta)))/d;
}
inline double habitat_abs_derivative(double g, double d, double eta) {
    if (eta == 0.0) return g > 0 ? 1.0 : (g < 0 ? -1.0 : 0.0);
    const double q = d*g;
    return q/std::hypot(q,eta);
}
}
#endif
