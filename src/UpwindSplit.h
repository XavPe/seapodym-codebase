#ifndef SEAPODYM_UPWIND_SPLIT_H
#define SEAPODYM_UPWIND_SPLIT_H

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace seapodym_upwind {

inline double parse_epsilon(const char* text, const char* setting)
{
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(text, &end);
    if (end == text || *end || errno == ERANGE || !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(setting) + " must be finite and nonnegative");
    }
    return value;
}

inline bool in_band(double velocity, double width)
{
    return std::abs(velocity) < width;
}

enum class Mode { QuadraticC1, QuarticC2 };

inline Mode parse_mode(const std::string& value)
{
    if (value == "quadratic_C1") return Mode::QuadraticC1;
    if (value == "quartic_C2") return Mode::QuarticC2;
    throw std::invalid_argument("upwind_smoothing_mode/@value must be quadratic_C1 or quartic_C2");
}

inline const char* mode_name(Mode mode)
{
    return mode == Mode::QuarticC2 ? "quartic_C2" : "quadratic_C1";
}

struct Split {
    double positive;
    double negative;
    double d_positive;
    double d_negative;
};

// Precondition: |velocity| < width. The caller retains the native stencil elsewhere.
inline Split inside_split(double velocity, double width, Mode mode)
{
    const double relative = velocity / width;
    if (mode == Mode::QuarticC2) {
        const double p = 1.0 + relative, m = 1.0 - relative;
        return {width*p*p*p*(3.0-relative)/16.0,
                -width*m*m*m*(3.0+relative)/16.0,
                p*p*(2.0-relative)/4.0, m*m*(2.0+relative)/4.0};
    }
    const double d_positive = (1.0 + relative) * 0.5;
    const double d_negative = (1.0 - relative) * 0.5;
    return {
        (width * d_positive) * d_positive,
        -(width * d_negative) * d_negative,
        d_positive,
        d_negative,
    };
}

} // namespace seapodym_upwind

#endif
