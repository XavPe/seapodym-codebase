#ifndef SEAPODYM_AAT_RUN_STATUS_H
#define SEAPODYM_AAT_RUN_STATUS_H

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>

namespace aat {

// Indices are the native SEAPODYM indices (no conversion to age or coordinates).
// A rejected evaluation has no completed objective: never publish its partial sum.
struct ZeroHabitatEvent {
    bool detected;
    int species, age, i, j, step, year, month, day;
    ZeroHabitatEvent() : detected(false), species(0), age(0), i(0), j(0),
        step(0), year(0), month(0), day(0) {}
};

inline bool finite_result(const double* values, int count) {
    for (int i = 0; i < count; ++i)
        if (!std::isfinite(values[i])) return false;
    return true;
}

inline void write_status(std::ostream& out, int id, const char* status,
                         const double* values, const ZeroHabitatEvent& event) {
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << id << '\t' << status;
    // Same component ordering as the existing native AAT table.
    for (int i = 0; i < 7; ++i) {
        out << '\t';
        if (values && std::isfinite(values[i])) out << values[i];
        else out << "NA";
    }
    if (event.detected)
        out << '\t' << event.species << '\t' << event.age << '\t' << event.i
            << '\t' << event.j << '\t' << event.step << '\t' << event.year
            << '\t' << event.month << '\t' << event.day;
    else
        out << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    // Flush every completed/rejected draw. An I/O failure must stop the campaign.
    out << std::endl;
}

} // namespace aat
#endif
