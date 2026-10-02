#ifndef SEAPODYM_AAT_EXPORT_H
#define SEAPODYM_AAT_EXPORT_H

#include <cmath>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace aat {
struct Parameter {
    std::string path, species;
    double lower, upper;
    std::string id() const { return path + "/@" + species; }
};
inline void field(const std::string& s) {
    if (s.empty() || s.find_first_of("\t\r\n") != std::string::npos ||
        s.find('\0') != std::string::npos)
        throw std::runtime_error("invalid AAT metadata field");
}
inline void validate(const std::vector<Parameter>& pars) {
    if (pars.empty()) throw std::runtime_error("no AAT parameters");
    std::set<std::string> ids;
    for (size_t i = 0; i < pars.size(); ++i) {
        const Parameter& p = pars[i];
        field(p.path); field(p.species);
        if (p.path[0] != '/' || !ids.insert(p.id()).second ||
            !std::isfinite(p.lower) || !std::isfinite(p.upper) ||
            !(p.lower < p.upper) || !std::isfinite(p.upper-p.lower))
            throw std::runtime_error("invalid or duplicate AAT parameter: " + p.id());
    }
}
// Exclusive creation also rejects existing symlinks. On failure remove only
// files created by this invocation, never files from a previous experiment.
class Files {
    std::vector<std::string> names;
    bool complete;
public:
    Files() : complete(false) { names.reserve(3); }
    ~Files() {
        if (!complete)
            for (size_t i = 0; i < names.size(); ++i) ::unlink(names[i].c_str());
    }
    void put(const std::string& name, const std::string& bytes) {
        names.push_back(name);
        int fd = ::open(name.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
        if (fd < 0) {
            names.pop_back();
            throw std::runtime_error("cannot create " + name + ": " + std::strerror(errno));
        }
        size_t offset = 0;
        while (offset < bytes.size()) {
            ssize_t n = ::write(fd, bytes.data()+offset, bytes.size()-offset);
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) {
                int saved = errno;
                ::close(fd);
                throw std::runtime_error("cannot write " + name + ": " + std::strerror(saved));
            }
            offset += static_cast<size_t>(n);
        }
        if (::fsync(fd) != 0) {
            int saved = errno;
            ::close(fd);
            throw std::runtime_error("cannot flush " + name + ": " + std::strerror(saved));
        }
        if (::close(fd) != 0)
            throw std::runtime_error("cannot close " + name + ": " + std::strerror(errno));
    }
    void commit() { complete = true; }
};
inline void write(const char* source, int seed, const std::vector<Parameter>& pars,
                  const std::vector<std::vector<double> >& rows) {
    validate(pars);
    field(source);
    if (rows.empty()) throw std::runtime_error("no AAT experiments");
    std::ostringstream table, meta;
    table.imbue(std::locale::classic()); meta.imbue(std::locale::classic());
    table << std::setprecision(std::numeric_limits<double>::max_digits10);
    meta << std::setprecision(std::numeric_limits<double>::max_digits10);
    table << "exp_id\tseed";
    meta << "format_version\t1\nseed\t" << seed << "\nn_aat\t" << rows.size()
         << "\nn_parameters\t" << pars.size() << "\nsource_parfile\t" << source
         << "\nsource_copy\taat_source_parfile.xml\n"
         << "index\tparameter_id\txml_path\tspecies_attribute\tlower\tupper\n";
    for (size_t i = 0; i < pars.size(); ++i) {
        table << '\t' << pars[i].id();
        meta << i+1 << '\t' << pars[i].id() << '\t' << pars[i].path << '\t'
             << pars[i].species << '\t' << pars[i].lower << '\t' << pars[i].upper << '\n';
    }
    table << '\n';
    for (size_t k = 0; k < rows.size(); ++k) {
        if (rows[k].size() != pars.size()) throw std::runtime_error("invalid AAT row size");
        table << k+1 << '\t' << seed;
        for (size_t i = 0; i < pars.size(); ++i) {
            double v = rows[k][i];
            if (!std::isfinite(v) || v < pars[i].lower || v > pars[i].upper)
                throw std::runtime_error("invalid AAT physical value: " + pars[i].id());
            table << '\t' << v;
        }
        table << '\n';
    }
    std::ifstream input(source, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open source parfile");
    std::string copy;
    char buffer[8192];
    while (input.read(buffer, sizeof buffer) || input.gcount())
        copy.append(buffer, static_cast<size_t>(input.gcount()));
    if (input.bad() || !input.eof()) throw std::runtime_error("cannot read source parfile");
    input.clear(); input.close();
    if (input.fail()) throw std::runtime_error("cannot close source parfile");
    if (!table || !meta) throw std::runtime_error("cannot format AAT export");
    Files files;
    files.put("aat_source_parfile.xml", copy);
    files.put("aat_parameters.tsv", table.str());
    // Metadata is written last; a complete export has all three closed files.
    files.put("aat_metadata.tsv", meta.str());
    files.commit();
}
} // namespace aat
#endif
