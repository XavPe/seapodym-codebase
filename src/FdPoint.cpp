// One deterministic forward evaluation at an exact internal coordinate.
// Optional seapodym_clte mode: no optimizer or adjoint; equations unchanged.
#include <fvar.hpp>
#include "SeapodymCoupled.h"
#include "AatExport.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

static std::string quoted(const std::string& s) {
    std::ostringstream o; o << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') o << '\\' << c;
        else if (c < 32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else o << c;
    }
    o << '"'; return o.str();
}
static void put(const std::string& path, const std::string& data) {
    aat::Files files; files.put(path, data + "\n"); files.commit();
}
static void help() {
    std::cout << "seapodym_clte --fd-point 0.2.0\n"
      << "--protocol | --version | --help\n"
      << "--describe --output metadata.json PARFILE\n"
      << "--evaluate --x-vector vector.txt --output result.json PARFILE\n"
      << "Optional: -mv/--mv BYTES -mg/--mg ENTRIES -mc/--mc BYTES (buffers, not RAM limits).\n"
      << "Keeps cached-forcing optimization initialization, includes bound penalty.\n"
      << "Forward only; no optimization and no large output fields.\n";
}
int seapodym_fd_point(int argc, char** argv) {
    try {
        if (argc == 2) {
            const std::string a=argv[1];
            if (a=="--protocol") {
                std::cout << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"tool_version\":\"0.2.0\"}\n";
                return 0;
            }
            if (a=="--help" || a=="-h") { help(); return 0; }
            if (a=="--version") { std::cout << "seapodym_clte --fd-point 0.2.0\n"; return 0; }
        }
        bool describe=false, evaluate=false;
        std::string output, vectorfile, parfile;
        long mv=200000000L, mg=100000L, mc=1000000L;
        for (int i=1;i<argc;++i) {
            const std::string a=argv[i];
            if (a=="--describe") describe=true;
            else if (a=="--evaluate") evaluate=true;
            else if (a=="--output" || a=="--x-vector" || a=="--mv" || a=="--mg" || a=="--mc" || a=="-mv" || a=="-mg" || a=="-mc") {
                if (++i>=argc) throw std::runtime_error("missing option value");
                if (a=="--output") output=argv[i];
                else if (a=="--x-vector") vectorfile=argv[i];
                else {
                    std::string v=argv[i]; size_t used=0;
                    long n=std::stol(v,&used);
                    if (used!=v.size() || n<=0) throw std::runtime_error("invalid buffer size");
                    if (a=="--mv" || a=="-mv") mv=n; else if (a=="--mg" || a=="-mg") mg=n; else mc=n;
                }
            } else if (a.empty() || a[0]=='-' || !parfile.empty()) throw std::runtime_error("unknown option or extra argument");
            else parfile=a;
        }
        if (describe==evaluate || output.empty() || parfile.empty() || (evaluate && vectorfile.empty()))
            throw std::runtime_error("choose exactly --describe or --evaluate, with required paths");
        // Refuse an existing output before allocating model state.
        std::ifstream exists(output.c_str());
        if (exists.good()) throw std::runtime_error("output already exists");
        std::ifstream input(parfile.c_str());
        if (!input.good()) throw std::runtime_error("parfile not readable");
        if (evaluate) {
            std::ifstream v(vectorfile.c_str());
            if (!v.good()) throw std::runtime_error("vector not readable");
        }
        gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
        gradient_structure::set_NO_DERIVATIVES();
        gradient_structure::set_GRADSTACK_BUFFER_SIZE(mg);
        gradient_structure::set_CMPDIF_BUFFER_SIZE(mc);
        gradient_structure gs(mv);
        SeapodymCoupled sc(parfile.c_str());
        if (sc.param->flag_coupling) throw std::runtime_error("forage coupling unsupported in FD v0.2");
        if (sc.param->nb_EEZ || sc.param->nb_tag_files || sc.param->tags_only)
            throw std::runtime_error("tagging/EEZ not validated for the FD v0.2 input profile");
        // Keep gcalc=true: same cached-forcing/density path as optimization.
        sc.param->set_gradcalc(true);
        sc.param->set_scalc(false); // Unlike native SA: keep bound penalty.
        const int n=sc.nvarcalc();
        if (n<=0) throw std::runtime_error("no free parameters");
        independent_variables x(1,n); adstring_array names(1,n);
        sc.xinit(x,names);
        std::ostringstream o; o.imbue(std::locale::classic());
        o << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (describe) {
            o << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"n_parameters\":" << n << ",\"parameters\":[";
            for (int i=1;i<=n;++i) {
                if (i>1) o << ',';
                o << "{\"index\":" << i << ",\"name\":" << quoted((char*)names[i])
                  << ",\"xml_path\":" << quoted(sc.param->fd_path(i))
                  << ",\"species\":" << quoted(sc.param->fd_species(i))
                  << ",\"lower\":" << sc.param->fd_lower(i)
                  << ",\"upper\":" << sc.param->fd_upper(i)
                  << ",\"physical\":" << sc.param->get_parval(i)
                  << ",\"internal\":" << x[i] << '}';
            }
            o << "],\"initialization\":\"native_xinit_cached_forcings\"}";
            put(output,o.str()); return 0;
        }
        std::ifstream v(vectorfile.c_str()); v.imbue(std::locale::classic());
        int supplied=0; v >> supplied;
        if (supplied!=n) throw std::runtime_error("vector dimension mismatch");
        for (int i=1;i<=n;++i) {
            double q; if (!(v>>q) || !std::isfinite(q) || q < -1.0 || q > 1.0)
                throw std::runtime_error("invalid internal coordinate");
            x[i]=q;
        }
        std::string extra; if (v>>extra) throw std::runtime_error("extra vector values");
        sc.prerun_model();
        const double total=sc.run_coupled((dvar_vector)x,false);
        const double c=sc.get_clike(), lf=sc.get_lflike(), stock=sc.get_stocklike();
        const double tags=sc.get_taglike(), larvae=sc.get_larvaelike(), spawn=sc.get_spawninglike();
        // Same deterministic bound transform as run_coupled; after forward only.
        const double penalty=value(sc.reset((dvar_vector)x));
        const double values[]={total,c,lf,stock,tags,larvae,spawn,penalty};
        for (double q:values) if (!std::isfinite(q)) throw std::runtime_error("nonfinite result");
        o << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"status\":\"ok\",\"n_parameters\":" << n
          << ",\"x\":[";
        for (int i=1;i<=n;++i) { if(i>1)o<<','; o<<x[i]; }
        o << "],\"physical\":[";
        for (int i=1;i<=n;++i) { if(i>1)o<<','; o<<sc.param->get_parval(i); }
        o << "],\"components\":{\"total\":" << total << ",\"captures\":" << c
          << ",\"LF\":" << lf << ",\"stock\":" << stock << ",\"marquage\":" << tags
          << ",\"larves\":" << larvae << ",\"gonades\":" << spawn << ",\"penalite\":" << penalty
          << "},\"initialization\":\"native_xinit_cached_forcings\",\"optimization\":false}";
        put(output,o.str());
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FD_POINT_ERROR: " << e.what() << '\n'; return 2;
    }
}
