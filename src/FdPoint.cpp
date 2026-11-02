// Exact-point FD/adjoint backend v0.3.0. One process, one forward, at most one reverse.
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
#include <cmath>

static std::string quoted(const std::string& s) {
    std::ostringstream o; o << '"';
    for (unsigned char c:s) {
        if (c=='"' || c=='\\') o << '\\' << c;
        else if (c<32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else o << c;
    }
    o << '"'; return o.str();
}
static void put(const std::string& path,const std::string& data) {
    aat::Files files; files.put(path,data+"\n"); files.commit();
}
static void help() {
    std::cout << "seapodym_clte --fd-point 0.3.0\n"
      << "--protocol | --version | --help\n"
      << "--describe --output metadata.json PARFILE\n"
      << "--evaluate --x-vector x.txt --output result.json PARFILE\n"
      << "--adjoint --component total|captures|LF|stock|marquage|larves|gonades|penalite --x-vector x.txt --output result.json PARFILE\n"
      << "--exclude-penalty : scalc=true, matching native Taylor; default includes penalty.\n"
      << "Buffers: --mv BYTES --mg ENTRIES --mc BYTES. Adjoint requires explicit --mg and --mc.\n"
      << "No optimization. Cached forcings (gcalc=true), one species/no tags/EEZ/coupling for adjoint.\n";
}
int seapodym_fd_point(int argc,char** argv) {
    try {
        if (argc==2) {
            const std::string a=argv[1];
            if (a=="--protocol") {
                std::cout << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"tool_version\":\"0.3.0\",\"capabilities\":[\"forward\",\"adjoint_components\",\"penalty_modes\"]}\n";
                return 0;
            }
            if (a=="--help" || a=="-h") { help();return 0; }
            if (a=="--version") { std::cout << "seapodym_clte --fd-point 0.3.0\n";return 0; }
        }
        bool describe=false,evaluate=false,adjoint=false,exclude_penalty=false,explicit_mg=false,explicit_mc=false;
        std::string output,vectorfile,parfile,component;
        long mv=200000000L,mg=100000L,mc=1000000L;
        for (int i=1;i<argc;++i) {
            const std::string a=argv[i];
            if (a=="--describe") describe=true;
            else if (a=="--evaluate") evaluate=true;
            else if (a=="--adjoint") adjoint=true;
            else if (a=="--exclude-penalty") exclude_penalty=true;
            else if (a=="--output" || a=="--x-vector" || a=="--component" || a=="--mv" || a=="--mg" || a=="--mc" || a=="-mv" || a=="-mg" || a=="-mc") {
                if (++i>=argc) throw std::runtime_error("missing option value");
                if (a=="--output") output=argv[i];
                else if (a=="--x-vector") vectorfile=argv[i];
                else if (a=="--component") component=argv[i];
                else {
                    const std::string value=argv[i];size_t used=0;
                    const long n=std::stol(value,&used);
                    if (used!=value.size() || n<=0) throw std::runtime_error("invalid buffer size");
                    if (a=="--mv" || a=="-mv") mv=n;
                    else if (a=="--mg" || a=="-mg") { mg=n;explicit_mg=true; }
                    else { mc=n;explicit_mc=true; }
                }
            } else if (a.empty() || a[0]=='-' || !parfile.empty()) throw std::runtime_error("unknown option or extra argument");
            else parfile=a;
        }
        if (int(describe)+int(evaluate)+int(adjoint)!=1 || output.empty() || parfile.empty() || (!describe && vectorfile.empty()))
            throw std::runtime_error("choose exactly --describe, --evaluate or --adjoint with required paths");
        const char* names[]={"total","captures","LF","stock","marquage","larves","gonades","penalite"};
        int target=-1;
        if (adjoint) {
            for (int i=0;i<8;++i) if (component==names[i]) target=i;
            if (target<0) throw std::runtime_error("unknown or missing adjoint component");
            if (!explicit_mg || !explicit_mc) throw std::runtime_error("adjoint requires explicit --mg and --mc buffers");
        } else if (!component.empty()) throw std::runtime_error("--component is reserved for --adjoint");
        if (std::ifstream(output.c_str()).good()) throw std::runtime_error("output already exists");
        if (!std::ifstream(parfile.c_str()).good()) throw std::runtime_error("parfile not readable");
        if (!describe && !std::ifstream(vectorfile.c_str()).good()) throw std::runtime_error("vector not readable");
        gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
        if (adjoint) gradient_structure::set_YES_DERIVATIVES();
        else gradient_structure::set_NO_DERIVATIVES();
        gradient_structure::set_GRADSTACK_BUFFER_SIZE(mg);
        gradient_structure::set_CMPDIF_BUFFER_SIZE(mc);
        gradient_structure gs(mv);
        SeapodymCoupled sc(parfile.c_str());
        if (sc.param->flag_coupling || sc.param->nb_EEZ || sc.param->nb_tag_files || sc.param->tags_only)
            throw std::runtime_error("coupling/tagging/EEZ not covered by this input profile");
        if (adjoint && (sc.param->get_nbspecies()!=1 || sum(sc.param->tag_like)!=0))
            throw std::runtime_error("adjoint components require one species and inactive tags");
        sc.param->set_gradcalc(true);
        sc.param->set_scalc(exclude_penalty);
        const int n=sc.nvarcalc();
        if (n<=0) throw std::runtime_error("no free parameters");
        independent_variables x(1,n);adstring_array parameter_names(1,n);
        sc.xinit(x,parameter_names);
        std::ostringstream o;o.imbue(std::locale::classic());
        o << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (describe) {
            o << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"n_parameters\":" << n << ",\"parameters\":[";
            for (int i=1;i<=n;++i) {
                if (i>1) o << ',';
                o << "{\"index\":" << i << ",\"name\":" << quoted((char*)parameter_names[i])
                  << ",\"xml_path\":" << quoted(sc.param->fd_path(i)) << ",\"species\":" << quoted(sc.param->fd_species(i))
                  << ",\"lower\":" << sc.param->fd_lower(i) << ",\"upper\":" << sc.param->fd_upper(i)
                  << ",\"physical\":" << sc.param->get_parval(i) << ",\"internal\":" << x[i] << '}';
            }
            o << "],\"initialization\":\"native_xinit_cached_forcings\"}";put(output,o.str());return 0;
        }
        std::ifstream vector(vectorfile.c_str());vector.imbue(std::locale::classic());int supplied=0;vector>>supplied;
        if (supplied!=n) throw std::runtime_error("vector dimension mismatch");
        for (int i=1;i<=n;++i) {
            double value;if (!(vector>>value) || !std::isfinite(value) || value < -1 || value > 1)
                throw std::runtime_error("invalid internal coordinate");
            x[i]=value;
        }
        std::string extra;if (vector>>extra) throw std::runtime_error("extra vector values");
        if (adjoint) { sc.fd_adjoint_target=target;sc.fd_selected.allocate(1,1); }
        sc.prerun_model();
        const double selected=sc.run_coupled((dvar_vector)x,false);
        dvector gradient(1,n);gradient.initialize();
        // No reset/getter performing AD work between selecting the scalar and gradcalc.
        if (adjoint) gradcalc(n,gradient);
        const double total=adjoint ? sc.fd_full_likelihood : selected;
        const double values[]={total,sc.get_clike(),sc.get_lflike(),sc.get_stocklike(),sc.get_taglike(),sc.get_larvaelike(),sc.get_spawninglike(),sc.fd_penalty_value};
        for (double v:values) if (!std::isfinite(v)) throw std::runtime_error("nonfinite cost");
        if (adjoint) for (int i=1;i<=n;++i) if (!std::isfinite(gradient[i])) throw std::runtime_error("nonfinite gradient");
        o << "{\"protocol\":\"seapodym-fd-point\",\"version\":1,\"status\":\"ok\",\"n_parameters\":" << n << ",\"x\":[";
        for (int i=1;i<=n;++i) { if (i>1) o << ',';o << x[i]; }
        o << "],\"physical\":[";
        for (int i=1;i<=n;++i) { if (i>1) o << ',';o << sc.param->get_parval(i); }
        o << "],\"components\":{";
        for (int i=0;i<8;++i) { if (i) o << ',';o << quoted(names[i]) << ':' << values[i]; }
        o << "},\"initialization\":\"native_xinit_cached_forcings\",\"optimization\":false,\"mode\":" << quoted(adjoint ? "adjoint" : "forward")
          << ",\"gcalc\":true,\"scalc\":" << (exclude_penalty ? "true" : "false")
          << ",\"penalty_mode\":" << quoted(exclude_penalty ? "exclude" : "include");
        if (adjoint) {
            o << ",\"adjoint_component\":" << quoted(component) << ",\"selected_cost\":" << selected << ",\"gradient\":[";
            for (int i=1;i<=n;++i) { if (i>1) o << ',';o << gradient[i]; }
            o << "],\"parameter_names\":[";
            for (int i=1;i<=n;++i) { if (i>1) o << ',';o << quoted((char*)parameter_names[i]); }
            o << "],\"gradient_coordinates\":\"internal\"";
        }
        o << '}';put(output,o.str());return 0;
    } catch (const std::exception& e) {
        std::cerr << "FD_POINT_ERROR: " << e.what() << '\n';return 2;
    }
}
