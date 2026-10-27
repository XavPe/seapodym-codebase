#include "Cmdopt.h"
#include <cstring>

int seapodym_fd_point(int argc, char** argv);

int seapodym_coupled(const char* parfile, Cmdopt* cmdopt);
int seapodym_phases(const char* parfile);

int main(int argc, char** argv) {


    // Dispatch before Cmdopt: an invalid FD request must never start optimization.
    if (argc > 1 && std::strncmp(argv[1], "--fd-", 5) == 0) {
        if (std::strcmp(argv[1], "--fd-point") != 0) {
            std::cerr << "Unknown FD option; use --fd-point --help\n";
            return 2;
        }
        return seapodym_fd_point(argc - 1, argv + 1);
    }

	Cmdopt cmdopt(argc, argv);

	if (cmdopt.cmp_regime!=-11)
		return seapodym_coupled(argv[argc-1], &cmdopt);
	else //not supported for now, will not work
		return seapodym_phases(argv[argc-1]);
}


