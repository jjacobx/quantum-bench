#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

void set_args(int argc, char *argv[], int &qreg_size, string &init_state, int &maxdim, double &cutoff, int &depth) {
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--nq") {
            qreg_size = atoi(argv[++i]);
        } else if (arg == "--init") {
            init_state = argv[++i];
        } else if (arg == "--maxd") {
            maxdim = atoi(argv[++i]);
        } else if (arg == "--cut") {
            cutoff = atof(argv[++i]);
        } else if (arg == "--dep") {
            depth = atoi(argv[++i]);
        } else {
            string message = "Error: Unknown argument '" + arg + 
                "'! Use: ./bin --nq $NQUBITS --cut $CUTOFF";
            throw invalid_argument(message);
        }
    }
}

void set_args(int argc, char *argv[], int &qreg_size, string &init_state, int &maxdim, double &cutoff) {
    int depth;
    set_args(argc, argv, qreg_size, init_state, maxdim, cutoff, depth);
}

int set_verbose() {
    const char *tmp = getenv("VERBOSE");
    string verbose_str(tmp ? tmp : "");
    if (verbose_str == "1")
        return 1;
    else
        return 0;
}
