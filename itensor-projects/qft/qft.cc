//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "../helpers/ops.h"
#include "../helpers/io.h"

using namespace itensor;
using namespace std;

#define QREG_DEFAULT 4
#define INIT_DEFAULT "|0..0>"
#define CUTOFF_DEFAULT 1E-4
#define MAXDIM_DEFAULT 1073741824

ITensor applyQFT_tensor(ITensor init);
MPS applyQFT_mps(MPS mps, double cutoff);

int main(int argc, char *argv[]) {
    int qreg_size = QREG_DEFAULT;
    string init_state = INIT_DEFAULT;
    double cutoff = CUTOFF_DEFAULT;
    int maxdim = MAXDIM_DEFAULT;

    set_args(argc, argv, qreg_size, init_state, maxdim, cutoff);
    int verbose = set_verbose();

    // cout << "Number of qubits: " << qreg_size << endl;
    // cout << "Init state: " << init_state << endl;
    // cout << "Cutoff: " << cutoff << endl;
    // cout << "Verbose: " << verbose << endl;

    srand(2140);


    // ============ Using standard tensors ============ //

    // auto init_tensor = initTensor(qreg_size, init_state);
    // auto result_tensor = applyQFT_tensor(init_tensor);
    // PrintData(result_tensor);


    // ============ Using MPS ============ //


    auto init_mps = initMPS(qreg_size, init_state);

    auto tstart = chrono::steady_clock::now();

    auto result_mps = applyQFT_mps(init_mps, cutoff);

    auto tstop = chrono::steady_clock::now();
    double tdiff = chrono::duration<double, milli>(tstop - tstart).count();
    cout << "Full simulation time: " << tdiff << " ms" << endl;

    // PrintData(result_mps);
    printfln("Norm: %f", norm(result_mps));
    printfln("Max link dim: %f", maxLinkDim(result_mps));
    printfln("Avg link dim: %f", averageLinkDim(result_mps));

    // auto measure_mps = MPS(InitState(spin_sites, "Up"));
    // PrintData(innerC(result_mps, measure_mps));
    // PrintData(ContractMPS(result_mps));

    return 0;
}


ITensor applyQFT_tensor(ITensor init) {
    IndexSet indices = inds(init);
    for (int i = 1; i <= order(init); i++)
        init *= makeQFT_STEP(indices, i);
    return init;
}

MPS applyQFT_mps(MPS mps, double cutoff) {
    SiteSet sites = SpinHalf(siteInds(mps));
    for (int i = 1; i <= length(mps); i++) {
        mps = applyMPO(popH(sites, i), mps, {"Cutoff=", cutoff});
        for (int j = i + 1; j <= length(mps); j++)
            mps = applyMPO(popCROT(sites, j, i, j - i), mps, {"Cutoff=", cutoff});
    }

    return mps;
}
