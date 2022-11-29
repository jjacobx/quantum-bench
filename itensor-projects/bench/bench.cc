//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "../helpers/ops.h"
#include "../helpers/io.h"

#include <chrono>

using namespace itensor;
using namespace std;

#define QREG_DEFAULT 4
#define INIT_DEFAULT "|0..0>"
#define CUTOFF_DEFAULT 1E-16
#define MAXDIM_DEFAULT 1073741824
#define DEPTH_DEFAULT 16

MPS applyRandomMPS(MPS mps, int depth, int maxdim, double cutoff);
vector<MPO> constructRandomMPOs(MPS mps, int depth);
Cplx wideOverlap(MPS left, vector<MPO> mpos, MPS right, int maxdim, double cutoff);

int main(int argc, char *argv[]) {
    int qreg_size = QREG_DEFAULT;
    string init_state = INIT_DEFAULT;
    int maxdim = MAXDIM_DEFAULT;
    double cutoff = CUTOFF_DEFAULT;
    int depth = DEPTH_DEFAULT;

    set_args(argc, argv, qreg_size, init_state, maxdim, cutoff, depth);
    int verbose = set_verbose();

    srand(2140);

    auto init_mps = initMPS(qreg_size, init_state);
    auto measure_mps = initMPS(qreg_size, "|0..0>");

    auto tstart = chrono::steady_clock::now();

    auto result_mps = applyRandomMPS(init_mps, depth, maxdim, cutoff);

    auto tstop = chrono::steady_clock::now();
    double tdiff = chrono::duration<double, milli>(tstop - tstart).count();
    cout << "Full simulation time: " << tdiff << " ms" << endl;

    // PrintData(result_mps);
    printfln("Norm: %f", norm(result_mps));
    printfln("Max link dim: %f", maxLinkDim(result_mps));
    printfln("Avg link dim: %f", averageLinkDim(result_mps));

    Cplx amp = innerC(init_mps, result_mps);
    cout << "Amplitude: " << amp << endl << endl;

    srand(2140);

    tstart = chrono::steady_clock::now();

    vector<MPO> random_circuit = constructRandomMPOs(init_mps, depth);
    amp = wideOverlap(init_mps, random_circuit, init_mps, maxdim, cutoff);

    tstop = chrono::steady_clock::now();
    tdiff = chrono::duration<double, milli>(tstop - tstart).count();

    cout << "Overlap time: " << tdiff << " ms" << endl;
    cout << "Amplitude: " << amp << endl;

    return 0;
}

MPS applyRandomMPS(MPS mps, int depth, int maxdim, double cutoff) {
    SiteSet sites = SpinHalf(siteInds(mps));
    for (int i = 0; i < depth; i++) {
        MPO rmpo = MPO(sites); // random MPO layer
        for (int j = 1; j <= length(mps); j++)
            rmpo = nmultMPO(prime(rmpo), popRAND(sites, j));
        mps = applyMPO(rmpo, noPrime(mps), {"MaxDim=", maxdim, "Cutoff=", cutoff});
        
        MPO empo = MPO(sites); // entangle MPO layer
        for (int j = 1 + i % 2; j < length(mps); j += 2)
            empo = nmultMPO(prime(empo), popCROT(sites, j, j + 1, 1));
        mps = applyMPO(empo, noPrime(mps), {"MaxDim=", maxdim, "Cutoff=", cutoff});
    }

    return mps;
}

vector<MPO> constructRandomMPOs(MPS mps, int depth) {
    SiteSet sites = SpinHalf(siteInds(mps));
    vector<MPO> mpos;

    for (int i = 0; i < depth; i++) {
        MPO rmpo = MPO(sites); // random MPO layer
        for (int j = 1; j <= length(mps); j++){
            rmpo = nmultMPO(prime(rmpo), popRAND(sites, j));
            rmpo.mapPrime(2, 1);
        }
        
        MPO empo = MPO(sites); // entangle MPO layer
        for (int j = 1 + i % 2; j < length(mps); j += 2) {
            empo = nmultMPO(prime(empo), popCROT(sites, j, j + 1, 1));
            empo.mapPrime(2, 1);
        }

        MPO layer = nmultMPO(prime(empo), rmpo);
        layer.mapPrime(2, 1);

        mpos.push_back(prime(layer, i));
    }

    return(mpos);
}

Cplx wideOverlap(MPS left, vector<MPO> mpos, MPS right, int maxdim, double cutoff) {
    SiteSet sites = SpinHalf(siteInds(right));
    int depth = mpos.size();
    left.prime(depth);

    if (depth == 0)
        return inner(left, right);
    if (depth == 1)
        return inner(left, mpos.at(0), right);

    int i, j; // remember values after loops exit
    for (i = 0; i < (depth - 1) / 2; i++)
        right = applyMPO(mpos.at(i), right, {"MaxDim=", maxdim, "Cutoff=", cutoff});
    for (j = depth - 1; j > (depth + 1) / 2; j--) 
        left = applyMPO(mpos.at(j), left, {"MaxDim=", maxdim, "Cutoff=", cutoff});
    
    return innerC(left, mpos.at(j), mpos.at(i), right);
}

// mps = applyMPO(popRAND(sites, j), mps, {"Cutoff=", cutoff});
// mps = applyMPO(popCROT(sites, j, j + 1, 1), mps, {"Cutoff=", cutoff});
