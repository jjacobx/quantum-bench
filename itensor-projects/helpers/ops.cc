#include <stdexcept>
#include <vector>
#include "ops.h"

using namespace std;


// ========================================================================= //
// ----------------------------- ITensor gates ----------------------------- //
// ========================================================================= //

ITensor makeH(Index const& s) {
    auto H = ITensor(s, prime(s));
    H.set(s = 1, prime(s) = 1, +1 / sqrt(2));
    H.set(s = 1, prime(s) = 2, +1 / sqrt(2));
    H.set(s = 2, prime(s) = 1, +1 / sqrt(2));
    H.set(s = 2, prime(s) = 2, -1 / sqrt(2));
    return(H);
}

ITensor makeCROT(Index const& s, Index const& t, int k) {
    auto CROT = ITensor(s, t, prime(s), prime(t));
    CROT.set(s = 1, t = 1, prime(s) = 1, prime(t) = 1, 1.0);
    CROT.set(s = 1, t = 2, prime(s) = 1, prime(t) = 2, 1.0);
    CROT.set(s = 2, t = 1, prime(s) = 2, prime(t) = 1, 1.0);
    CROT.set(s = 2, t = 2, prime(s) = 2, prime(t) = 2, exp(1_i * Pi / (1 << k)));
    return CROT;
}

ITensor makeQFT_STEP(IndexSet indices, int which) {
    auto QFT_STEP = makeH(prime(indices[which - 1], which - 1));
    for (int i = which; i < indices.size(); i++)
        QFT_STEP *= makeCROT(prime(indices[which - 1], i), prime(indices[i], which - 1), i);
    return QFT_STEP;
}


// ========================================================================= //
// ------------------------------- MPO gates ------------------------------- //
// ========================================================================= //

MPO popH(SiteSet sites, int target) {
    auto ampo = AutoMPO(sites);
    ampo += sqrt(2), "Sx", target;
    ampo += sqrt(2), "Sz", target;
    return toMPO(ampo);
}

MPO popSX(SiteSet sites, int target) {
    auto ampo = AutoMPO(sites);
    ampo += 0.5 * (1 + Cplx_i), "projUp", target;
    ampo += 0.5 * (1 - Cplx_i), "S+", target;
    ampo += 0.5 * (1 - Cplx_i), "S-", target;
    ampo += 0.5 * (1 + Cplx_i), "projDn", target;
    return toMPO(ampo);
}

MPO popSY(SiteSet sites, int target) {
    auto ampo = AutoMPO(sites);
    ampo += 0.5 * (1 + Cplx_i), "projUp", target;
    ampo += -0.5 * (1 + Cplx_i), "S+", target;
    ampo += 0.5 * (1 + Cplx_i), "S-", target;
    ampo += 0.5 * (1 + Cplx_i), "projDn", target;
    return toMPO(ampo);
}

MPO popSW(SiteSet sites, int target) {
    auto ampo = AutoMPO(sites);
    ampo += 0.5 * (1 + Cplx_i), "projUp", target;
    ampo += 1 / sqrt(2), "S+", target;
    ampo += -Cplx_i / sqrt(2), "S-", target;
    ampo += 0.5 * (1 + Cplx_i), "projDn", target;
    return toMPO(ampo);
}

MPO popRAND(SiteSet sites, int target) {
    int r = rand() % 3;
    if (r == 0)
        return popSX(sites, target);
    else if (r == 1)
        return popSY(sites, target);
    else
        return popSW(sites, target);
}

MPO popCNOT(SiteSet sites, int control, int target) {
    auto ampo = AutoMPO(sites);
    ampo += 1, "projUp", control;
    ampo += 2, "projDn", control, "Sx", target;
    return toMPO(ampo);
}

MPO popSWAP(SiteSet sites, int target1, int target2) {
    auto ampo = AutoMPO(sites);
    ampo += 1, "projUp", target1, "projUp", target2;
    ampo += 1, "S+", target1, "S-", target2;
    ampo += 1, "S-", target1, "S+", target2;
    ampo += 1, "projDn", target1, "projDn", target2;
    return toMPO(ampo);
}

MPO popCROT(SiteSet sites, int control, int target, int k) {
    auto ampo = AutoMPO(sites);
    ampo += 1, "projUp", control;
    ampo += 1, "projDn", control, "projUp", target;
    ampo += exp(Cplx_i * Pi / (1 << k)), "projDn", control, "projDn", target;
    return toMPO(ampo);
}


// ========================================================================= //
// ----------------------------- Init functions ---------------------------- //
// ========================================================================= //

IndexSet initIndices(int len) {
    vector<Index> idx;
    for (int i = 1; i <= len; i++) {
        string s = "s";
        s += to_string(i);
        idx.push_back(Index(2, s));
    }

    return IndexSet(idx);
}

ITensor initTensor(int len, string form) {
    IndexSet idx = initIndices(len);
    ITensor init = ITensor(idx);
    
    if (form == "|0..0>") {
        vector<int> nz0(len, 1);
        init.set(nz0, 1.0);
    } else if (form == "|1..1>") {
        vector<int> nz1(len, 2);
        init.set(nz1, 1.0);
    } else if (form == "|+..+>") {
        vector<int> nz0(len, 1);
        init.set(nz0, 1.0);
        for (int i = 0; i < len; i++)
            init *= makeH(idx[i]);
    } else if (form == "|-..->") {
        vector<int> nz1(len, 2);
        init.set(nz1, 1.0);
        for (int i = 0; i < len; i++)
            init *= makeH(idx[i]);
    } else if (form == "|GHZn>") {
        vector<int> nz0(len, 1);
        vector<int> nz1(len, 2);
        init.set(nz0, 1.0 / sqrt(2));
        init.set(nz1, 1.0 / sqrt(2));
    } else if (form == "|Wn>") {
        for (int i = 0; i < len; i++) {
            vector<int> nz(len, 1);
            nz.at(i) = 2;
            init.set(nz, 1.0 / sqrt(len));
        }
    } else {
        throw invalid_argument("Unknown form, please use one of the following: "
                               "'|0..0>', '|1..1>', '|+..+>', '|-..->', '|GHZn>', '|Wn>'");
    }

    return noPrime(init);
}

MPS initMPS(int len, string form) {
    auto sites = SpinHalf(len, {"ConserveQNs=", false});
    MPS init;

    if (form == "|0..0>") {
        init = MPS(InitState(sites, "Up"));
    } else if (form == "|1..1>") {
        init = MPS(InitState(sites, "Dn"));
    } else if (form == "|+..+>") {
        init = MPS(InitState(sites, "Up"));
        for (int i = 1; i <= len; i++)
            init = applyMPO(popH(sites, i), init);
    } else if (form == "|-..->") {
        init = MPS(InitState(sites, "Dn"));
        for (int i = 1; i <= len; i++)
            init = applyMPO(popH(sites, i), init);
    } else if (form == "|GHZn>") {
        init = MPS(InitState(sites, "Up"));
        init.plusEq(MPS(InitState(sites, "Dn")));
        init.normalize();
    } else if (form == "|Wn>") {
        init = MPS(InitState(sites, "Up"));
        auto xplus = AutoMPO(sites);
        for (int i = 1; i <= len; i++)
            xplus += 1, "Sx", i;
        init = applyMPO(toMPO(xplus), init);
        init.normalize();
    } else {
        throw invalid_argument("Unknown form, please use one of the following: "
                               "'|0..0>', '|1..1>', '|+..+>', '|-..->', '|GHZn>', '|Wn>'");
    }

    return init;
}


// ========================================================================= //
// ----------------------------- Other helpers ----------------------------- //
// ========================================================================= //

ITensor ContractMPS(MPS mps) {
    auto con = mps(1);
    for (int i = 2; i <= length(mps); i++)
        con *= mps(i);
    
    return con;
}
