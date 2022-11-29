#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

// ITensor gates
ITensor makeH(Index const& s);
ITensor makeCROT(Index const& s, Index const& t, int k);
ITensor makeQFT_STEP(IndexSet indices, int which);

// MPO gates
MPO popH(SiteSet sites, int target);
MPO popRAND(SiteSet sites, int target);
MPO popCNOT(SiteSet sites, int control, int target);
MPO popSWAP(SiteSet sites, int control, int target);
MPO popCROT(SiteSet sites, int control, int target, int k);

// Init methods
ITensor initTensor(int len, string form);
MPS initMPS(int len, string type);

// Other helpers
ITensor ContractMPS(MPS mps);
