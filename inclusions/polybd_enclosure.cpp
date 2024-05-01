#include "polybd_enclosure.h"
#include "pde_bounds.h"

Scalar memoizedPow(unsigned, unsigned){
    // I think we can use array cache here easily
}

Scalar get(const PolyBdConstants &pbc, size_t i) {
    return pbc.C / memoizedPow(i + 1, pbc.s);
}

PolyBdConstants bPbC(const PolyBd &pb) {
    auto ret = nonlinearPbC(pb);
    return PolyBdConstants { ret.C / V(pb.M()), ret.s + p};
}

size_t L(const PolyBdConstants &a, const PolyBdConstants &b);

bool operator<=(const PolyBdConstants& a, const PolyBdConstants &b) {
    return true;
}
bool operator<(const PolyBdConstants& a, const PolyBdConstants &b) {
    return true;
}

PolyBd enclose(Scalar t, const PolyBd &X, PolyBdVectorField vf, PolyBdConstants &out_pbc) {
    vf.C = X.C; vf.s = X.s;
    PolyBd Z; Z.C = X.C; Z.s = X.s; // TODO intial guess 
    bool validated = false;
    while(!validated) {
        Z.mainModes = enclose(t, X.mainModes, vf.vf);
        validated = true;
        auto b = bPbC(Z);
        auto g = exponentialPbC(Z);
        if(b.s > X.s) {
            for(size_t i = X.M(); i < L(b, X); ++i) {
                if(get(g, i) > get(Z, i)){
                    validated = false;
                    break;
                }
            }
        } else if (b.s == X.s) {
            validated = !(b <= Z);
        } else { // b.s < X.S
            size_t L_ = L(X, b);
            size_t p = X.M() > L_ ? X.M() : L_;
            validated = b.s >= Z.s && get(b, p) <= get(Z, p);
        }
        

        // TODO ponder 207-209; compare with what is done for near modes and my kinda natural candidate
        if(!validated) {
            Z.C = ; //TODO
        }
        else if(b < X) {
            Z.C = X.C * memoizedPow(X.M() + 1, Z.s - X.s);
        } 
        // TODO decpower
    }
    // TODO refinement, see also 231
    // TODO set out pbc
}