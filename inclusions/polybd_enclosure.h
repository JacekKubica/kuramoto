#pragma once
#include "diss_enclosure.h"


////////////////////////////////////////////////////////////////////////////////
///
/// DATATYPES
///
////////////////////////////////////////////////////////////////////////////////

struct PolyBd {
    Vector mainModes;
    size_t s;
    Scalar C;
    size_t M() const { return mainModes.dimension(); }
};

struct PolyBdVectorField {
    VectorField vf;
    VectorField perturbation;

    Vector operator()(const Vector &v) {
        perturbation; // = calcPerturbation(v); // TODO
        return vf(v) + perturbation(v);
    }

    Scalar C;
    size_t s;
};

struct PolyBdSolver {
    PolyBdVectorField polybd_vf;

    Scalar getStep() {
        return 1e-4;
    }
};

struct PolyBdConstants {
    PolyBdConstants(const PolyBd &pb) : C(pb.C), s(pb.s) {}
    PolyBdConstants(Scalar C, size_t s) : C(C), s(s) {}
    Scalar C;
    size_t s;
};


////////////////////////////////////////////////////////////////////////////////
///
/// FUNCTIONS USED IN IMPLEMENTATION 
///
////////////////////////////////////////////////////////////////////////////////

PolyBdConstants bPbC(const PolyBd &pb);
Scalar memoizedPow(unsigned, unsigned);
Scalar get(const PolyBdConstants &pbc, size_t i);
bool operator<=(const PolyBdConstants& a, const PolyBdConstants &b);
size_t L(const PolyBdConstants &a, const PolyBdConstants &b);


////////////////////////////////////////////////////////////////////////////////
///
/// IMPLEMENTATION
///
////////////////////////////////////////////////////////////////////////////////

PolyBd enclose(Scalar t, const PolyBd &X, PolyBdVectorField &vf, PolyBdConstants &out_pbc);