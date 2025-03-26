#pragma once
#include "capd/capdlib.h"

class AbsDiffInclVecField {
};

// implementuje ^, addytywne
class DiffInclVectorField {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;

    DiffInclVectorField(MapType &selector, MapType &perturbation) : selector(selector), perturbation(perturbation) {}

    VectorType operator()(ScalarType t, const VectorType &v) {
        return selector(t, v) + perturbation(t, v);
    }

    MatrixType derivative(ScalarType t, const VectorType &v) {
        return selector.derivative(t, v) + perturbation.derivative(t, v);
    }

    ScalarType getCurrentTime() { 
        throw std::logic_error("not implemented");
    }

    MapType& getSelector() { return selector; }
    MapType& getPerturbation() { return perturbation; }
    virtual void calcPerturbations(const VectorType &Z) {
        // default implementation doesn't recalculate perturbations
    }
private:
    MapType selector, perturbation;
};
