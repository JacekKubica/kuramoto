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

    DiffInclVectorField(const MapType &selector, const MapType &perturbation) : selector(selector), perturbation(perturbation) {}

    VectorType operator()(ScalarType t, const VectorType &v) const {
        return selector(t, v) + perturbation(t, v);
    }
    VectorType mainModesCall(ScalarType t, const VectorType &v) {
        return (*this)(t, v);
    }

    MatrixType derivative(ScalarType t, const VectorType &v) {
        return selector.derivative(t, v) + perturbation.derivative(t, v);
    }

    ScalarType getCurrentTime() { 
        throw std::logic_error("not implemented");
    }

    MapType& getSelector() { return selector; }
    MapType& getPerturbation() { return perturbation; }
private:
    MapType selector, perturbation;
};
