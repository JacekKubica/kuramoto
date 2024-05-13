#pragma once
#include "capd/capdlib.h"
#include "DiffInclVectorField.h"


void node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i) {
        out[i] = param[i];
    }
}

class GalProjVectorField {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;

    GalProjVectorField(const MapType &selector, // formula for the terms left in projection
                       const MapType &perturbation_formula, // formula for the terms cut out
                       const MapType &full_map, // formula with all terms
                       const MapType &nonlinear, // formula with all nonlinear terms
                       const VectorType &linear_diagonal) // vector of diag linear coefficients
        : main_modes(selector, MapType(node_perturb, selector.dimension(), selector.dimension(), selector.dimension())),
          full_map(full_map),
          linear_diagonal(linear_diagonal) {}

    VectorType operator()(ScalarType t, const VectorType &v) {
        return full_map(t, v);
    }

    MatrixType derivative(ScalarType t, const VectorType &v) {
        full_map.derivative(t, v);
    }

    ScalarType getCurrentTime() { 
        throw std::logic_error("not implemented");
    }

    MapType& getSelector() { return main_modes.getSelector(); }
    MapType& getPerturbation() { return main_modes.getPerturbation(); }
    const VectorType& getLinearDiagonal() {
        return linear_diagonal;
    }
    const MapType& getNonlinear() {
        return nonlinear;
    }
    void setPerturbation(const VectorType &Z) {
        auto perturb_vector = perturbation_formula(Z);
        for(size_t i = 0; i < perturb_vector.dimension(); ++i) {
            main_modes.getPerturbation().setParameter(i, perturb_vector[i]);
        }
    }
private:
    DiffInclVectorField main_modes;
    MapType full_map, perturbation_formula, nonlinear;
    VectorType linear_diagonal;
    
};
