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
          perturbation_formula(perturbation_formula),
          full_map(full_map),
          nonlinear(nonlinear),
          linear_diagonal(linear_diagonal) {}

    VectorType operator()(ScalarType t, const VectorType &v){
        return v.dimension() == full_map.dimension() ? full_map(t, v) : main_modes(t, v);
    }

    MatrixType derivative(ScalarType t, const VectorType &v){
        return v.dimension() == full_map.dimension() ? full_map.derivative(t, v) : main_modes.derivative(t, v);
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
        COUT(perturb_vector);
        for(size_t i = 0; i < perturb_vector.dimension(); ++i) {
            main_modes.getPerturbation().setParameter(i, perturb_vector[i] - perturb_vector[i].mid());
            main_modes.getSelector().setParameter(i, perturb_vector[i].mid());
        }
    }
protected:
    DiffInclVectorField main_modes;
    MapType full_map, perturbation_formula, nonlinear;
    VectorType linear_diagonal;
};

