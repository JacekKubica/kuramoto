#pragma once
#include "capd/capdlib.h"
#include "DiffInclVectorField.h"

int fact(int j) {
    int ret = 1;
    for(int i = 1; i <= j; ++i) {
        ret *= i;
    }
    return ret;
}

capd::interval ksderivative(int i, int j, int deg, capd::IVector *derivatives) {
    capd::interval ret = 0;
    for(int k = 0; k < deg + 1; ++k) {
        ret += fact(deg) / fact(k) / fact(deg - k) * derivatives[k][i] * derivatives[deg - k][j];            
    }
    return ret;
}
int m = 8; int M = 2 * m;

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
    typedef capd::IMap MapType;

    GalProjVectorField(const MapType &selector, // formula for the terms left in projection
                       const MapType &perturbation_formula, // formula for the terms cut out
                       const MapType &full_map, // formula with all terms
                       const MapType &nonlinear, // formula with all nonlinear terms
                       const VectorType &linear_diagonal, // vector of diag linear coefficients
                       bool higher_order=false)
        : main_modes(selector, MapType(node_perturb, selector.dimension(), selector.dimension(), selector.dimension())),
          perturbation_formula(perturbation_formula),
          full_map(full_map),
          nonlinear(nonlinear),
          linear_diagonal(linear_diagonal),
          higher_order(higher_order) {}

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
        if(!higher_order){
            auto perturb_vector = perturbation_formula(Z);
            for(size_t i = 0; i < perturb_vector.dimension(); ++i) {
                main_modes.getPerturbation().setParameter(i, perturb_vector[i] - perturb_vector[i].mid());
                main_modes.getSelector().setParameter(i, perturb_vector[i].mid());
            }
            return;
        }
        // else
        capd::IVector vs[32];
        for(int i = 0; i < 32; ++i) {
            vs[i] = capd::IVector(full_map.dimension());
        }
        const int DEG = 5;
        full_map.setOrder(DEG);
        full_map.computeODECoefficients(vs, DEG);   
        
        capd::IVector values(m);
        capd::IVector perturbs(m);
        double h = 1e-4;
        double hpow = 1;
        for(int deg = 0; deg < DEG - 1; ++deg){
            capd::interval temp = 0;
            for(int i = 0; i < m; ++i){
                perturbs[i] = 0;
                int k = i + 1;
                for(int j = m - k; j + k < M; ++j){
                    values[i] += 2 * k * ksderivative(j, k + j, deg, vs) / factorial(deg) * hpow;
                    perturbs[i] += 2 * k * ksderivative(j, k + j, deg + 1, vs) / factorial(deg + 1) * capd::interval(0, hpow * h);
                }
            }
            hpow *= h;
        }
        for(size_t i = 0; i < perturbs.dimension(); ++i) {
            main_modes.getPerturbation().setParameter(i, perturbs[i] - perturbs[i].mid());
            main_modes.getSelector().setParameter(i, perturbs[i].mid());
            main_modes.getSelector().setParameter(i, main_modes.getSelector().getParameter(i) + values[i]);
        }
    }
protected:
    DiffInclVectorField main_modes;
    MapType full_map, perturbation_formula, nonlinear;
    VectorType linear_diagonal;
    bool higher_order;
};

