#pragma once
#include "capd/capdlib.h"
#include <exception>
#include "ExtendedMultiMap.h"
#include "PolynomialBound.h"
#include "KSDissipativeOperator.h"
#include "ExtendedMultiMapContainer.h"


// for projection and perturbation construction
void node_full(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        out[i] = k * k * (1 - param[0] * k * k) * in[i];
        for(int j = 0; j < i; ++j){
            out[i] -= k * in[j] * in[i - 1 - j];
        }
        for(int j = 0; j + k < dimIn; ++j){
            out[i] += 2 * k * in[j] * in[j + k];
        }
        out[i] += param[i + 1];
    }
}
void node_nonlinear(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        for(int j = 0; j < i; ++j){
            out[i] -= k * in[j] * in[i - 1 - j];
        }
        for(int j = 0; j + k < dimIn; ++j){
            out[i] += 2 * k * in[j] * in[j + k];
        }
        out[i] += param[i];
    }
}
void node_linear(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        out[i] = k * k * (1 - param[0] * k * k) * in[i];
    }
}
void node_parameter_placeholder(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        out[i] = param[i];
    }
}


// ASSUMPTION: it contains data from only one map
template <int m, int M>
class KSExtendedMultiMapContainer : ExtendedMultiMapContainer<IMap> {
public:
    typedef ExtendedMultiMapContainer<capd::IMap> Base;
    typedef typename Base::ScalarType ScalarType;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::MapType MapType;    
    typedef typename VectorType::size_type size_type; 
    KSExtendedMultiMapContainer(capd::interval niu,
                                const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb);

    const MapType& getNonlinear() const { return emm.getNonlinear(); }

    VectorType diagonalCoefficients() const {
        return emm.diagonalCoefficients();
    }

    VectorType operator()(const VectorType & X) const { return emm(X); }
// private:
    MapType projection, perturb, linear, nonlinear;  // we need to keep those because IMultiMap keeps pointers
    ExtendedMultiMap<MapType> emm;
};


template <int m, int M>
KSExtendedMultiMapContainer<m, M>::KSExtendedMultiMapContainer(
                                  capd::interval niu,
                                  const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb)
    : Base(MapType(node_full, m, m, m + 1),
           MapType(node_parameter_placeholder, m, m, m),
           MapType(node_linear, m, m, 1),
           MapType(node_nonlinear, m, m, m))
{
    
    projection.setParameter(0, niu);
    linear.setParameter(0, niu);
    KSDissipativeOperator<m, M> ksDis(niu);
    for(int i = 0; i < m; ++i){
        auto pe = ksDis.projection_error(i, pb);
        auto mid = pe.mid();
        perturb.setParameter(i, pe - mid);
        projection.setParameter(i + 1, mid);
        nonlinear.setParameter(i, (i + 1) * (2 * ksDis.IS(i, pb) - ksDis.FS(i, pb)));
    }
}