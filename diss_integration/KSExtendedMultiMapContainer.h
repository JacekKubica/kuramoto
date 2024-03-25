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
        out[i] = capd::autodiff::Node(0);
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


template <int M>
class IPbMap {
public:
    using Pb = capd::pdes::PolynomialBound<capd::interval, unsigned, M>;
    typedef capd::IMatrix MatrixType;
    typedef capd::IVector VectorType; 
    IPbMap(const capd::IMap &f) : f(f) {}
    IPbMap(capd::IMap &&f) : f(f) {}
    Pb operator()(const Pb &pb) {
        return Pb(M, pb.C(), pb.exponent(), f(pb.mainModes).begin()));
    }  
private:
    capd::IMap f;
};

// ASSUMPTION: it contains data from only one map
template <int m, int M>
class KSExtendedMultiMapContainer : public ExtendedMultiMapContainer<capd::IMap> {
public:
    typedef ExtendedMultiMapContainer<IPbMap<M>> Base;
    typedef typename Base::ScalarType ScalarType;
    typedef typename Base::VectorType VectorType;
    typedef capd::IMap MapType;
    typedef typename MapType::MatrixType MatrixType;
    typedef typename VectorType::size_type size_type; 
    KSExtendedMultiMapContainer(capd::interval niu,
                                const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb,
                                bool only_projection=false);
    
    void recalcPerturbations(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb);
private:
    KSDissipativeOperator<m, M> operator;
    bool only_projection;
};

template <int m, int M>
void KSExtendedMultiMapContainer<m, M>::recalcPerturbations(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) {
    KSDissipativeOperator<m, M> ksDis(niu);
    for(int i = 0; i < m; ++i){
        auto pe = ksDis.projection_error(i, pb);
        auto mid = pe.mid();
        if(only_projection) { pe = 0; mid = 0; }
        perturb.setParameter(i, pe - mid);
        projection.setParameter(i + 1, mid);
        nonlinear.setParameter(i, pe); // ???
        //nonlinear.setParameter(i, (i + 1) * (2 * ksDis.IS(i, pb) - ksDis.FS(i, pb)));
        // if(only_projection) nonlinear.setParameter(i, 0);

        // for (size_t i = 0; i < M; i++)
        // {
        //     COUT("HERE");
        //     COUT(2 * ksDis.IS(i, pb) - ksDis.FS(i, pb));
        //     COUT(ksDis.projection_error(i, pb));
        // }  
    }

}

template <int m, int M>
KSExtendedMultiMapContainer<m, M>::KSExtendedMultiMapContainer(
                                  capd::interval niu,
                                  const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb,
                                  bool only_projection)
    : Base(MapType(node_full, m, m, m + 1, 4), // 4 is maxDerOrder
           MapType(node_parameter_placeholder, m, m, m),
           MatrixType(m, m),
           MapType(node_nonlinear, m, m, m)), // nonlinear has the same perturbation as full
           only_projection(only_projection),
           operator(niu)
{   
    projection.setParameter(0, niu);
    linear[i][i] = ksDis.lambda(i);
    recalcPerturbations(pb);
    // linear.setParameter(0, niu);
}