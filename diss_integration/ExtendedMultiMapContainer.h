#pragma once
#include "capd/capdlib.h"

template <class FType, class GType = FType, class LType = FType, class NType = FType>
class ExtendedMultiMapContainer {
public:
    typedef ExtendedMultiMap<FType, GType> Base; // not really base
    typedef typename Base::ScalarType ScalarType;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::MapType MapType;    
    typedef typename VectorType::size_type size_type; 
    ExtendedMultiMapContainer(const FType &projection,
                              const GType &perturb, 
                              const Ltype &linear, 
                              const NType &nonlinear)
        : projection(projection),
          perturb(perturb),
          linear(linear),
          nonlinear(nonlinear)
          emm(projection, perturb, linear, nonlinear) {}


    const MapType& getNonlinear() const { return nonlinear; }
    const MapType& getLinear() const { return linear; }
    const MapType& getProjection() const { return projection; }
    const MapType& getPerturb() const { return perturb; }

    VectorType diagonalCoefficients() const {
        return emm.diagonalCoefficients();
    }

    VectorType operator()(const VectorType & X) const { return emm(X); }
// private:
    MapType projection, perturb, linear, nonlinear;  // we need to keep those because IMultiMap keeps pointers
    ExtendedMultiMap<MapType> emm;
};


template <class MapT>
struct ComposedMap {
    typedef typename capd::IVector VectorType;
    MapT f; capd::IMatrix A; capd::IVector V;
    ComposedMap(const MapT &f, const capd::IMatrix &A, const capd::IVector &V)
        : f(f), A(A), V(V) {}
    
    VectorType operator()(const VectorType & Y) const { 
        return A * f(capd::matrixAlgorithms::inverse(A) * Y + V);
    }
};


capd::IVector moved_to_zero(capd::IVector V) {
    for(auto &x: V) {
        x = x - x.mid();
    }
    return V;
}

capd::IVector centralized(capd::IVector V) {
    for(auto &x: V) {
        x = x.mid();
    }
    return V;
}

capd::DVector centralized_double(const capd::IVector &V) {
    capd::DVector R(V.dimension());
    for(int i = 0; i < V.dimension(); ++i) {
        R[i] = V[i].mid().leftBound();
    }
    return R;
}

capd::DMatrix centralized_double(const capd::IMatrix &A) {
    capd::DMatrix R(A.numberOfRows(), A.numberOfColumns());
    for(int i = 0; i < A.numberOfRows(); ++i){
        for(int j = 0; j < A.numberOfColumns(); ++j){
            R[i][j] = A[i][j].mid().leftBound();
        }
    }
    return R;
}


capd::IMatrix change_of_basis_inverse(const capd::DVector &V, const capd::DMatrix &Derivative, int m){
    int M = V.dimension();
    capd::DVector eigenRealPart(m);
    capd::DVector eigenImPart(m);
    capd::DMatrix vectorRealPart(m, m);
    capd::DMatrix vectorImPart(m, m);

    capd::DMatrix temp(m, m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            temp[i][j] = Derivative[i][j];
        }
    }

    capd::alglib::computeEigenvaluesAndEigenvectors(temp, eigenRealPart, eigenImPart, vectorRealPart, vectorImPart);
    capd::IMatrix A = capd::IMatrix::Identity(M);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            A[i][j] = vectorRealPart[i][j];
        }
    }
    return A;
}

ExtendedMultiMapContainer<ComposedMap<capd::IMap>, ComposedMap<capd::IMap>,
                          ComposedMap<capd::IMatrix>, ComposedMap<capd::IJet>> changeVariables(
    const capd::IMatrix &A,
    const capd::IVector &V,
    const ExtendedMultiMapContainer<capd::IMap>& emm,
    int M
) {
    capd::IMatrix new_linear_part(M, M);
    const int max_der_order = 4;
    capd::IJet fmz(M, M, max_der_order - 1);
    fm(centralized(V), fmz); // calculate the jet
    auto mp = fmz.first(1);
    for(int k = 0; k < M; ++k){
        for(int i = 0; i < M; ++i){
            new_linear_part[i][k] = fmz(mp)[i];
            fmz(mp)[i] = 0;
        }
        fmz.hasNext(mp);
    }

    // change of basis actually return inverse
    auto Ainv = change_of_basis_inverse(centralized_double(V), centralized_double(new_linear_part), M); 
    auto A    = capd::matrixAlgorithms::inverseMatrix(Ainv);

    return ExtendedMultiMapContainer(ComposedMap(emm.getProjection(), A, V),
                                     ComposedMap(emm.gerPerturb(), A, V),
                                     ComposedMap(emm.getNonlinear(), A, V),
                                     ComposedMap(emm.getLinear(), A, V))
}