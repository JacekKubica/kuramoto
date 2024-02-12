#pragma once
#include "capd/capdlib.h"

template <class FMapType, class GMapType = FMapType, class LMapType = typename FMapType::MatrixType, class NMapType = FMapType>
class ExtendedMultiMapContainer {
public:
    typedef ExtendedMultiMap<FMapType, GMapType> Base; // not really base
    typedef typename Base::ScalarType ScalarType;
    typedef typename Base::VectorType VectorType;
    typedef typename VectorType::size_type size_type; 
    ExtendedMultiMapContainer(FMapType projection,
                              GMapType perturb, 
                              LMapType linear, 
                              NMapType nonlinear)
        : projection(projection),
          perturb(perturb),
          linear(linear),
          nonlinear(nonlinear),
          emm(this->projection, this->perturb, this->linear, this->nonlinear) {}

    const FMapType& getProjection() const { return projection; }
    const GMapType& getPerturb() const { return perturb; }
    const LMapType& getLinear() const { return linear; }
    const NMapType& getNonlinear() const { return nonlinear; }

    VectorType diagonalCoefficients() const {
        return emm.diagonalCoefficients();
    }

    VectorType evalNonlinearPlusPerturb(const VectorType &X) const { return nonlinear(X) + perturb(X); }
    VectorType evalLinearNondiagonal(const VectorType &X) const { return emm.evalLinearNondiagonal(X); }
    VectorType operator()(const VectorType & X) const { return emm(X); }
// private:
    FMapType projection;
    GMapType perturb; 
    LMapType linear;
    NMapType nonlinear;  // we need to keep those because IMultiMap keeps pointers
    ExtendedMultiMap<FMapType, GMapType, LMapType, NMapType> emm;
};


template <class MapT>
struct ComposedMap {
    typedef MapT FunctionType;
    typedef typename capd::IMatrix MatrixType;
    typedef typename capd::IVector VectorType;
    typedef typename VectorType::size_type size_type;
    
    MapT f; capd::IMatrix A; capd::IVector V;

    size_type dimension() { return f.dimension().first; }

    ComposedMap(const MapT &f, const capd::IMatrix &A, const capd::IVector &V)
        : f(f), A(A), V(V) {}
    ComposedMap(const MapT &f, const capd::IMatrix &A)
        : f(f), A(A), shift(false) {}
    
    VectorType operator()(const VectorType & Y) const { 
        auto ret = capd::matrixAlgorithms::inverseMatrix(A) * Y;
        if(shift) ret += V;
        return A * f(ret);
        // return A * f(capd::matrixAlgorithms::inverseMatrix(A) * Y + V);
    }
private:
    bool shift = true;
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

capd::IMatrix change_of_basis_inverse(const capd::DVector &V, const capd::DMatrix &derivative, int m){
    capd::DVector eigenRealPart(m);
    capd::DVector eigenImPart(m);
    capd::DMatrix vectorRealPart(m, m);
    capd::DMatrix vectorImPart(m, m);

    capd::alglib::computeEigenvaluesAndEigenvectors(derivative, eigenRealPart, eigenImPart, vectorRealPart, vectorImPart);
    capd::IMatrix A = capd::IMatrix::Identity(m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            A[i][j] = vectorRealPart[i][j];
        }
    }
    return A;
}

ExtendedMultiMapContainer<ComposedMap<capd::IMap>, ComposedMap<capd::IMap>,
                          capd::IMatrix, ComposedMap<capd::IJet>> changeVariables(
    capd::IMatrix &A, // we return change of basis
    const capd::IVector &V,
    const ExtendedMultiMapContainer<capd::IMap>& emm,
    unsigned toLinearize
) {
    const unsigned m = A.dimension().first;
    capd::IMatrix new_linear_part(m, m);
    const int max_der_order = 4;
    capd::IJet fmz(m, m, max_der_order - 1);
    // capd::IMap f(node_no_params2, M, M, 1, max_der_order); f.setParameter(0, 0.75);
    // f(centralized(V), fMz);
    emm.getProjection()(centralized(V), fmz); // calculate the jet
    auto mp = fmz.first(1);
    for(int k = 0; k < m; ++k){
        for(int i = 0; i < m; ++i){
            if(i < toLinearize || i == k){
                new_linear_part[i][k] = fmz(mp)[i];
                fmz(mp)[i] = 0;
            }
        }
        fmz.hasNext(mp);
    }
    // COUT("!!!");
    // COUT(fmz(V - centralized(V)) + new_linear_part * (V - centralized(V)));
    // COUT(emm(V));
    // COUT("???");

    // change of basis actually return inverse
    auto Ainv = change_of_basis_inverse(centralized_double(V), centralized_double(new_linear_part), m); 
    A         = capd::matrixAlgorithms::inverseMatrix(Ainv);

    return ExtendedMultiMapContainer<ComposedMap<capd::IMap>, ComposedMap<capd::IMap>,
                          capd::IMatrix, ComposedMap<capd::IJet>>
                          (ComposedMap<capd::IMap>(emm.getProjection(), A, V),
                           ComposedMap<capd::IMap>(emm.getPerturb(), A, V),
                           A * new_linear_part * Ainv,
                           ComposedMap<capd::IJet>(fmz, A));
                        //    ComposedMap<capd::IMatrix>(new_linear_part, A, V),
                        //    ComposedMap<capd::IJet>(fmz, A, V));
}