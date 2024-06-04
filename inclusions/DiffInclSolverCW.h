#pragma once
#include "DiffInclSolver.h"

template <typename VectorField=DiffInclVectorField,
          typename EnclosurePolicy=DiffInclusionEnclosurePolicy>
class DiffInclSolverCW : public MDiffInclSolver<VectorField, EnclosurePolicy>{
public:
    using Base = MDiffInclSolver<VectorField, EnclosurePolicy>;
    typedef typename Base::ScalarType ScalarType;
    typedef typename Base::VectorType VectorType;
    typedef typename Base::MatrixType MatrixType;
    typedef typename Base::MapType MapType;
    typedef capd::IMaxNorm Norm;
    
    DiffInclSolverCW(const VectorField &vf) : Base(vf) {}
    virtual VectorType perturbationsEffects(ScalarType currentTime,
                                            const VectorType &initial,
                                            const VectorType &W1,
                                            const VectorType &W2) override;
};

#include "DiffInclSolverCW.cpp"
