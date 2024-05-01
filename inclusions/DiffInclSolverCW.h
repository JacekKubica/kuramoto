#pragma once
#include "DiffInclSolver.h"

class DiffInclSolverCW : public MDiffInclSolver
{
public:
    typedef MDiffInclSolver::ScalarType ScalarType;
    typedef MDiffInclSolver::VectorType VectorType;
    typedef MDiffInclSolver::MatrixType MatrixType;
    typedef MDiffInclSolver::MapType MapType;
    typedef MDiffInclSolver::EnclosurePolicy EnclosurePolicy;
    typedef capd::IMaxNorm Norm;
    
    DiffInclSolverCW(const DiffInclVectorField &vf) : MDiffInclSolver(vf) {}
    virtual VectorType perturbationsEffects(ScalarType currentTime,
                                            const VectorType &initial) override;
};

#include "DiffInclSolverCW.cpp"
