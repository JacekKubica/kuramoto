#pragma once
#include "capd/capdlib.h"
#include "Estimations.h"


template <int m, int M>
class LogNorm{
public:
    LogNorm(capd::interval C, int s, capd::interval niu, const capd::IVector &V, const capd::IMatrix &A)
        : C(C), s(s), niu(niu), V(V), A(A), Ainv(capd::matrixAlgorithms::inverseMatrix(A)), estimations(Estimations<m, M>(C, s, niu, V, A))
        {}

    void compute();
    capd::interval max_lognorm();
private:
    capd::interval C, niu, eps;
    int s;
    capd::IVector V;
    capd::IMatrix A, Ainv;
    Estimations<m, M> estimations;

    capd::interval S_ND(int i);

    capd::interval Gamma(int i);
};

#include "LogNorm.cpp"
