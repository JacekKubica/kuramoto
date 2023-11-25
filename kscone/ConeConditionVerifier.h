#pragma once
#include "capd/capdlib.h"
#include "Estimations.h"


template <int m, int M>
class ConeConditionVerifier{
public:
    ConeConditionVerifier(capd::interval C, int s, capd::interval niu, capd::interval eps, const capd::IVector &V, const capd::IVector &Q, const capd::IMatrix &A);

    bool verify();
    interval calc_eps();
private:
    capd::interval C, niu, eps;
    int s;
    capd::IVector V;
    capd::IVector Q;
    capd::IMatrix A, Ainv;
    capd::IMatrix Anew; // will depend on A
    Estimations<m, M> estimations;

    capd::interval S_ND(int i);

    capd::interval Gamma(int i);

    capd::interval f(int i){
        return (S_ND(i) + 2 * estimations.S(1)) / i;
    }
};

#include "ConeConditionVerifier.cpp"
