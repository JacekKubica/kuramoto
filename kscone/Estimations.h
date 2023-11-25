#pragma once
#include "capd/capdlib.h"


template <int m, int M>
class Estimations{
public:
    Estimations(capd::interval C, int s, capd::interval niu, const capd::IVector &V, const capd::IMatrix &A)
        : C(C), s(s), niu(niu), V(V), A(A), Ainv(capd::matrixAlgorithms::inverseMatrix(A)){}

    capd::interval S(int l);

    capd::interval K(int l, int n);

    capd::interval dF(int i, int j);
    capd::interval dFnew(int i, int j);

void print(){
    for (size_t i = 1; i <= M; i++) {
        for (size_t j = 1; j <= M; j++) {
            std::cout << dFnew(i, j) << ' ';
        }
        std::cout << '\n';
    }
}

private:
    capd::interval C, niu;
    int s;
    capd::IVector V;
    capd::IMatrix A, Ainv;

    capd::interval tail_approx(double l, int e){
        return capd::interval(1) / (e * (capd::interval(l) ^ e));
    }

    capd::interval mode(int k){
        auto max = (C / (capd::interval(k) ^ s)).rightBound();
        return capd::interval(-max, max);
    }

    capd::interval a(int k){
        return (k <= M) ? V[k] : mode(k);
    }
};

#include "Estimations.cpp"
