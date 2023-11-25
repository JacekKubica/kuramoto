
#include "capd/capdlib.h"
#include <iostream>
#include "ConeConditionVerifier.h"
#include "LogNorm.h"

using namespace capd;

const double niu = 0.75;

bool verify_cone_conditions(){
    const int m = 8, M = 16, s = 16;
    const double eps = 0.0329203, C = 1.28151e007;

    IVector V(M + 1);
    V[1] = interval(-0.075, 0.075);
    V[2] = -7.031250e-004 + 7.035838e-004 * interval(-1, 1);
    V[3] = 0.000000e+000 + 1.223300e-005 * interval(-1, 1);
    V[4] = -2.248670e-008 + 6.419194e-008 * interval(-1, 1);
    V[5] = 0.000000e+000 + 5.342963e-010 * interval(-1, 1);
    V[6] = -8.851774e-013 + 2.151082e-012 * interval(-1, 1);
    V[7] = 0.000000e+000 + 1.630102e-014 * interval(-1, 1);
    V[8] = -1.661427e-017 + 6.735586e-017 * interval(-1, 1);
    V[9] = 0 + 4.19092e-019 * interval(-1, 1);
    V[10] = -4.40581e-022 + 1.59964e-021 * interval(-1, 1);
    V[11] = 9.73596e-024 * interval(-1, 1);
    V[12] = -6.79843e-027 + 2.76354e-025 * interval(-1, 1);
    V[13] = 2.8995e-023 * interval(-1, 1);
    V[14] = -1.61085e-031 + 3.20088e-021 * interval(-1, 1);
    V[15] = 2.96259e-019 * interval(-1, 1);
    V[16] = 1.30238e-017 * interval(-1, 1);


    DVector Q(M + 1);
    for(int i = 0; i <= m; ++i) Q[i] = 1;
    for(int i = m; i <= M; ++i) Q[i] = -1;

    const auto A = IMatrix::Identity(9);

    return ConeConditionVerifier<m, M>(C, s, niu, eps, V, Q, A).verify();
}


template <int m, int M>
double dF(int i, int j, const DVector &a){
    if(i > j) return -2 * i * a[i - j] + 2 * i * a[i + j];
    else if(i == j) return i * i * (1 - niu * i * i) + 2 * i * a[2 * i];
    return 2 * i * (a[j - i] + a[j + i]);
}

template <int m, int M>
IMatrix change_of_basis(const DVector &V){
    DMatrix temp(m, m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            temp[i][j] = dF<m, M>(i + 1, j + 1, V);
        }
    }

    DVector eigenRealPart(3);
    DVector eigenImPart(3);
    DMatrix vectorRealPart(3, 3);
    DMatrix vectorImPart(3, 3);

    alglib::computeEigenvaluesAndEigenvectors(temp, eigenRealPart, eigenImPart, vectorRealPart, vectorImPart);
    IMatrix A(m + 1, m + 1);
    A[0][0] = 1;
    for (size_t i = 1; i <= m; i++) {
        for (size_t j = 1; j <= m; j++) {
            A[i][j] = vectorRealPart[i-1][j-1];
        }
    }
    return capd::matrixAlgorithms::inverseMatrix(A);
}

void compute_lognorm(){
    const int m = 3, M = 10, s = 10;
    const double C = 9077.32, delta = 5.652457e-002;

    DVector _V {0.712361, -0.123239, 0.0101786};

    auto A = change_of_basis<m, M>(_V);

    IVector V(M + 1);

    V[1] = interval(_V[0] - delta, _V[0] + delta);
    V[2] = interval(_V[1] - 2.077854e-002, _V[1] + 2.077854e-002);
    V[3] = interval(_V[2] - 2.672539e-003, _V[2] + 2.672539e-003);

    V[4] = -0.000699649 + 0.00027784 * interval(-1, 1);
    V[5] = 4.20236e-005 + 2.10249e-005 * interval(-1, 1);
    V[6] = -2.34933e-006 + 1.35644e-006 * interval(-1, 1);
    V[7] = 1.24839e-007 + 8.07514e-008 * interval(-1, 1);
    V[8] = -6.41761e-009 + 4.55835e-009 * interval(-1, 1);
    V[9] = 3.21464e-010 + 4.44867e-010 * interval(-1, 1);
    V[10] = -1.58037e-011 + 8.01665e-010 * interval(-1, 1);

    LogNorm<m, M>(C, s, niu, V, A).compute();
}

int main(){
    std::cout << verify_cone_conditions() << '\n';
    compute_lognorm();
}
