#include <iostream>

template <int m, int M>
ConeConditionVerifier<m, M>::ConeConditionVerifier(capd::interval C, int s, capd::interval niu, capd::interval eps, const capd::IVector &V, const capd::IVector &Q, const capd::IMatrix &A):
    C(C), s(s), niu(niu), eps(eps), V(V), Q(Q), A(A), Ainv(capd::matrixAlgorithms::inverseMatrix(A)), estimations(Estimations<m, M>(C, s, niu, V, A))
{
    auto Anew = capd::IMatrix(m + 1, m + 1);
    for(int i = 1; i <= m; ++i)
        for(int j = 1; j <= m; ++j)
            Anew[i][j] = capd::abs(A[i][j]) * j;
}

template <int m, int M>
bool ConeConditionVerifier<m, M>::verify(){
    std::cout.precision(17);
    const int MAX_ITER = 100;
    for(int i = 1; i <= M; ++i){
        std::cout << Gamma(i) << '\n';
        if(Gamma(i).leftBound() < eps) return false;
    }
    for(int i = M + 1; i < MAX_ITER; ++i){
        std::cout << Gamma(i) << '\n';
        if(f(i) < 2 * (capd::interval(niu) * i * i * i - i)) return true;
        if(Gamma(i).leftBound() < eps) return false;
    }
    return false;
}

template <int m, int M>
capd::interval ConeConditionVerifier<m, M>::calc_eps(){
    auto ret = Gamma(1);
    for(int i = 2; i <= M; ++i){
        ret = min(Gamma(i), ret);
    }
    return ret;
}


template <int m, int M>
capd::interval ConeConditionVerifier<m, M>::S_ND(int i){
    capd::interval res = 0;
    if(i <= M){
        for(int j = 1; j <= M; ++j){
            if(j == i) continue;
            auto temp = Q[j] * estimations.dFnew(j, i) + Q[i] * estimations.dFnew(i, j); // if it does not work, we can be more subtle here
            res += capd::abs(temp).rightBound();
        }
    }

    if(i <= m){
        for(int k = 1; k <= m; ++k){
            res += 2 * k * capd::abs(A[i][k]) * (estimations.S(M + 1 - k) + estimations.S(M + 1 + k));
            res += 2 * capd::abs(Ainv[k][i]) * (estimations.K(M + 1, -k) + estimations.K(M + 1, k));
        }
    }
    else if(i <= M){
        res += 2 * i * (estimations.S(M + 1 - i) + estimations.S(M + 1 + i));
        res += 2 * (estimations.K(M + 1, -i) + estimations.K(M + 1, i));
    }
    else{
        auto linfty = capd::vectalg::MaxNorm<capd::IVector, capd::IMatrix>();
        auto l1 = capd::vectalg::SumNorm<capd::IVector, capd::IMatrix>();
        res += (2 * l1(Anew) + 2 * i * linfty(Ainv)) * (estimations.S(i - m) + estimations.S(i + 1));
        res += 2 * i * estimations.S(i + m + 1) + (8 * i - 2) * estimations.S(1);
        res += 2 * (estimations.K(m + 1, i) + estimations.K(1, 0));
    }
    return res;
}

template <int m, int M>
capd::interval ConeConditionVerifier<m, M>::Gamma(int i){
    if(i <= M){
        return 2 * capd::abs(estimations.dFnew(i, i)).leftBound() - S_ND(i);
    }
    return 2 * capd::abs(i * i * (1 - niu * i * i)) - S_ND(i);
}
