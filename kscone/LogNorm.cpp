
template <int m, int M>
void LogNorm<m, M>::compute(){
    std::cout.precision(17);
    for(int i = 1; i <= M; ++i){
        std::cout << Gamma(i) << '\n';
    }
    for(int i = M + 1; i < 20; ++i){
        std::cout << Gamma(i) << '\n';
    }
}


template <int m, int M>
capd::interval LogNorm<m, M>::max_lognorm(){
    auto maximal = Gamma(1);
    for(int i = 2; i <= M; ++i){
        maximal = max(Gamma(i), maximal);
    }
    return maximal;
}

template <int m, int M>
capd::interval LogNorm<m, M>::S_ND(int i){
    capd::interval res = 0;
    if(i <= M){
        for(int j = 1; j <= M; ++j){
            if(j == i) continue;
            res += capd::abs(estimations.dFnew(i, j)).rightBound();
        }
    }

    if(i <= m){
        for(int k = 1; k <= m; ++k){
            res += 2 * k * capd::abs(A[i][k]) * (estimations.S(M + 1 - k) + estimations.S(M + 1 + k));
        }
    }
    else if(i <= M){
        res += 2 * i * (estimations.S(M + 1 - i) + estimations.S(M + 1 + i));
    }
    else{
        auto linfty = capd::vectalg::MaxNorm<capd::IVector, capd::IMatrix>();
        res += 2 * i * linfty(Ainv) * (estimations.S(i - m) + estimations.S(i + 1));
        res += 2 * i * estimations.S(i + m + 1) + 4 * i * estimations.S(1);
    }
    return res;
}

template <int m, int M>
capd::interval LogNorm<m, M>::Gamma(int i){
    if(i <= M){
        // COUT(estimations.dFnew(i, i));
        return estimations.dFnew(i, i).rightBound() + S_ND(i);
    }
    return i * i * (1 - niu * i * i) + S_ND(i);
}
