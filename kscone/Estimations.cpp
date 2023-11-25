template <int m, int M>
capd::interval Estimations<m, M>::S(int l){
    capd::interval res = 0;
    for(int i = l; i <= M; ++i) res += capd::max(capd::abs(V[i].leftBound()), V[i].rightBound());
    return res + C * tail_approx(capd::max(M, l - 1), s - 1);
}

template<int m, int M>
capd::interval Estimations<m, M>::K(int l, int n){
    capd::interval res = 0;
    for(int i = l; i <= M - n; ++i) res += capd::max(capd::abs(V[i + n].leftBound()), V[i + n].rightBound());
    auto r = capd::max(l, M - n + 1);
    return res + C * (tail_approx(r + n - 1, s - 2) - tail_approx(r + n, s - 1));
}

template <int m, int M>
capd::interval Estimations<m, M>::dF(int i, int j){
    if(i > j) return 2 * i * (-a(i - j) + a(i + j));
    else if(i == j) return i * i * (1 - niu * i * i) + 2 * i * a(2 * i);
    else if(i < j) return 2 * i * (a(j - i) + a(j + i));
}

template <int m, int M>
capd::interval Estimations<m, M>::dFnew(int i, int j){
    capd::interval res = 0;
    if(i <= m && j <= m){
        for(int k = 1; k <= m; ++k){
            for(int l = 1; l <= m; ++l){
                res += A[i][k] * dF(k, l) * Ainv[l][j];
            }
        }
    }
    else if(i <= m && j > m){
        for(int k = 1; k <= m; ++k){
            res += A[i][k] * dF(k, j);
        }
    }
    else if(i > m && j <= m){
        for(int l = 1; l <= m; ++l){
            res += dF(i, l) * Ainv[l][j];
        }
    }
    else
        res = dF(i, j);

    return res;
}
