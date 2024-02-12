
// TODO should we take right bounds for all those constants?

// D_1(k<2M) = \max \set{(k + 1)^(s - 1) \abs{FS(k)} \mid M \leq k < 2M}
template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::D1ksmall(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    capd::interval res = 0;
    for(int k = M; k < 2 * M; ++k){
        res = capd::max(res, take_power(k + 1, pb.exponent() - 1) * capd::abs(KSDissipativeOperator<m, M>::FS(k, pb)));
    }
    return res.rightBound();
}

// D_1(k \geq 2M) = C (\frac{2^{s+1}}{2M + 1} \sum_{n = 0}^{M - 1} \abs{a_n}
// + \frac{C4^s}{(2M + 1)^{s + 1}} + \frac{C2^s}{(s - 1)M^s})
template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::D1klarge(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    capd::interval res = 0;
    auto s = pb.exponent();
    auto C = pb.C();
    for(int n = 0; n < M; ++n){
        res += capd::abs(pb[n]);
    }
    res *= take_power(2, s + 1) / (2 * M + 1);
    res += C * take_power(4, s) / take_power(2 * M + 1, s + 1);
    res += C * take_power(2, s) / ((s - 1) * take_power(M, s));
    return (C * res).rightBound();
}


template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::D1(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    return capd::max(KSDissipativeOperator<m, M>::D1ksmall(pb), KSDissipativeOperator<m, M>::D1klarge(pb)).rightBound();
}


template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::FS(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    capd::interval res = 0;
    auto s = pb.exponent();
    auto C = pb.C();
    // FS(k) = \sum_{n = 0}^{k-1} a_n a_{k-n-1}
    if(k < M){
        for(int n = 0; n < k; ++n){
            res += pb[n] * pb[k - 1 - n];
        }
    }
    // FS(k) \subset 2 \sum_{k-M \leq n < \frac{k-1}{2}} a_n a_{k - n - 1} + o(k) a_{\frac{k - 1}{2}}^2
    // + 2C \sum_{n = 0}^{k - M - 1} \frac{\abs{a_n}}{(k-n)^s} \interval{-1, 1}
    else if(k >= M && k < 2 * M){
        for(int n = k - M; n < (k - 1) / 2; ++n){
            res += 2 * pb[n] * pb[k - 1 - n];
        }
        res += (k % 2) * pb[(k - 1) / 2] * pb[(k - 1) / 2];
        for(int n = 0; n < k - M; ++n){
            res += 2 * C * capd::abs(pb[n]) / take_power(k - n, s) * capd::interval(-1, 1);
        }
    }

    // FS(k) = \frac{C}{(k + 1)^{s-1}}(\frac{2^{s+1}}{2M + 1} \sum_{n = 0}{M - 1} \abs{a_n}
    // + \frac{C4^s}{(2M + 1)^{s + 1}} + \frac{C2^s}{(s - 1)M^s}) \interval{-1, 1}
    else if(k >= 2 * M){
        res = KSDissipativeOperator<m, M>::D1klarge(pb);
        res /= take_power(k + 1, s - 1);
        res *= capd::interval(-1, 1);
    }
    return res;
}


// D_2 =  \frac{C}{M + 1} (\frac{C}{(M + 1)^{s - 1} (s - 1)} + \sum_{n = 0}^{M - 1} \abs{a_n})
template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::D2(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    capd::interval res = 0;
    auto s = pb.exponent();
    auto C = pb.C();
    res += C / (take_power(M + 1, s - 1) * (s - 1));
    for(int n = 0; n < M; ++n){
        res += capd::abs(pb[n]);
    }
    res *= C / (M + 1);
    return res.rightBound();
}


template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::IS(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb)  const{
    capd::interval res = 0;
    auto s = pb.exponent();
    auto C = pb.C();
    // IS(k) \subset \sum_{n = 0}^{M - k - 2} a_n a_{n + k + 1} + C \sum_{n = M - k -1 }^{M-1} \frac{\abs{a_n}}{(k+2+n)^s} \interval{-1, 1}
    // + \frac{C^2}{(k + M + 2)^s (s-1) M^{s-1}} \interval{-1, 1}
    if(k < M){
        for(int i = 0; i < M - k - 1; ++i) res += pb[i] * pb[i + k + 1];
        for(int i = M - k - 1; i < M; ++i){
            res += C * capd::abs(pb[i]) / take_power(k + i + 2, s) * capd::interval(-1, 1);
        }
        res += C * C / (take_power(k + M + 2, s) * (s - 1) * take_power(M, s - 1)) * capd::interval(-1, 1);
        return res;
    }
    //IS(k) \subset \frac{C}{(k + 1)^{s - 1} (M + 1)} (\frac{C}{(M + 1)^{s - 1} (s - 1)} + \sum_{n = 0}^{M - 1} \abs{a_n}) \interval{-1, 1}
    if(k >= M){ // TODO probably no need for this one
        res = D2(pb);
        res /= take_power(k + 1, s - 1);
        res *= capd::interval(-1, 1);
        return res;
    }
}

template <int m, int M>
capd::pdes::PolynomialBound<capd::interval, int, M> KSDissipativeOperator<m, M>::N(const capd::IVector &z) const{
    auto T_plus_main = T;
    for(int i = 0; i < z.dimension(); ++i) T_plus_main[i] = z[i];
    capd::pdes::PolynomialBound<capd::interval, int, M> res;
    for(int i = 0; i < M; ++i){
        res[i] = (i + 1) * (2 * IS(i, T_plus_main) - FS(i, T_plus_main));
    }
    res.exponent() = T.exponent() - 2;
    res.C() = D1(T_plus_main) + 2 * D2(T_plus_main);
    return res;
}

template <int m, int M>
capd::pdes::PolynomialBound<capd::interval, int, M> KSDissipativeOperator<m, M>::guess_far_tail(const capd::IVector &x, const capd::pdes::PolynomialBound<capd::interval, int, M> &T0) const{
    capd::pdes::PolynomialBound<capd::interval, int, M> res;

    // res.C() = 3 * T0.C();
    // res.exponent() = T0.exponent();

    capd::interval A = 0;
    for(int i = 0; i < m; ++i){
        A += capd::abs(x[i]).rightBound();
    }
    for(int i = m; i < M; ++i){
        A += capd::abs(T0[i]).rightBound();
        res[i] = T0[i];
    }
    if(T0.C() == 0){
        res.C() = capd::interval(25);
        res.exponent() = 10;
        return res;
    }

    res.C() = capd::max(3 * D1ksmall(T0) * (2 * M + 1) / (take_power(2, T0.exponent() + 1) * A), 1.01 * T0.C()); // TODO idk if this 1.01 should be here
    int s_max = 1;
    while((M + 1) * (M + 1) * (2 * M + 1) * (niu - 1 / ((M + 1) * (M + 1))) / A >= take_power(2, s_max + 1) && s_max < T0.exponent()) ++s_max;
    res.exponent() = s_max;
    res.C() *= take_power(M + 1, s_max - T0.exponent());
    return res;
}

template <int m, int M>
capd::IVector KSDissipativeOperator<m, M>::P(const capd::IVector &x) const{
    capd::IVector res(m); for(int i = 0; i < m; ++i) res[i] = 0; //TODO is it zeroed at cstr?
    for(int i = 0; i < m; ++i){
        int k = i + 1;
        res[i] = k * k * (1 - niu * k * k) * x[i];
        for(int j = 0; j < i; ++j){
            res[i] -= k * x[j] * x[i - 1 - j];
        }
        for(int j = 0; j + k < m; ++j){
            res[i] += 2 * k * x[j] * x[j + k];
        }
    }
    return res;
}

template <int m, int M>
capd::interval KSDissipativeOperator<m, M>::projection_error(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const{
    int s = pb.exponent(); auto C = pb.C();
    capd::interval res = 0;
    for(int n = m - k - 1; n < M - k - 1; ++n){
        res += 2 * (k + 1) * pb[n] * pb[n + k + 1];
    }
    for(int i = M - k - 1; i < M; ++i){
        res += C * 2 * (k + 1) * capd::abs(pb[i]) / take_power(k + i + 2, s) * capd::interval(-1, 1);
    }
    res += 2 * (k + 1) * C * C  / (take_power(k + M + 2, s) * (s - 1) * take_power(M, s - 1)) * capd::interval(-1, 1);

    return res;
}

template <int m, int M>
void KSDissipativeOperator<m, M>::nodeP(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        out[i] = k * k * (1 - param[0] * k * k) * in[i];
        for(int j = 0; j < i; ++j){
            out[i] -= k * in[j] * in[i - 1 - j];
        }
        for(int j = 0; j + k < dimIn; ++j){
            out[i] += 2 * k * in[j] * in[j + k];
        }
        out[i] += param[i + 1];
    }
}

template <int m, int M>
void KSDissipativeOperator<m, M>::node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        out[i] = param[i];
    }
}

template <int m, int M>
void KSDissipativeOperator<m, M>::update_perturbation(const capd::IVector &z){
    auto T_plus_main = T;
    for(int i = 0; i < z.dimension(); ++i) T_plus_main[i] = z[i];
    for(int i = 0; i < m; ++i){
        auto pe = projection_error(i, T_plus_main);
        auto mid = pe.mid();
        perturb.setParameter(i, pe - mid);
        pm.setParameter(i + 1, mid);
    }
}
