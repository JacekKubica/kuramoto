#pragma once
#include "Series.h"
#include "Equations.h"
#include "capd/capdlib.h"
#include <iostream>
#include <vector>


namespace KS {
    capd::interval intpow(int a, int b){
        return power(capd::interval(a), b);
    }

    capd::interval lambda(int k, capd::interval nu){
        return k * k * (1 - nu * k * k);
    }

    // k is dim as in equation;
    // it is not a Galerkin projection anymore,
    // we incorporate all tilde N terms here
    Equation gal_proj_k(int k, int dim, capd::interval nu){
        auto lambda_vec = std::vector<int>(dim + 1, 0);
        lambda_vec[k - 1] = 1; 
        auto ret = Equation(lambda_vec, lambda(k, nu));
        for(int i = 1; i < k; ++i){
            auto temp = std::vector<int>(dim + 1, 0);
            temp[i - 1] += 1;
            temp[k - i - 1] += 1;
            ret = ret + Equation(temp, -k);
        }
        for(int i = 1; i + k <= dim; ++i){
            auto temp = std::vector<int>(dim + 1, 0);
            temp[i - 1] += 1;
            temp[k + i - 1] += 1;
            ret = ret + Equation(temp, 2 * k);
        }
        
        return ret;
    }

    Series tilde_Nk_finite_part(int k, int M, int s) {
        Equation val;
        Equation for_derivatives;
        std::vector<Equation> partials(2 * M + 1);
        Equation hsum;
        Equation last_partial;

        if(k > 2 * M) 
            return Series(Equation(), Equation(), Equation(std::vector<int>(2 * M + 1, 0), 0).partials());
        if(k < 1 || k > 2 * M + 1) throw std::logic_error("wrong k");
        for(int i = 2 * M - k + 1; i <= M; ++i) { // tilde N terms,
                                                        // we overestimate VS and HS this way
            auto temp = std::vector<int>(2 * M + 1, 0);
            temp[2 * M] += 1;
            partials[i - 1] = Equation(temp, 2 * k / power(capd::interval(k), s));
            hsum = hsum + partials[i - 1];
            temp[i - 1] += 1;
            val = val + Equation(temp, 2 * k / power(capd::interval(k), s)); // temp[dim] holds only constant
            temp[2 * M] = 0;
            last_partial = last_partial + Equation(temp, 2 * k);
            for(int j = i; j < M; ++j){ // this is overestimation for simplicity, because otherwise
                                    // bound of variables >2M would not be valid,
                                    // as changes might not preserve monotonicity
                                    // on the main variables; it does not change too much
                temp[j - 1] = 0; temp[j] = 1;
                last_partial = last_partial + Equation(temp, 2 * k);
            }
        }

        hsum = hsum + last_partial;
        partials[2 * M] = last_partial;

        return Series(val, hsum, partials);
    }



    Equations gal_proj(int dim, capd::interval nu){
        std::vector<Equation> v(dim + 1);
        for (size_t i = 0; i < dim; i++) {
            v[i] = gal_proj_k(i + 1, dim, nu);
        }
        return Equations(v, dim + 1);
    }

    std::vector<Equation> gal_proj_v(int dim, capd::interval nu){
        std::vector<Equation> v(dim + 1);
        for (size_t i = 0; i < dim; i++) {
            v[i] = gal_proj_k(i + 1, dim, nu);
        }
        return v;
    }

    Equation ISFgtM(int k, int M, int s){
        auto temp = std::vector<int>(2 * M + 1, 0);
        temp[2 * M] = 2; // D^2
        if(1 <= k && k <= M){
            capd::interval denom = power(capd::interval(k + M + 1), s) * power(capd::interval(M), s - 1) * (s - 1);
            return Equation(temp, 1 / denom);
        }
        if (k > M) {
            capd::interval denom = power(capd::interval(k), s - 2) * power(capd::interval(M), s - 1) * (s - 2);
            return Equation(temp, 1 / denom);
        }
        throw std::logic_error("wrong k");
    }

    Equation FSFgtM(int k, int M, int s){
        if (1 <= k && k <= 2 * M) return Equation(); // zero, as those are incorporated into main modes
        auto temp = std::vector<int>(2 * M + 1, 0);
        temp[2 * M] = 2; // D^2
        if(k > 2 * M){
            capd::interval coeff = 1 / power(capd::interval(k), s - 2) * 
                                   (intpow(4, s) / intpow(2 * M + 1, s - 1) +
                                    intpow(2, s) / (intpow(M, s) * (s - 1)));
            return Equation(temp, coeff);
        }
        throw std::logic_error("wrong k");
    }

    Equation HSNkgtM(int k, int M, int s){
        auto temp = std::vector<int>(2 * M + 1, 0);
        temp[2 * M] = 1; // D^1
        if(1 <= k && k < 2 * M + 1){
            return Equation(temp, 4 * k / (intpow(M, s - 1) * (s - 1)));
        }
        if (k >= 2 * M + 1){
            return Equation(temp, 6 * k / (intpow(M, s - 1) * (s - 1)));
        }
        throw std::logic_error("wrong k");
    }

    Series NkgtM(int k, int M, int s){
        return Series(ISFgtM(k, M, s) + FSFgtM(k, M, s), HSNkgtM(k, M, s), Equation(std::vector<int>(2 * M + 1, 0), 0).partials()); // all partials = 0
    }


    Equation tildeNk_val(int k, int M, int s){
        if(1 <= k && k <= 2 * M) return Equation();
        if(k > 2 * M){
            Equation ret;
            for(int i = 0; i < M; ++i){
                auto temp = std::vector<int>(2 * M + 1, 0);
                temp[2 * M] = 1; // D^1
                temp[i] = 1;
                ret = ret + Equation(temp, (intpow(2, s) + 2) / (intpow(k, s - 2) * (s - 1)));
            }
            return ret;
        }
        throw std::logic_error("wrong k");
    }

    Equation tildeNkHS(int k, int M, int s){
        if(1 <= k && k <= 2 * M) return Equation();
        if(k > 2 * M){
            Equation ret;
            for(int i = 0; i < M; ++i){
                auto temp = std::vector<int>(2 * M + 1, 0);
                temp[i] = 1;
                ret = ret + Equation(temp, 4 * k);
            }
            auto temp = std::vector<int>(2 * M + 1, 0);
            temp[2 * M] = 1;
            ret = ret + Equation(temp, (intpow(2, s) + 2) / (intpow(k, s - 2) * (s - 1)));
            return ret;
        }
        throw std::logic_error("wrong k");
    }


    Series tildeNk(int k, int M, int s){
        std::vector<Equation> partials(2 * M + 1);
        return Series(tildeNk_val(k, M, s), tildeNkHS(k, M, s), partials);
    }

    Series F_remainder(int k, int M, int s) {
        return NkgtM(k, M, s) + tilde_Nk_finite_part(k, M, s) + tildeNk(k, M, s);
    }

    std::vector<Series> F_remainders(int M, int s){
        std::vector<Series> v;
        for(int i = 0; i < 2 * M + 1; ++i) {
            v.push_back(F_remainder(i + 1, M, s));
        }
        return v;
    }

    Equation VStildeNk(int k, int M, int s){
        if(1 <= k && k <= M){
            auto temp = std::vector<int>(2 * M + 1, 0);
            temp[2 * M] = 1;
            return Equation(temp, ((2 / intpow(2 * M - k - 1, s - 2) + 1 / intpow(2 * M + k - 1, s - 2)) / (s - 1)));
        }
        if(k > M){
            // we overestimate here to not differentiate between k <= 2M and k > 2M
            Equation ret;
            for(int i = 0; i < M; ++i){
                auto temp = std::vector<int>(2 * M + 1, 0);
                temp[i] = 1;
                ret = ret + Equation(temp, 4 * k);
            }
            return ret;
        }
        throw std::logic_error("wrong k");
    }

    Equation VSNkgeqM(int k, int M, int s){
        if(1 <= k && k <= M) return Equation();
        if(M < k && k < 2 * M + 1){
            auto temp = std::vector<int>(2 * M + 1, 0);
            temp[2 * M] = 1;
            return Equation(temp, 2 / (s - 1) * (1 / intpow(M, s - 2) + k / intpow(M, s - 1) + 1 / intpow(k, s - 2)));
        }
        if(k >= 2 * M + 1){ // it is a bound for > 2M + 1 but it is bounds from above one to 2M + 1 and
                            // for this k we want already tail bound
            auto temp = std::vector<int>(2 * M + 1, 0);
            temp[2 * M] = 1;
            return Equation(temp, 2 / (s - 1) * (1 / intpow(M, s - 2) + 2 * k / intpow(M, s - 1) + 1 / intpow(k, s - 2)));
        }
        throw std::logic_error("wrong k");
    }

    std::vector<Equation> VS_bounds(int M, int s){
        std::vector<Equation> v;
        for(int i = 0; i < 2 * M + 1; ++i) {
            v.push_back(VStildeNk(i + 1, M, s) + VSNkgeqM(i + 1, M, s));
        }
        return v;
    }

    Equations fullKS(int M, int s, capd::interval nu) {
        return Equations(gal_proj_v(2 * M, nu), F_remainders(M, s), VS_bounds(M, s), 2 * M + 1);      
    }
};