#pragma once
#include "Equations.h"
#include "KS.h"

void check_bif(const std::vector<Bound> &bounds, const std::vector<Bound> &derivbounds, int m, int M,
               capd::interval mu_lower, capd::interval mu_bif, capd::interval mu_upper,
               capd::interval K, capd::interval cmuupper, capd::interval Cmuuper, int s, int p,
               capd::interval l, int omega, capd::interval gamma_minus, capd::interval gamma_plus,
               const std::vector<Bound> &derivbounds2, capd::interval zeta, capd::interval cmubif){
    
    
    int n1 = bounds[m - 1].degree;
    capd::interval alpha1 = bounds[m - 1].coeff;
    if(zeta * (1 - zeta * zeta * cmubif / cmuupper) 
        + power(K, n1) * power(KS::lambda(m, mu_upper), (n1 - 3) / capd::interval(2))
                        * power(-cmuupper, (n1 - 1) / capd::interval(2)) < 0){
                    

                    std::cout << "COND 49 SATISFIED\n";
    } else {
        std::cout << "COND 49 FAILED\n";
        COUT(zeta * (1 - zeta * zeta)  
                + power(K, n1) * power(KS::lambda(m, mu_upper), (n1 - 3) / capd::interval(2))
                               * power(-cmuupper, (n1 - 1) / capd::interval(2)));
    }


    bool cond = true;
    for(int k = m; k < M; ++k) {
        int nk = bounds[k].degree;
        capd::interval alphakbyktop = bounds[k].coeff;
        if(KS::lambda(k + 1, mu_upper) / power(capd::interval(k + 1), s)
           + alphakbyktop * power(K, nk) * power(Cmuuper, nk) /* / power(capd::interval(k + 1), p) */ < 0){
                continue;
           }
        cond = false;
        COUT(k);
        COUT(KS::lambda(k + 1, mu_upper) / power(capd::interval(k + 1), s)
           + alphakbyktop * power(K, nk) * power(Cmuuper, nk - omega));
        std::cout << "COND 50 FAILED\n";
    }
    if(cond){
        std::cout << "COND 50 SATISFIED\n";
    }

    capd::interval beta1 = derivbounds[m - 1].coeff;
    capd::interval d1 = derivbounds[m - 1].degree;
    if(-1 + beta1 * power(KS::lambda(m, mu_lower), omega - 1) < 0){
        std::cout << "COND 51 SATISFIED\n";
    } else {
        std::cout << "COND 51 FAILED\n";
    }

    cond = true;
    for(int k = m; k < M; ++k) {
        int nk = bounds[k].degree;
        capd::interval betaktimesk = derivbounds[k].coeff;
        if(2 * KS::lambda(k + 1, mu_bif)
           - betaktimesk * power(capd::abs(KS::lambda(m, mu_lower)), d1) < -l){
                continue;
           }
        cond = false;
        std::cout << "COND 53 FAILED\n";
        COUT(k);
        COUT(betaktimesk * capd::abs(KS::lambda(m, mu_lower)));
    }
    if(cond) {
        std::cout << "COND 53 SATISFIED\n";
    }


    capd::interval beta12 = derivbounds2[m - 1].coeff;
    capd::interval d12 = derivbounds2[m - 1].degree;
    if((2 - 6 * gamma_minus * gamma_minus) - beta12 * capd::abs(power(KS::lambda(m, mu_upper), (d12 - capd::interval(2)) / 2)
                    / power(-cmuupper, capd::interval(d12) / 2)) > 0)
    {
        std::cout << "COND 54 SATISFIED\n";
    } else {
        COUT(beta12 * capd::abs(power(KS::lambda(m, mu_upper), (d12 - capd::interval(2)) / 2)
                    / power(-cmuupper, capd::interval(d1) / 2)));
        COUT(capd::abs(power(KS::lambda(m, mu_upper), (d12 - capd::interval(2)) / 2)
                    / power(-cmuupper, capd::interval(d12) / 2)));
        std::cout << "COND 54 FAILED\n";
    }
    
    if((2 - 6 * gamma_plus * gamma_plus) + beta1 * capd::abs(power(KS::lambda(m, mu_upper), (d1 - capd::interval(2)) / 2)
                    / power(-cmuupper, capd::interval(d1) / 2)) < 0)
    {
        std::cout << "COND 55 SATISFIED\n";
    } else {
        std::cout << "COND 55 FAILED\n";
        COUT(beta1); COUT(power(KS::lambda(m, mu_upper), (d1 - capd::interval(2)) / 2)); // todo fix in paper
        COUT(power(-cmuupper, capd::interval(d1) / 2));
        COUT(beta1 * capd::abs(power(KS::lambda(m, mu_upper), (d1 - capd::interval(2)) / 2)
                    / power(-cmuupper, (d1 - capd::interval(0)) / 2)));
        COUT((2 - 6 * gamma_plus * gamma_plus) + beta1 * capd::abs(power(KS::lambda(m, mu_upper), (d1 - capd::interval(2)) / 2)
                    / power(-cmuupper, (d1 - capd::interval(0)) / 2)));
    }

    cond = true;
    for(int k = m; k < M; ++k) {
        int nk = bounds[k].degree;
        capd::interval betaktimesk = derivbounds[k].coeff;
        if(2 * KS::lambda(k + 1, mu_upper)
           + betaktimesk * Cmuuper < -l){
                continue;;
           }
        
        COUT("COND 57 FAILED");
        cond = false;
    }
    if(cond) {
        std::cout << "COND 57 SATISFIED\n";
    }

    if(
        capd::min(gamma_minus * (1 - gamma_minus * gamma_minus), gamma_plus * (1 - gamma_plus * gamma_plus))
            - alpha1 * 
                capd::abs(power(KS::lambda(m, mu_upper), (n1 - capd::interval(3)) / 2)
                    / power(-cmuupper, (n1 - capd::interval(1)) / 2)) > 0
    ) {
        std::cout << "COND 58 SATISFIED\n";
    } else {
        std::cout << "COND 58 FAILED\n";
    }


    if(
        capd::max(capd::abs(2 - 6 * gamma_plus * gamma_plus), capd::abs(2 - 6 * gamma_minus * gamma_minus)) 
        * KS::lambda(m, mu_upper) + beta1 * power(Cmuuper, d1) < l / 4
        // 10 + beta1 * 
        //      capd::abs(power(KS::lambda(m, mu_upper), (d1 - capd::interval(3)) / 2)
        //             / power(-cmuupper, (d1 - capd::interval(1)) / 2)) < 1 / (4 * KS::lambda(m, mu_upper))
    ) {
        std::cout << "COND 59 SATISFIED\n";
    } else {
        COUT(capd::max(capd::abs(2 - 6 * gamma_plus * gamma_plus), capd::abs(2 - 6 * gamma_minus * gamma_minus)) 
        * KS::lambda(m, mu_upper) + beta1 * power(Cmuuper, d1));
        std::cout << "COND 59 FAILED\n";
        // 10 * KS::lambda(m, mu_upper) + beta1 * power(Cmuuper, d1);
    }


    if(m > 1) { // unstable directions present
        std::cout << "**********\nUNSTABLE DIRECTIONS CONDS\n**********\n";
        cond = true;
        for(int k = 0; k < m - 1; ++k) {
        int nk = bounds[k].degree;
        capd::interval alphakbyktop = bounds[k].coeff;
        if(KS::lambda(k + 1, mu_bif) / power(capd::interval(k + 2), s) // k + 2 because we want unstable dirs to be / 2^s
           - alphakbyktop * power(K, nk) * power(Cmuuper, nk - omega) /* / power(capd::interval(k + 1), p) */ > 0){
                continue;
           }
            cond = false;
            std::cout << "COND 47 FAILED\n";
            COUT(k + 1); COUT(mu_bif);
            COUT(KS::lambda(k + 1, mu_bif) / power(capd::interval(k + 2), s));
            COUT(alphakbyktop * power(K, nk) * power(Cmuuper, nk - omega));
        }
        if(cond){
            std::cout << "COND 47 SATISFIED\n";
        }


        cond = true;
        for(int k = 0; k < m - 1; ++k) {
            int nk = bounds[k].degree;
            capd::interval betaktimesk = derivbounds[k].coeff;
            if(2 * KS::lambda(k + 1, mu_lower) // put 2 * in paper
            - betaktimesk * power(capd::abs(KS::lambda(m, mu_lower)), d1) > l){
                    continue;
            }
            cond = false;
            std::cout << "COND 52 FAILED\n";
            COUT(k);
            COUT(betaktimesk * capd::abs(KS::lambda(m, mu_lower)));
        }
        if(cond) {
            std::cout << "COND 52 SATISFIED\n";
        }


        cond = true;
        for(int k = 0; k < m - 1; ++k) {
            int nk = bounds[k].degree;
            capd::interval betaktimesk = derivbounds[k].coeff;
            if(2 * KS::lambda(k + 1, mu_bif)
            - betaktimesk * Cmuuper > l){
                    continue;;
            }
            COUT(2 * KS::lambda(k + 1, mu_bif)
            - betaktimesk * Cmuuper);
            COUT("COND 56 FAILED");
            cond = false;
        }
        if(cond) {
            std::cout << "COND 56 SATISFIED\n";
        }
    }


}