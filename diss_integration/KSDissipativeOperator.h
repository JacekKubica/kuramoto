#pragma once
#include "capd/capdlib.h"
#include <exception>
#include "PolynomialBound.h"


template <int m, int M>
class KSDissipativeOperator {
public:
    KSDissipativeOperator(capd::interval niu) : niu(niu) { }

    capd::interval V() const { return niu - capd::interval(1) / ((M + 1) * (M + 1)); }
    capd::interval lambda(int l) const { int k = l + 1; return k * k * (1 - niu * k * k); }
    capd::interval E(int r, capd::interval::BoundType h) const {
        return exp(h * lambda(M)) * take_power(M + 1, r);
    }
    int p() const { return 4; }
    capd::IVector P(const capd::IVector &x) const;
    capd::pdes::PolynomialBound<capd::interval, unsigned, M> guess_far_tail(const capd::IVector &x, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &T0) const;
    void verify_step(const capd::interval &h, int sb, int sT0) const {
        if(-4 * h * niu * take_power(M + 1, 4) + 2 * h * take_power(M + 1, 2) + sb - sT0 > 0 ||
           4 * niu * take_power(M + 1, 2) < 1) { throw std::runtime_error("step not verified"); }
    }

    
    static capd::interval take_power(int a, int b){
        return capd::pdes::PolynomialBound<capd::interval, unsigned, M>::takePower(a, b);
    }
    capd::interval D1klarge(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval D1ksmall(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval D1(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval D2(const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval FS_minus_galerkin(int k, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval FS(int k, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval IS_minus_galerkin(int k, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval IS(int k, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;
    capd::interval projection_error(int k, const capd::pdes::PolynomialBound<capd::interval, unsigned, M> &pb) const;

private:
    capd::interval niu;

};

#include "KSDissipativeOperator.cpp"
