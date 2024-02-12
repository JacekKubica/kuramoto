#pragma once
#include "capd/capdlib.h"
#include <exception>
#include "PolynomialBound.h"


template <int m, int M>
class KSDissipativeOperator {
public:
// private:
    static void nodeP(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam);
    static void node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam);
public:
    KSDissipativeOperator(capd::interval niu, const capd::pdes::PolynomialBound<capd::interval, int, M> &T) : niu(niu), T(T), pm(nodeP, m, m, m + 1), perturb(node_perturb, m, m, m), mm{pm, perturb} { pm.setParameter(0, niu); }

    capd::interval V() const { return niu - capd::interval(1) / ((M + 1) * (M + 1)); }
    capd::interval lambda(int l) const { int k = l + 1; return k * k * (1 - niu * k * k); }
    capd::interval E(int r, capd::interval::BoundType h) const {
        return exp(h * lambda(M)) * take_power(M + 1, r);
    }
    int p() const { return 4; }
    capd::IVector P(const capd::IVector &x) const;
    capd::pdes::PolynomialBound<capd::interval, int, M> N(const capd::IVector &z) const;
    capd::pdes::PolynomialBound<capd::interval, int, M> guess_far_tail(const capd::IVector &x, const capd::pdes::PolynomialBound<capd::interval, int, M> &T0) const;
    void update_perturbation(const capd::IVector &z);
    capd::IMultiMap& perturbated() { return mm; }
    void verify_step(const capd::interval &h, int sb, int sT0) const {
        if(-4 * h * niu * take_power(M + 1, 4) + 2 * h * take_power(M + 1, 2) + sb - sT0 > 0 ||
           4 * niu * take_power(M + 1, 2) < 1) { throw std::runtime_error("step not verified"); } // TODO valid error here
    }
// private:
    capd::interval niu;
    const capd::pdes::PolynomialBound<capd::interval, int, M> &T;
    capd::IMap pm, perturb;
    capd::IMultiMap mm;

    static capd::interval take_power(int a, int b){
        return capd::pdes::PolynomialBound<capd::interval, int, M>::takePower(a, b);
    }
    capd::interval D1klarge(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval D1ksmall(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval D1(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval D2(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval FS(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval IS(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;
    capd::interval projection_error(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) const;

public:
    void dump(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) {
        for(int i = 0; i < 2 * M + 3; ++i)
            std::cout << "FS" << i + 1 << ": " << FS(i, pb) << '\n';
        for(int i = 0; i < 2 * M; ++i)
            std::cout << "IS" << i + 1 << ": " << IS(i, pb) << '\n';
        std::cout << "D2: " << D2(pb) << '\n';
        std::cout << "D1klarge: " << D1klarge(pb) << '\n';
        for(int i = 0; i < m; ++i)
            std::cout << "projection_error" << i + 1 << ": " << projection_error(i, pb) << '\n';
    }
};

#include "KSDissipativeOperator.cpp"
