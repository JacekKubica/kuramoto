#pragma once
#include "capd/capdlib.h"
#include "PolynomialBound.h"
#include "TailValidator.h"
#include "DiffInclusionDiss.h"


template<int m, int M, class Scalar, class Exponent, class VectorType, class SetType, class ParamType, class DissipativeOperator, class Norm>
class DissipativeIntegrator{
public:
    using BoundType = typename Scalar::BoundType;

    DissipativeIntegrator(BoundType h, const VectorType &x, const capd::pdes::PolynomialBound<Scalar, Exponent, M> &T0, const ParamType &p);
    bool take_step();
    void print() {
        std::cout << "main modes: " << VectorType(x) << '\n';
        std::cout << "tail:" << '\n';
        std::cout << tv.T0 << '\n';
    }

    void print_centered(const VectorType &x_vect){
        for (size_t i = 0; i < m; i++) {
            std::cout << "x" << i << " = " << x_vect[i].mid() << " +/- " << x_vect[i].rightBound() - x_vect[i].mid() << '\n';
        }
    }

    void print_centered(){
        auto x_vect = VectorType(x);
        for (size_t i = 0; i < m; i++) {
            std::cout << "x" << i + 1 << " = " << x_vect[i].mid() << " +/- " << x_vect[i].rightBound() - x_vect[i].mid() << '\n';
        }
        tv.print_centered();
    }

    VectorType current_values() {
        VectorType v(m + M);
        auto x_vect = VectorType(x);
        for (size_t i = 0; i < m; i++) {
            v[i] = x_vect[i];
        }
        tv.fill_values(v);
        return v;
    }
// private:
    SetType x;
    VectorType W2;
    TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm> tv;
    bool incl_enclosure();
    bool enclosure_with_tail();

    static Scalar take_power(int a, const Exponent &b){
        return capd::pdes::PolynomialBound<Scalar, Exponent, M>::takePower(a, b);
    }
};

#include "DissipativeIntegrator.cpp"
