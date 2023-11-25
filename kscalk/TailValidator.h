#pragma once
#include "capd/capdlib.h"
#include "PolynomialBound.h"
#include "DiffInclusionDiss.h"


template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm> 
class TailValidator{
public:
    TailValidator(const typename Scalar::BoundType h, const VectorType &x, const capd::pdes::PolynomialBound<Scalar, Exponent, M> &T0, const ParamType &p);
    bool validate_far_tail(bool update=false);
    bool validate_tail(bool gen_new);
    capd::pdes::PolynomialBound<Scalar, Exponent, M>& get_T() { return T; }
    const capd::pdes::PolynomialBound<Scalar, Exponent, M>& get_T0() const { return T0; }
    capd::pdes::PolynomialBound<Scalar, Exponent, M> get_Th() const;
    typename Scalar::BoundType get_h() const { return h; }
    const DissipativeOperator& get_diss_op() const { return diss_op; }
    VectorType& get_W2() { return W2; }
    capd::diffIncl::DiffInclusionDiss<capd::IMultiMap, DissipativeOperator, m>& get_diff_inclusion_cw_solver(){
        return diff_inclusion_cw_solver; } 
    void reset(const capd::pdes::PolynomialBound<Scalar, Exponent, M> &new_T0, const VectorType &x, const VectorType &W2);
    void decpower(int dec) { T.C() /= take_power(M + 1, dec); T.exponent() -= dec; set_N_b_and_g(); }
    void print() {
        std::cout << "T0: " << T0 << '\n';
        std::cout << "T: " << T << '\n';
    }

    void print_centered_T(){
        for(int i = m; i < M; ++i){
            std::cout << "x" << i << " = " << T[i].mid() << " +/- " << T[i].rightBound() - T[i].mid() << '\n';
        }
    }

    void print_centered(){
        for(int i = m; i < M; ++i){
            std::cout << "x" << i << " = " << T0[i].mid() << " +/- " << T0[i].rightBound() - T0[i].mid() << '\n';
        }
        std::cout << "C = " << T0.C() << '\n';
        std::cout << "s = " << T0.exponent() << '\n';
    }

    void fill_values(VectorType &v) const {
        for(int i = m; i < M; ++i){
            v[i] = T0[i];
        }
    }

    void test(){
        std::cout << "TEST: " << '\n';
        for(int i = m; i < M; ++i){
            std::cout << T0[i] << '\n';
            std::cout << T[i] << '\n';
            std::cout << '\n';
        }
    }

    void what_failed_in_validation();

// private:
    capd::pdes::PolynomialBound<Scalar, Exponent, M> T0, T, N, b;
    typename Scalar::BoundType h;
    VectorType W2;
    bool validated, far_tail_validated;
    std::vector<bool> kvalidated;
    DissipativeOperator diss_op;
    capd::diffIncl::DiffInclusionDiss<capd::IMultiMap, DissipativeOperator, m> diff_inclusion_cw_solver;
    const double DG = 0.1, D2 = 1.01; 

    Scalar eval_g(int i) const;
    static void inflate(Scalar &x, typename Scalar::BoundType c);
    static bool check_far_tail_inclusion(const capd::pdes::PolynomialBound<Scalar, Exponent, M>&, const capd::pdes::PolynomialBound<Scalar, Exponent, M>&);

    void set_N_b_and_g();
    void update_tail(bool only_far_tail=false);

    static Scalar take_power(int a, const Exponent &b){
        return capd::pdes::PolynomialBound<Scalar, Exponent, M>::takePower(a, b);
    }
};

#include "TailValidator.cpp"
