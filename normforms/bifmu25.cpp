#include "Equations.h"
#include "KS.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "check_bif.h"
using namespace std;



const int M = 6;
const int M2p1 = 2 * M + 1;


const capd::interval nu = capd::interval(25) / 100 - capd::interval(2) / 10000;
const capd::interval nu_bif = capd::interval(25) / 100;
const capd::interval nu_lower = capd::interval(25) / 100 + capd::interval(1) / 100;
const capd::interval zeta = 1.2;
const int omega = 3;
const int s = 6;


capd::interval remove_coeff(int i, const std::vector<int> &v) {
    auto ret = KS::lambda(i + 1, nu);
    for(int i = 0; i < v.size(); ++i){
        ret -= v[i] * KS::lambda(i + 1, nu);
    }
    return 1 / ret;
}

Equations remove(int variable, const std::vector<int> &term, const Equations &system) {
    std::vector<int> l(M2p1, 0), c(M2p1, 0);
    std::copy(term.begin(), term.end(), c.begin());
    Monomial m(c);
    l[variable] = 1;
    auto change = Equation(l, 1) + Equation(c, -system[variable][m] * remove_coeff(variable, c));
    return system.apply(variable, change, true);
}

capd::interval coefficient_next_to_b2pow3(const capd::interval &mu) {
    return 16 / KS::lambda(4, mu);
}

int main(int argc, char const *argv[])
{
    capd::interval cmuupper = coefficient_next_to_b2pow3(nu);
    capd::interval cmubif = coefficient_next_to_b2pow3(nu_bif);
    capd::interval Cmuuper = sqrt(KS::lambda(2, nu) / (-cmuupper));// sqrt(-KS::lambda(1, nu) * KS::lambda(2, nu) / 4);
    const auto c = Constant(Cmuuper);
    std::vector<Bound> v;
    v.push_back(Bound(power(capd::interval(1) / 2, s), omega, &c));
    v.push_back(Bound(zeta, 1, &c));
    for(int i = 3; i < M2p1; ++i){
        v.push_back(Bound(power(capd::interval(1) / i, s), omega, &c));
    }
    auto set = Set(v, Bound(1, omega, &c));

    auto system = KS::fullKS(M, s, nu);
    // auto system = KS::gal_proj(2 * M, nu);
    auto new_system = system;
    new_system = remove(1, {0, 1, 0, 1}, new_system);
    new_system = remove(3, {0, 2}, new_system);
    new_system = remove(1, {0, 2, 0, 0, 0, 1}, new_system);
    new_system = remove(5, {0, 3}, new_system);
    // new_system = remove(1, {0, 4}, new_system);


    new_system = new_system.nonlinear_part();
    
    std::vector<int> c1(M2p1, 0);
    c1[1] = 3;
    new_system[1].remove(Monomial(c1));

    auto r = new_system.vals(set);
    auto hsums = new_system.hsums(set);
    COUT(new_system);
    auto vsums = new_system.vsums(set);
    
    auto deriv_bounds = new_system.deriv_bounds(set);
    for (int i = 0; i < new_system.get_no_of_variables(); i++)
    {
        COUT(r[i]);   
    }
    for(int i = 0; i < new_system.get_no_of_variables(); ++i){
        cout << "i = " << i << ", ";
        COUT(hsums[i]);
    }
    for(int i = 0; i < new_system.get_no_of_variables(); ++i){
        cout << "i = " << i << ", ";
        COUT(vsums[i]);
    }
    

    COUT(new_system.check_geom_args(capd::interval(1) / 2, set));
    int m_ = 2; int M_ = M2p1;
    capd::interval mu_lower = nu_lower; capd::interval mu_upper = nu;
    capd::interval K = (1 + 8 * mu_upper) / (-1 + 16 * mu_upper);
    COUT(K)
    int s = 6; int p = 4;
    capd::interval l = 0.066; /* int omega */; capd::interval gamma_minus = 1 / capd::interval(2);
    capd::interval gamma_plus = 1 / sqrt(capd::interval(2));
    gamma_minus = 1 / capd::interval(15);
    gamma_plus = 1 / 1.1;
    
    std::vector<Bound> v2;

    v2.push_back(Bound(power(capd::interval(1) / 2, s), omega, &c));
    v2.push_back(Bound(gamma_minus, 1, &c));
    for(int i = 3; i < M2p1; ++i){
        v2.push_back(Bound(power(capd::interval(1) / i, s), omega, &c));
    }
    auto set2 = Set(v2, Bound(1, omega, &c));

    auto db2 = new_system.deriv_bounds(set2);
    check_bif(r, deriv_bounds, m_, M_, mu_lower, nu_bif, mu_upper, K, cmuupper,
              Cmuuper, s, p, l, omega, gamma_minus, gamma_plus,
              db2, zeta, cmubif);
    return 0;
}
