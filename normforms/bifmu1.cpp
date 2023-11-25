#include "Equations.h"
#include "KS.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "check_bif.h"
using namespace std;



const int M = 4;
const int M2p1 = 2 * M + 1;


const capd::interval nu = capd::interval(99) / 100;
const capd::interval nu_lower = capd::interval(11) / 10;
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

int main(int argc, char const *argv[])
{
    std::vector<int> c1(M2p1, 0), c2(M2p1, 0);
    c1[2] = 1; c2[0] = 1; c2[1] = 1;
    auto change0 = Equation(c1, 1) + Equation(c2, 6 / (KS::lambda(3, nu) - KS::lambda(1, nu) - KS::lambda(2, nu)));
    
    std::fill(c1.begin(), c1.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    c1[0] = 1; c2[0] = 1; c2[1] = 1;
    auto change1 = Equation(c1, 1) + Equation(c2, 2 / (KS::lambda(2, nu)));
    
    std::fill(c1.begin(), c1.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    c1[1] = 1; c2[0] = 2;
    auto change2 = Equation(c1, 1) + Equation(c2, 2 / (KS::lambda(2, nu) - 2 * KS::lambda(1, nu)));


    capd::interval Cmuuper = sqrt(-KS::lambda(1, nu) * KS::lambda(2, nu) / 4);
    const auto c = Constant(Cmuuper);
    std::vector<Bound> v;
    v.push_back(Bound(zeta, 1, &c));

    for(int i = 2; i < M2p1; ++i){
        v.push_back(Bound(power(capd::interval(1) / i, s), omega, &c));
    }
    auto set = Set(v, Bound(1, omega, &c));

    auto system = KS::fullKS(M, s, nu);
    // auto system = KS::gal_proj(2 * M, nu);

    auto new_system = system;
    // // // auto new_system = system.apply(2, change0, true);
    new_system = new_system.apply(0, change1, true).apply(1, change2, true);
    std::vector<int> monomial_v(2 * M + 1, 0);


    std::fill(c1.begin(), c1.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    c1[0] = 1; c2[0] = 2; c2[2] = 1;
    Monomial m(c2);
    new_system = new_system.apply(0, 
        Equation(c1, 1) 
        + Equation(c2, -new_system[0][m] * remove_coeff(0, c2)), true);

    
    std::fill(c1.begin(), c1.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    c1[2] = 1; c2[0] = 3;
    m = Monomial(c2);
    
    new_system = new_system.apply(2, 
        Equation(c1, 1) 
        + Equation(c2, -new_system[2][m] / (KS::lambda(3, nu) - 3 * KS::lambda(1, nu))), true);


    //  on the first mode removing linear part and x_0^3 which we treat theoretically
    new_system = new_system.nonlinear_part();
    std::fill(c1.begin(), c1.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    c1[0] = 3;
    new_system[0].remove(Monomial(c1));


    // COUT(new_system.get_series()[0].get_value());
    // new_system[0].print_lowest_order_terms(set);
    // new_system[0].hsum().print_lowest_order_terms(set);
    // new_system.get_series()[0].get_value().print_lowest_order_terms(set);
    // new_system.get_series()[0].get_hsum().print_lowest_order_terms(set);
    // COUT(new_system);
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

    // for(auto &series: new_system.get_all_series()){
    //     series.print(set);
    //     break;
    // }
    

    COUT(new_system.check_geom_args(capd::interval(1) / 2, set));
    int m_ = 1; int M_ = M2p1;
    capd::interval mu_lower = nu_lower; capd::interval mu_upper = nu;
    capd::interval K = sqrt(capd::interval(2)); capd::interval cmuupper = 4 / KS::lambda(2, nu);
    capd::interval cmubif = 4 / KS::lambda(2, 1);
    /*capd::interval Cmuuper; */ int s = 6; int p = 4;
    capd::interval l = 6; /* int omega */; capd::interval gamma_minus = 1 / capd::interval(2);
    capd::interval gamma_plus = 1 / sqrt(capd::interval(2));

    gamma_minus = 1 / capd::interval(20);
    gamma_plus = 1 / 1.03;
    
    std::vector<Bound> v2;

    v2.push_back(Bound(gamma_minus, 1, &c));

    for(int i = 2; i < M2p1; ++i){
        v2.push_back(Bound(power(capd::interval(1) / i, s), omega, &c));
    }
    auto set2 = Set(v2, Bound(1, omega, &c));

    auto db2 = new_system.deriv_bounds(set2);
    check_bif(r, deriv_bounds, m_, M_, mu_lower, 1, mu_upper, K, cmuupper,
              Cmuuper, s, p, l, omega, gamma_minus, gamma_plus,
              db2, zeta, cmubif);
    return 0;
}
