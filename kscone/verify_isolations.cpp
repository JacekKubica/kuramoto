#include <iostream>
#include "capd/capdlib.h"
#include <random>
 
using namespace capd;
using namespace std;
using capd::autodiff::Node;

const int max_der_order = 4;

#include "PolynomialBound.h"
#include "../../../kscone/LogNorm.h"
#include "new_linear_utils.h"

template <int M>
using PolyBd<M> = capd::pdes::PolynomialBound<capd::interval, int, M>;

int main(int argc, char const *argv[])
{
    IVector target{
        interval(0.17392, 0.17392) + interval(-0.00717178, 0.00717178),
        interval(-0.00510524, -0.00510524) + interval(-0.000486048, 0.000486048),
        interval(7.13344e-05, 7.13344e-05) + interval(-6.40317e-05, 6.40317e-05),
        interval(-1.10565e-07, -1.10565e-07) + interval(-1.25043e-06, 1.25043e-06),
        interval(1.83035e-10, 1.83035e-10) + interval(-1.82125e-08, 1.82125e-08),
        interval(9.4895e-12, 9.4895e-12) + interval(-6.29512e-10, 6.29512e-10),
        interval(1.18624e-13, 1.18624e-13) + interval(-1.18034e-11, 1.18034e-11),
        interval(-9.93583e-15, -9.93583e-15) + interval(-2.47285e-13, 2.47285e-13),
        interval(-6.60653e-16, -6.60653e-16) + interval(-5.50441e-15, 5.50441e-15),
    }
    PolyBd<9> target_pb(interval(-0.00634106, 0.00634106), 12, target.begin());

    calc_lognorm<m, M>(target_pb.C(), target_pb.exponent(), niu, N, A, Ainv, V);


    return 0;
}
