#include "DissipativeIntegrator.h"
#include "KSDissipativeOperator.h"
#include "capd/capdlib.h"
#include "PolynomialBound.h"
#include <iostream>
using namespace capd;


int main(){
    const int m = 8;
    const int M = 16;
    const capd::interval niu = interval(75) / interval(100);

    IVector x(m);
    x[0] = interval(0.075) + interval(-1e-10, 1e-10);
    x[1] = -7.031250e-004 + 7.035838e-004 * interval(-1, 1);
    x[2] = 0.000000e+000 + 1.223300e-005 * interval(-1, 1);
    x[3] = -2.248670e-008 + 6.419194e-008 * interval(-1, 1);
    x[4] = 0.000000e+000 + 5.342963e-010 * interval(-1, 1);
    x[5] = -8.851774e-013 + 2.151082e-012 * interval(-1, 1);
    x[6] = 0.000000e+000 + 1.630102e-014 * interval(-1, 1);
    x[7] = -1.661427e-017 + 6.735586e-017 * interval(-1, 1);


    capd::pdes::PolynomialBound<interval, int, M> T0(M, 1.28151e7, 16);
    for(int i = 0; i < m; ++i) T0[i] = 0;
    T0[8] = 0 + 4.19092e-019 * interval(-1, 1);
    T0[9] = -4.40581e-022 + 1.59964e-021 * interval(-1, 1);
    T0[10] = 9.73596e-024 * interval(-1, 1);
    T0[11] = -6.79843e-027 + 2.76354e-025 * interval(-1, 1);
    T0[12] = 2.8995e-023 * interval(-1, 1);
    T0[13] = -1.61085e-031 + 3.20088e-021 * interval(-1, 1);
    T0[14] = 2.96259e-019 * interval(-1, 1);
    T0[15] = 1.30238e-017 * interval(-1, 1);

    double h = 0.0005;
    DissipativeIntegrator<m, M, interval, int, IVector, InclRect2Set, interval, KSDissipativeOperator<m, M>, IMaxNorm()> di(h, x, T0, niu); //max norm?
    int i = 0;
    for (i = 0; h * i < 12.4; i++) {
        if(!di.take_step()) break;
    }
    std::cout << i << '\n';
    di.print();
    di.print_centered();
}
