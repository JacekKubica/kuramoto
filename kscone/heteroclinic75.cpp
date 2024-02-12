#include "heteroclinic_proof.h"


int main() {
    constexpr int m_trapping = 3, M_trapping = 10;
    const interval niu = capd::interval(75) / 100;
    
    DVector _V {0.712361, -0.123239, 0.0101786};
    IVector V(M_trapping);
    V[0] = _V[0]; 
    V[1] = _V[1]; 
    V[2] = _V[2]; 
    for(int i = m_trapping; i < M_trapping; ++i){
        V[i] = 0;
    }


    std::cout << "verifying target\n";
    IVector R {interval(-0.0528214, 0.0507603),interval(-0.00574555, 0.00498893),interval(-0.00351964, 0.00282908),interval(-0.00104017, 0.000720668),interval(-7.122e-05, 7.04373e-05),interval(-5.60354e-06, 5.60354e-06),interval(-4.45228e-07, 5.17671e-07),interval(-3.83262e-08, 3.75635e-08),interval(-2.93988e-09, 2.88138e-09),interval(-1.92939e-10, 2.02882e-10)};
    auto R_pb = capd::pdes::PolynomialBound<capd::interval, int, M_trapping>(
        M_trapping,
        interval(-202.882, 202.882),
        12,
        R.begin()
    ); // trapping region in "good" coordinates
    auto target_set = verify_trapping_region<m_trapping, M_trapping>(niu, V, R_pb);
    COUT(target_set);


    std::cout << "verifying origin isolation\n";
    constexpr int m = 8, M = 16;
    IVector N {interval(-0.273971, 0.260674),interval(-0.0402309, 0.0410396),interval(-0.0310602, 0.028675),interval(-0.00154983, 0.00148876),interval(-0.000139405, 0.000145124),interval(-2.92129e-05, 2.89208e-05),interval(-1.92877e-06, 1.87148e-06),interval(-1.46937e-07, 1.59239e-07),interval(-1.78103e-08, 1.76322e-08),interval(-1.30469e-09, 1.21606e-09),interval(-1.01927e-10, 1.01927e-10),interval(-9.3813e-12, 9.3813e-12),interval(-7.17293e-13, 6.48707e-13),interval(-5.58679e-14, 5.20726e-14),interval(-4.45646e-15, 4.32409e-15),interval(-3.49577e-16, 3.53108e-16)};
    auto N_pb = capd::pdes::PolynomialBound<capd::interval, int, M>(
        M, interval(-4.26881e+08, 4.26881e+08), 20, N.begin());
    verify_isolating_segment<m, M>(niu, IVector(M), N_pb);

    
    auto pb = N_pb;
    pb[0] = pb[0].rightBound(); // we integrate from the right boundary
    capd::pdes::PolynomialBound<interval, int, M> T0(M, pb.C(), pb.exponent());
    for(int i = 0; i < m; ++i) T0[i] = 0; for(int i = m; i < M; ++i) T0[i] = pb[i];
    IVector x(m); for(int i = 0; i < m; ++i) x[i] = pb[i];
    x[0] = x[0].rightBound();
    double h = 0.0005;
    DissipativeIntegrator<m, M, interval, int, IVector, InclRect2Set, interval, KSDissipativeOperator<m, M>, IMaxNorm()> di(h, x, T0, niu); //max norm?
    int i = 0;
    try{
        for (i = 0; h * i < 14; i++) {
            if(!di.take_step()) throw std::runtime_error("could not make step");
            // if(i % 1000 == 0){
            //     std::cout << i << '\n';
            //     di.print();
            //     di.print_centered();
            // } 
        }
    }
    catch(...) {
        COUT("ERR");
    }
    di.print_centered();
    std::cout << i << '\n';
    auto calculated_t = di.tv.T0;
    di.tv.print();   
    for(int j = 0; j < m; ++j) {
        COUT(j);
        COUT(subsetInterior(IVector(di.x)[j], target_set[j]));
    }
    for(int j = m; j < M; ++j) {
        COUT(j);
        COUT(subsetInterior(calculated_t[j], target_set[j]));
    }
    auto exp_dif = calculated_t.exponent() - target_set.exponent();
    COUT((exp_dif >= 0 && calculated_t.C().subsetInterior(takePower<m, M>(M, exp_dif) * target_set.C())));   
}