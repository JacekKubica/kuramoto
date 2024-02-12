#include <iostream>
#include "capd/capdlib.h"
 
using namespace capd;
using namespace std;
using capd::autodiff::Node;

const int max_der_order = 4;

#include "../kscalk/PolynomialBound.h"
#include "LogNorm.h"
#include "new_linear_utils.h"

#include "../kscalk/DissipativeIntegrator.h"
#include "../kscalk/KSDissipativeOperator.h"
#include "ConeConditionVerifier.h"

template <int m, int M>
const capd::pdes::PolynomialBound<capd::interval, int, M> 
    verify_trapping_region(interval niu, const IVector &V,
                           const capd::pdes::PolynomialBound<capd::interval, int, M> &R_pb) {

    auto R = R_pb.mainModes();
    auto sys = System<m, M>(niu);
    IMap fm = sys.getfm();
    IMap fM = sys.getfM();
    IJet fmz(M, M, max_der_order - 1);
    IJet fmz_remainder(M, M, max_der_order);
    auto V_moved = moved_to_zero(V);

    fm(centralized(V), fmz); // calculate the jet
    IMatrix new_linear_part(M, M);
    auto mp = fmz.first(1);
    for(int k = 0; k < m; ++k){
        for(int i = 0; i < m; ++i){
            new_linear_part[i][k] = fmz(mp)[i];
            fmz(mp)[i] = 0;
        }
        for(int i = m; i < M; ++i){
            new_linear_part[i][k] = 0;
        }
        fmz.hasNext(mp);
    }
    fun f(fmz, fM, zero_perturb<m, M>(), centralized(V));

    // change of basis actually return inverse
    auto Ainv = change_of_basis(centralized_double(V), centralized_double(new_linear_part), m); 
    auto A    = capd::matrixAlgorithms::inverseMatrix(Ainv);

    const int first_stable = 0;
    bool isolation_main = check_isolation(
            R, new_linear_part, f, first_stable, A, Ainv, -1, false
          );
    bool isolation_tail = check_far_tail_isolation<m, M>(
        niu, R_pb
    );
    auto lognorm = calc_lognorm<m, M>(R_pb.C(), R_pb.exponent(), niu, R, A, Ainv, V);
    COUT(isolation_main);
    COUT(isolation_tail);
    COUT(lognorm);

    auto temp = Ainv * R + V;
    return capd::pdes::PolynomialBound<capd::interval, int, M>(M,
                                                               R_pb.C(), R_pb.exponent(),
                                                               temp.begin());
}

template <int m, int M>
const capd::pdes::PolynomialBound<capd::interval, int, M> 
    verify_isolating_segment(interval niu, const IVector &V,
                           const capd::pdes::PolynomialBound<capd::interval, int, M> &N_pb) {

    auto N = N_pb.mainModes();
    auto sys = System<m, M>(niu);
    IMap f = sys.get_full();

    // change of basis actually return inverse
    auto Ainv = IMatrix::Identity(M);
    auto A    = IMatrix::Identity(M);
    IMatrix new_linear_part(M, M);
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < M; j++)
        {
            new_linear_part[i][j] = 0;
        }
        
    }

    const int first_stable = 1;
    bool isolation_main = check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv, -1, false, false
          );
    bool isolation_tail = check_far_tail_isolation<m, M>(
        niu, N_pb
    );
    auto conecond = calc_conecond<m, M>(N_pb.C(), N_pb.exponent(), niu, N, A, Ainv, V);

    COUT(isolation_main);
    COUT(isolation_tail);
    COUT(conecond);

    auto temp = Ainv * N + V;
    return capd::pdes::PolynomialBound<capd::interval, int, M>(M,
                                                               N_pb.C(), N_pb.exponent(),
                                                               temp.begin());
}