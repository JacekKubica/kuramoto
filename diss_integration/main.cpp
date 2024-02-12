#include "MyDissipativeSolver.h"
#include "KSExtendedMultimapContainer.h"
#include "PolynomialBound.h"
using namespace capd;


#include <algorithm>
#include <vector>

// f is ExtentedMultiMapContainer
capd::IVector findIsolation(auto f, const capd::IVector &initial, int firstDiss=0) {
    auto B = initial;

    std::vector<bool> validated(initial.dimension());
    while(!std::all_of(validated.begin(), validated.end(), [](auto b) { return b; })) {
        // COUT(B);
        // COUT(f.getLinear());
        for(auto b: validated) std::cout << int(b);
        std::cout << "\n";
        auto temp = B;
        for(int i = 0; i < initial.dimension(); ++i) {
            // COUT(i);
            auto dir = i < firstDiss ? 1 : -1;
            temp[i] = B[i].rightBound();
            // COUT(temp);
            // COUT((f.getLinear() * temp +  f.evalNonlinearPlusPerturb(temp))[i]);
            validated[i] = (f.getLinear() * temp +  f.evalNonlinearPlusPerturb(temp))[i] * dir > 0;
            
            temp[i] = B[i].leftBound();
            // COUT(temp);
            // COUT((f.getLinear() * temp +  f.evalNonlinearPlusPerturb(temp))[i]);
            validated[i] = validated[i] && 
                           (f.getLinear() * temp +  f.evalNonlinearPlusPerturb(temp))[i] * dir < 0;

            temp[i] = B[i];
        }
        auto N = f.evalNonlinearPlusPerturb(B);
        auto lambdas = f.diagonalCoefficients();
        for (int i = firstDiss; i < initial.dimension(); ++i) {
            if(!validated[i]){
                B[i].setRightBound(2 * (N[i].rightBound() / (-lambdas[i])).rightBound());
                B[i].setLeftBound(2 * (N[i].leftBound() / (-lambdas[i])).leftBound());
            }
        }
    }
    return B;
}

void node_no_params(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        out[i] = k * k * (1 - param[0] * k * k) * in[i];
        for(int j = 0; j < i; ++j){
            out[i] -= k * in[j] * in[i - 1 - j];
        }
        for(int j = 0; j + k < dimIn; ++j){
            out[i] += 2 * k * in[j] * in[j + k];
        }
    }
}

void test(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    out[0] = in[0] - in[0] * in[0];
}

int main() {
    constexpr int M = 16, m = 8;;
    // const capd::interval niu = capd::interval(99) / 100;
    const capd::interval niu = capd::interval(75) / 100;
    
    

    capd::IVector N(M);
    capd::IVector N_(m);
    // N[0] = N_[0] = capd::interval(-0.0482822, 0.0482822);
    N[0] = N_[0] = 0.075 * capd::interval(-1, 1);
    capd::pdes::PolynomialBound<capd::interval, unsigned, M> pb(M, 2.31307e+08 * capd::interval(-1, 1), 20, N.begin());
    
    // TODO what happens when pb changes?
    KSExtendedMultiMapContainer<m, M> f(niu, pb, true);
    COUT(findIsolation(f, N_, 1));

    N_ = {capd::interval(-7.500000e-002, 7.500000e-002),    capd::interval(-1.406709e-003, 4.587986e-007),    capd::interval(-1.223300e-005, 1.223300e-005),    capd::interval(-8.667864e-008, 4.170524e-008),    capd::interval(-5.342963e-010, 5.342963e-010),    capd::interval(-3.036260e-012, 1.265905e-012),    capd::interval(-1.630102e-014, 1.630102e-014),    capd::interval(-8.397013e-017, 5.074158e-017)};
    for(int i = 0; i < m; ++i) {
        COUT(i);
        auto y = N_;
        y[i] = N_[i].leftBound();
        COUT((f.getLinear()*y + f.getNonlinear()(y))[i]);
        COUT((f.getLinear()*y)[i]);
        COUT(-f.getNonlinear()(y)[i] / f.diagonalCoefficients()[i]);
        COUT(f(y)[i]);
        COUT("=");
        y[i] = N_[i].rightBound();
        COUT((f.getLinear()*y + f.getNonlinear()(y))[i]);
        COUT((f.getLinear()*y)[i]);
        COUT(-f.getNonlinear()(y)[i] / f.diagonalCoefficients()[i]);
        COUT(f(y)[i]);
    }

    // MyDissipativeSolver<KSExtendedMultiMapContainer<m, M>> mds(f);
    // auto x = mds.enclosure(1e-10, N_);
    // // IVector x(m);
    // // x[0] = N_[0];
    // // x[1] = 7e-4 * capd::interval(-1, 1);
    // // x[2] = 1.2e-5 * capd::interval(-1, 1);
    // // x[1] = 6e-8 * capd::interval(-1, 1);
    // x[0] = capd::interval(-7.500000e-002,7.500000e-002);    x[1] = capd::interval(-1.406709e-003,4.587986e-007);    x[2] = capd::interval(-1.223300e-005,1.223300e-005);    x[3] = capd::interval(-8.667864e-008,4.170524e-008);    x[4] = capd::interval(-5.342963e-010,5.342963e-010);    x[5] = capd::interval(-3.036260e-012,1.265905e-012);    x[6] = capd::interval(-1.630102e-014,1.630102e-014);    x[7] = capd::interval(-8.397013e-017,5.074158e-017);

    // auto f_ = capd::IMap(node_no_params, m, m, 1); f_.setParameter(0, niu);
    
    // for(int i = 0; i < m; ++i) {
    //     COUT(i);
    //     auto y = x;
    //     y[i] = x[i].rightBound();
    //     COUT(f_(y)[i]);
    //     COUT(f(y)[i]);
    //     y[i] = x[i].leftBound();
    //     COUT(f_(y)[i]);
    //     COUT(f(y)[i]);
        
    // }
    // COUT(x);



 
}




// simple test to understand a bit about jets
// auto t = capd::IMap(test, 1, 1, 0, 100);
// auto t = capd::IMap("var:x;fun:x*(1+x);", 5);
// capd::IJet tz(1, 1, 4);
// capd::IVector tx {capd::interval(1 - 0.1, 1 + 0.1)};
// t(tx, tz);
// auto p = tz.first(1);
// capd::IVector tlin(1);
// tlin[0] = 0;
// tlin[0] = tz(p)[0];
// tz(p)[0] = 0;
// COUT(t(tx));
// COUT(tz({capd::interval(- 0.1, 0.1)}));
// COUT(tlin * capd::interval(- 0.1, 0.1) + tz({capd::interval(- 0.1, 0.1)}));    


// stable isolation test
// const capd::interval niu(1.1);
// capd::IVector R(M);
// capd::IVector R_(m);
// // N[0] = N_[0] = capd::interval(-0.0482822, 0.0482822);
// R[0] = R_[0] = 0.01 * capd::interval(-1, 1);
// R[1] = R_[1] = 0;
// R[2] = R_[2] = 0;
// capd::pdes::PolynomialBound<capd::interval, unsigned, M> pb(M, 9077.32 * capd::interval(-1, 1), 10, R.begin());
// KSExtendedMultiMapContainer<m, M> f(niu, pb, true);
// // capd::IMatrix A = capd::IMatrix::Identity(m);
// // auto emm2 = changeVariables(A, R_, f, 3);
// // R_ = A * (R_ - centralized(R_));
// MyDissipativeSolver<decltype(f)> mds(f);
// auto x = mds.enclosure(1e-10, R_);


// COUT(x);    
// for(int i = 0; i < m; ++i) {
//     COUT(i);
//     auto y = x;
//     y[i] = x[i].rightBound();
//     COUT(f(y)[i]);
//     y[i] = x[i].leftBound();
//     COUT(f(y)[i]);
// }


// pointless to try to do it by using enclosure
// capd::IVector R(M);
// capd::IVector R_(m);
// // N[0] = N_[0] = capd::interval(-0.0482822, 0.0482822);
// R[0] = R_[0] = 0.712361 + 0.01 * capd::interval(-1, 1);
// R[1] = R_[1] = -0.123239;
// R[2] = R_[2] = 0.0101786;
// capd::pdes::PolynomialBound<capd::interval, unsigned, M> pb(M, 9077.32 * capd::interval(-1, 1), 10, R.begin());
// KSExtendedMultiMapContainer<m, M> f(niu, pb, true);
// capd::IMatrix A = capd::IMatrix::Identity(m);
// auto emm2 = changeVariables(A, R_, f, 3);
// auto R__ = A * (R_ - centralized(R_));
// // COUT(emm2(R_));
// // COUT(emm2.getLinear() * R_ + emm2.evalNonlinearPlusPerturb(R_));
// MyDissipativeSolver<decltype(emm2)> mds(emm2);
// auto x = mds.enclosure(1e-10, R__);

// COUT(x);    
// for(int i = 0; i < m; ++i) {
//     COUT(i);
//     auto y = x;
//     y[i] = x[i].rightBound();
//     auto N = emm2.evalNonlinearPlusPerturb(y);
//     auto L = emm2.getLinear() * y;
//     auto lambdas = emm2.diagonalCoefficients();
//     COUT(N[i]);
//     COUT(L[i]);
//     auto t = L + N;
//     COUT(t[i]);
//     COUT(-N[i] / lambdas[i] - y[i]);
//     y[i] = x[i].leftBound();
//     N = emm2.evalNonlinearPlusPerturb(y);
//     L = emm2.getLinear() * y;
//     t = L + N;

// }

// COUT(capd::matrixAlgorithms::inverseMatrix(A) * x + centralized(R_));   