
#include "capd/capdlib.h"
#include <iostream>
#include "ConeConditionVerifier.h"
#include "LogNorm.h"
#include "../kscalk/KSDissipativeOperator.h"
#include <exception>
using capd::autodiff::Node;

using namespace capd;



const double niu = 0.75;

void _ks_lin(Node /*t*/, Node in[], int dimIn, Node out[], int /*dimOut*/, Node params[], int noParams)
{
  for(int i = 0; i < dimIn; ++i) {
    out[i] = (i + 1) * (i + 1) * (1 - niu * (i + 1) * (i + 1)) * in[i];
  }
}


template <int m, int M>
double dF(int i, int j, const DVector &a, const capd::interval &niu){
    if(i > j) return -2 * i * a[i - j] + 2 * i * a[i + j];
    else if(i == j) return i * i * (1 - niu.rightBound() * i * i) + 2 * i * a[2 * i];
    return 2 * i * (a[j - i] + a[j + i]);
}

template <int m, int M>
IMatrix change_of_basis(const DVector &V, const capd::interval &niu, IMatrix &new_linear_part){
    DMatrix temp(m, m);
    DMatrix temp_(m + 1, m + 1);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            temp[i][j] = dF<m, M>(i + 1, j + 1, V, niu);
            temp_[i + 1][j + 1] = temp[i][j];
        }
    }

    DVector eigenRealPart(m);
    DVector eigenImPart(m);
    DMatrix vectorRealPart(m, m);
    DMatrix vectorImPart(m, m);

    alglib::computeEigenvaluesAndEigenvectors(temp, eigenRealPart, eigenImPart, vectorRealPart, vectorImPart);
    IMatrix A(m + 1, m + 1);
    A[0][0] = 1;
    for (size_t i = 1; i <= m; i++) {
        for (size_t j = 1; j <= m; j++) {
            A[i][j] = vectorRealPart[i-1][j-1];
        }
    }

    auto Ainv = capd::matrixAlgorithms::inverseMatrix(A);
    auto bigger_new_linear_part = Ainv * IMatrix(temp_) * A;
    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < m; ++j) {
            new_linear_part[i][j] = bigger_new_linear_part[i + 1][j + 1];
        }     
    } 
    COUT(A * IMatrix(temp_) * Ainv);
    COUT(Ainv * IMatrix(temp_) * A);
    return A;
}

template <class mapType>
void check_isolation(const IVector &x, const mapType &f,
                     int firstStable, const IMatrix &A, const IMatrix &Ainv,
                     const IVector &shift, const IMatrix linear_part) {
    // Indexing from 0

    std::cout << "Checking isolations\n";
    // auto Ainv = capd::matrixAlgorithms::inverseMatrix(A);
    for(int i = 0; i < x.dimension(); ++i) {
        COUT(i);

        interval coeff = i >= firstStable ? 1 : -1; // whether direction is stable / unstable

        std::cout << "Right isolation\n";
        auto temp = x;
        temp[i] = x[i].rightBound();
        auto linv = linear_part * temp;
        auto vf = A * f(Ainv * temp + shift);
        COUT(vf); COUT(linv);
        COUT(vf[i] - linv[i]);
  
        std::cout << "Left isolation\n";
        temp = x;
        temp[i] = x[i].leftBound();
        linv = linear_part * temp;
        vf = A * f(Ainv * temp + shift);
        COUT(vf); COUT(linv);
        COUT(vf[i] - linv[i]);
    }
}

// approx_fixed_point indexed from 1
template <int m, int M>
void find_isolations(const DVector &approx_fixed_point, const capd::interval &niu, int s) {
    IMatrix new_linear_part(m, m);
    auto A_ = change_of_basis<m, M>(approx_fixed_point, niu, new_linear_part);
    COUT(new_linear_part);
    IMatrix A(m, m);
    for(int i = 0; i < m; ++i){
        for (int j = 0; j < m; ++j)
        {
            A[i][j] = A_[i + 1][j + 1];
        }
    }
    auto Ainv = capd::matrixAlgorithms::inverseMatrix(A);
    
    capd::pdes::PolynomialBound<interval, int, M> T0(M + 1, 0, s);
    // tail is 0 at the beginning
    KSDissipativeOperator<m, M> gal_proj_ks({niu, T0});
    
    double delta = 6e-2;
    auto sym_interval = [](double x) { return interval(-x, x);};
    IVector isolation_candidate = {sym_interval(delta), sym_interval(delta / 1000), sym_interval(delta / 100)};
    IVector shift(m); for(int i = 1; i <= m; ++i) {shift[i - 1] = approx_fixed_point[i];}
    // auto map = IMap(_ks_lin, m, m, 0);
    auto map = gal_proj_ks.perturbated();
    check_isolation(isolation_candidate, map,
                    0, Ainv, A, shift, new_linear_part);
}

int main(){
    find_isolations<3, 10>({0, 0.712361, -0.123239, 0.0101786}, 0.75, 12);
}
