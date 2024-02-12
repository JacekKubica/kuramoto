#pragma once
#include <random>
#include <exception>
#include "ConeConditionVerifier.h"
#include "../kscalk/KSDissipativeOperator.h"

const int MAX_ITER = 10000;
const double D = 1.01;
const double D2 = .99;

template <int m, int M>
struct System {
    interval niu;

    System(interval niu) : niu(niu) {}

    static void far_modes(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
        for(int i = 0; i < m; ++i) out[i] = Node(0);
        for(int i = m; i < dimOut; ++i){
            int k = i + 1;
            out[i] = k * k * (1 - param[0] * k * k) * in[i];
            for(int j = 0; j < i; ++j){
                out[i] -= k * in[j] * in[i - 1 - j];
            }
            for(int j = 0; j + k < dimIn; ++j){
                out[i] += 2 * k * in[j] * in[j + k];
            }
            // out[i] += param[i + 1];
        }
    }

    static void main_modes(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
        for(int i = 0; i < m; ++i){
            int k = i + 1;
            out[i] = k * k * (1 - param[0] * k * k) * in[i];
            for(int j = 0; j < i; ++j){
                out[i] -= k * in[j] * in[i - 1 - j];
            }
            for(int j = 0; j + k < dimIn; ++j){
                out[i] += 2 * k * in[j] * in[j + k];
            }
            // out[i] += param[i + 1];
        }
        for(int i = m; i < dimOut; ++i) out[i] = Node(0);
    }

    static void full(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
        for(int i = 0; i < dimOut; ++i){
            int k = i + 1;
            out[i] = k * k * (1 - param[0] * k * k) * in[i];
            for(int j = 0; j < i; ++j){
                out[i] -= k * in[j] * in[i - 1 - j];
            }
            for(int j = 0; j + k < dimIn; ++j){
                out[i] += 2 * k * in[j] * in[j + k];
            }
            // out[i] += param[i + 1];
        }
    }

    IMap getfm() {
        auto f = IMap(main_modes, M, M, M, max_der_order);
        f.setParameter(0, niu);
        return f;
    }
    IMap getfM() {
        auto f = IMap(far_modes, M, M, M, max_der_order);
        f.setParameter(0, niu);
        return f;
    }

    IMap get_full() {
        auto f = IMap(full, M, M, M, max_der_order);
        f.setParameter(0, niu);
        return f;
    }
};

IVector moved_to_zero(IVector V) {
    for(auto &x: V) {
        x = x - x.mid();
    }
    return V;
}

IVector centralized(IVector V) {
    for(auto &x: V) {
        x = x.mid();
    }
    return V;
}

DVector centralized_double(IVector V) {
    DVector R(V.dimension());
    for(int i = 0; i < V.dimension(); ++i) {
        R[i] = V[i].mid().leftBound();
    }
    return R;
}

DMatrix centralized_double(const IMatrix &A) {
    DMatrix R(A.numberOfRows(), A.numberOfColumns());
    for(int i = 0; i < A.numberOfRows(); ++i){
        for(int j = 0; j < A.numberOfColumns(); ++j){
            R[i][j] = A[i][j].mid().leftBound();
        }
    }
    return R;
}

IMatrix change_of_basis(const DVector &V, const DMatrix &Derivative, int m){
    int M = V.dimension();
    DVector eigenRealPart(m);
    DVector eigenImPart(m);
    DMatrix vectorRealPart(m, m);
    DMatrix vectorImPart(m, m);

    DMatrix temp(m, m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            temp[i][j] = Derivative[i][j];
        }
    }

    alglib::computeEigenvaluesAndEigenvectors(temp, eigenRealPart, eigenImPart, vectorRealPart, vectorImPart);
    IMatrix A = IMatrix::Identity(M);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            A[i][j] = vectorRealPart[i][j];
        }
    }
    return A;
}

template <class FType, class VType>
bool check_isolation(VType &x, const IMatrix &L, const FType &f,
                     int firstStable, const IMatrix &A, const IMatrix Ainv,
                     int end_variable=-1, bool improve=false, bool log=false) {
    if(log) std::cout << "Checking isolations\n" << x << "\n";
    double D = 1.01;
    bool flag = true;
    for(int i = 0; i < ((end_variable == -1) ? x.dimension() : end_variable); ++i) {
        interval direction = (i < firstStable) ? -1 : 1; // outwards (-1) for unstable, inwards for stable
        
        if(log) {COUT(i);
            std::cout << "Right isolation\n";}
        auto temp = x;
        temp[i] = temp[i].rightBound();
        auto vf = A * L * Ainv * temp + A * f(Ainv * temp) * direction;
        if(log) {COUT(vf[i]);
            COUT((vf[i] < 0));}
        if(!(vf[i] < 0)) 
        {
            if(improve) {
                x[i] = interval(x[i].leftBound(), x[i].rightBound() * D);
                flag = false;
            }
            else return false;
        }
        if(log) std::cout << "Left isolation\n";
        temp = x;
        temp[i] = temp[i].leftBound();
        vf = A * L * Ainv * temp + A * f(Ainv * temp) * direction;
        if(log) {COUT(vf[i]);
            COUT((vf[i] > 0));}
        if(!(vf[i] > 0)) {
            if(improve) {
                x[i] = interval(x[i].leftBound() * D, x[i].rightBound());
                flag = false;
            }
            else return false;
        }
    }
    return flag;
}

struct fun {
    IJet fmz;
    IMap fM, perturb;
    IVector shift;

    fun(const IJet &fmz, const IMap &fM, const IMap &perturb, const IVector &shift) 
        : fmz(fmz), fM(fM), perturb(perturb), shift(shift) {}

    IVector operator()(const IVector &x) const { // x will be A^-1 V
        return fmz(x) + fM(x + shift) + perturb(x + shift);
    }
};


template <int m, int M>
capd::interval takePower(int a, int b){
    return capd::pdes::PolynomialBound<capd::interval, int, M>::takePower(a, b);
}

template <int m, int M>
capd::interval projection_error(int k, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) {
    int s = pb.exponent(); auto C = pb.C();
    capd::interval res = 0;
    for(int n = m - k - 1; n < M - k - 1; ++n){
        res += 2 * (k + 1) * pb[n] * pb[n + k + 1];
    }
    for(int i = M - k - 1; i < M; ++i){
        res += C * 2 * (k + 1) * capd::abs(pb[i]) / takePower<m, M>(k + i + 2, s) * capd::interval(-1, 1);
    }
    res += 2 * (k + 1) * C * C  / (takePower<m, M>(k + M + 2, s) * (s - 1) * takePower<m, M>(M, s - 1)) * capd::interval(-1, 1);

    return res;
}

template <int m, int M>
void node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        out[i] = param[i];
    }
}

template <int m, int M>
IMap get_perturbation_ks(const capd::pdes::PolynomialBound<capd::interval, int, M> &pb){
    IMap perturb(node_perturb<m, M>, M, M, M);
    for(int i = 0; i < M; ++i){
        auto pe = projection_error<m, M>(i, pb);
        perturb.setParameter(i, pe);
    }
    return perturb;
}


template <int m, int M>
IMap zero_perturb(){
    IMap perturb(node_perturb<m, M>, M, M, M);
    for(int i = 0; i < M; ++i){
        perturb.setParameter(i, 0);
    }
    return perturb;
}

template <int m, int M, class FType>
void main_modes_isolation_guess(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, int first_stable=0) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist0(0, m); // random coord to inflate
    std::uniform_int_distribution<std::mt19937::result_type> dist1(1, 2); // random direction to inflate

    for (size_t i = 0; i < MAX_ITER; i++)
    {
        auto temp = N;
        int k = dist0(rng);
        temp[k] = dist1(rng) == 1 ? interval(temp[k].leftBound() * D, temp[k].rightBound()) 
                                  : interval(temp[k].leftBound(), temp[k].rightBound() * D);
        auto b = check_isolation(
            temp, new_linear_part, f, first_stable, A, Ainv, m, false, false
        );
        if(b) N = temp;
    }

    for(int i = 0; i < first_stable; ++i) {
        N[i] /= 1.2;
    }
    for(int i = first_stable; i < m; ++i) {
        if(first_stable > 0) N[i] /= 1.4;
        else N[i] /= 2;
    }
    auto b = check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv, m, false, true
    );
    std::cout << "Initial guess:\n";
    COUT(N);
    if(!b) throw std::runtime_error("bad initial candidate");
}


template <int m, int M, class FType>
void near_tail_isolation_guess(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, const interval &niu, int first_stable=0) {
    
    for(int i = m; i < M; ++i) {
        N[i] = 2 * f(N)[i] / ((i + 1) * (i + 1) * (1 - niu * (i + 1) * (i + 1))) * interval(-1, 1);
    }
    
    auto b = check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv, -1
    );
    std::cout << "Near tail guess:\n";
    COUT(N);
    if(!b) throw std::runtime_error("bad near tail candidate");                   
}

template <int m, int M, class FType>
void first_refinement(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, int first_stable=0) {

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist1(1, 2); // random coord to deflate
    std::uniform_int_distribution<std::mt19937::result_type> dist2(1, M - 1);
    for (size_t i = 0; i < MAX_ITER * 2; i++)
    {
        auto temp = N;
        int k = dist2(rng);
        temp[k] = dist1(rng) == 1 ? interval(temp[k].leftBound() * D2, temp[k].rightBound()) 
                                  : interval(temp[k].leftBound(), temp[k].rightBound() * D2);
        auto b = check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv
        );
        if(b) N = temp;
    }
    std::cout << "First refinement:\n";
    COUT(N);
}


template <int m, int M, class FType>
void make_lognorm_negative(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, capd::interval niu,
                            const IVector &V) {
    interval lognorm;
    IMatrix Afromone = IMatrix::Identity(M + 1);
        for(int i = 1; i <= M; ++i){
            for(int j = 1; j <= M; ++j){
                Afromone[i][j] = A[i - 1][j - 1];
        }
    }

    auto cond = false;
    while(!cond) {
        for(int i = 0; i < M; ++i) {
            N[i] *= D2;
        }
        auto target_set = Ainv * N + centralized(V);
        IVector target_set_from_one(M + 1);
        for(int i = 1; i <= M; ++i){
            target_set_from_one[i] = target_set[i - 1];
        }

        // no far tail, hence C = 0
        lognorm = LogNorm<m, M>(0, 10, niu, target_set_from_one, Afromone).max_lognorm();
        cond = lognorm < 0;
    }
    std::cout << "Lognorm based refinement:\n";
    COUT(N);
    COUT(lognorm);
}

template <int m, int M>
interval calc_lognorm(const interval &C, int s, const interval &niu, const IVector &N,
                      const IMatrix &A, const IMatrix &Ainv, const IVector &V){
    IMatrix Afromone = IMatrix::Identity(M + 1);
        for(int i = 1; i <= M; ++i){
            for(int j = 1; j <= M; ++j){
                Afromone[i][j] = A[i - 1][j - 1];
        }
    }

    auto target_set = Ainv * N + centralized(V);
    IVector target_set_from_one(M + 1);
    for(int i = 1; i <= M; ++i){
        target_set_from_one[i] = target_set[i - 1];
    }
    return LogNorm<m, M>(C, s, niu, target_set_from_one, Afromone).max_lognorm();
}

template <int m, int M>
interval calc_conecond(const interval &C, int s, const interval &niu, const IVector &N,
                      const IMatrix &A, const IMatrix &Ainv, const IVector &V){
    IMatrix Afromone = IMatrix::Identity(M + 1);
        for(int i = 1; i <= M; ++i){
            for(int j = 1; j <= M; ++j){
                Afromone[i][j] = A[i - 1][j - 1];
        }
    }

    auto target_set = Ainv * N + centralized(V);
    IVector target_set_from_one(M + 1);
    for(int i = 1; i <= M; ++i){
        target_set_from_one[i] = target_set[i - 1];
    }
    auto Q = IVector(M + 1);
    for(int i = 2; i < M + 1; ++i) {
        Q[i] = -1;
    }
    return ConeConditionVerifier<m, M>(C, s, niu, 0, target_set, Q, A).calc_eps();
}

template <int m, int M, class FType>
void fix_isolations(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, capd::interval niu,
                            const IVector &V, int first_stable=0) {
    
    while(!check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv, -1, true
          )) {}
    
    std::cout << "Isolation fixed:\n";
    COUT(N);
}

struct fun2 {
    IMap f, perturb;
    IVector operator()(const IVector &V) const {
        return f(V) + perturb(V);
    }
};

template <int m, int M, class FType>
capd::pdes::PolynomialBound<capd::interval, int, M> make_polybd(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, capd::interval niu,
                            const IVector &V, int min_s, int first_stable=0, int s_max=16) {

    int s = min_s;
    auto f_cpy = f;
    
    capd::pdes::PolynomialBound<capd::interval, int, M> pb(M, 0, s, N.begin());
    bool found = false;
    // we look for the biggest s
    while(true){
        if(s > s_max){
            if(!found) throw std::runtime_error("could not fit polynomial bounds");
            else return pb;
        }
        auto pb_cpy = pb;
        capd::interval C = 2 * pb[M - 1] * takePower<m, M>(M + 1, s) * interval(-1, 1);
        pb_cpy.C() = C;
        pb_cpy.exponent() = s;
        auto N_cpy = N;
        f_cpy.perturb = get_perturbation_ks<M, M>(pb_cpy);
        
        const int MAX_ITER = 20;
        for(int i = 0; i < MAX_ITER; ++i){
            check_isolation(
                N_cpy, new_linear_part, f_cpy, first_stable, A, Ainv, -1, true
            );
        }
        for(int i = 0; i < M; ++i) {
            pb_cpy[i] = N_cpy[i];
        }
        pb_cpy.C() = pb[M - 1] * takePower<m, M>(M, s) * interval(-1, 1);
        f_cpy.perturb = get_perturbation_ks<M, M>(pb);

        bool cond = (first_stable > 0) ? (calc_conecond<m, M>(pb.C(), pb.exponent(), niu, N, A, Ainv, V) > 0)
                                       : (calc_lognorm<m, M>(pb.C(), pb.exponent(), niu, N, A, Ainv, V) < 0);
        if(check_isolation(
            N_cpy, new_linear_part, f_cpy, first_stable, A, Ainv, -1, false) && cond) 
        {    
            found = true;
            N = N_cpy;
            pb = pb_cpy;
        }
        s += 1;
    }
}


struct Ret {
    interval lognorm, C, s;
    IMatrix A, Ainv;
    IVector target;
    IVector target_good_coords;

    void print() {
        COUT(C);
        COUT(s);
        for(const auto &x: target){
            std::cout << x.mid()  << " + " << (x.rightBound() - x.mid()) * interval(-1, 1) << "\n"; 
        }
        COUT(target_good_coords);
        cout << "\n{";
        for(const auto &x: target_good_coords) {
            std::cout << "interval(" << x.leftBound() << ", " << x.rightBound() << "),";
        }
        cout << "}\n";
        COUT(lognorm);
    }
};

template <int m, int M, class FType>
void make_eps_positive(IVector &N, const IMatrix &new_linear_part, const FType &f,
                            const IMatrix &A, const IMatrix &Ainv, capd::interval niu,
                            const IVector &V) {
    interval eps;

    auto cond = false;
    while(!cond) {
        for(int i = 0; i < M; ++i) {
            N[i] *= D2;
        }
        auto target_set = Ainv * N + centralized(V);
        IVector target_set_from_one(M + 1);
        for(int i = 1; i <= M; ++i){
            target_set_from_one[i] = target_set[i - 1];
        }

        // no far tail, hence C = 0
        eps = calc_conecond<m, M>(0, 10, niu, N, A, Ainv, V);
        cond = eps > 0;
    }
    std::cout << "Cone conditions based refinement:\n";
    COUT(N);
    COUT(eps);
}

template <int m, int M>
bool check_far_tail_isolation(capd::interval niu, const capd::pdes::PolynomialBound<capd::interval, int, M> &pb) {
    auto zero_pb = pb; for(int i = 0; i < M; ++i) zero_pb[i] = 0; // class interface compatibility
    auto nonlinear = KSDissipativeOperator<M, M>(niu, zero_pb).N(pb.mainModes());
    return (nonlinear.C() / ((M + 1) * (M + 1) * (1 - niu * (M + 1) * (M + 1)))).subsetInterior(pb.C() * interval(-1, 1));
}

template <int m, int M, class FType>
void extend_direction_zero(capd::pdes::PolynomialBound<capd::interval, int, M> &pb, 
                           const IMatrix &new_linear_part, const FType &f,
                           const IMatrix &A, const IMatrix &Ainv, capd::interval niu, const IVector& V,
                           int first_stable=0) {

    auto f_cpy = f;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist1(1, 2); // random coord to inflate
    std::uniform_int_distribution<std::mt19937::result_type> dist2(0, 2 * M - 1);
    bool found = false;
    for (size_t i = 0; i < MAX_ITER * 2; i++)
    {
        auto temp = pb;
        auto N = temp.mainModes();
        int k = dist2(rng);
        double infldefl = (dist1(rng) == 1) ? D : D2;
        k -= M;
        if(k < 0) k = 0;
        if(k != 0)
            temp[k] = dist1(rng) == 1 ? interval(temp[k].leftBound() * infldefl, temp[k].rightBound()) 
                                    : interval(temp[k].leftBound(), temp[k].rightBound() * infldefl);
        if(k == 0)
            interval(temp[k].leftBound() * D, temp[k].rightBound());
        
        for(int i = 0; i < M; ++i) temp[i] = N[i];
        f_cpy.perturb = get_perturbation_ks<M, M>(temp);

        auto b = check_isolation(
            N, new_linear_part, f, first_stable, A, Ainv
        );
        auto b2 = calc_lognorm<m, M>(pb.C(), pb.exponent(), niu, N, A, Ainv, V) < 0;
        if(b && b2) { found = true; pb = temp; }
    }
    if(!found) throw std::runtime_error("conditions do not hold");
    if(!check_far_tail_isolation<M, M>(niu, pb)) throw std::runtime_error("no far tail isolation");

    std::cout << "Final refinement:\n";
    COUT(pb);
}