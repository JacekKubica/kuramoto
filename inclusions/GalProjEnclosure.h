#pragma once
#include "capd/capdlib.h"
#include <stdexcept>


class DissipativeEnclosure {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;

private:
    constexpr static double DISS_CONST = -0.01;
    constexpr static double INFL_CONST = 0.01;
    static void set_b_g(ScalarType N, ScalarType linear, ScalarType X, ScalarType t, ScalarType &b, ScalarType &g) {
        b = -N / linear;
        // TODO is this correct? Should it not be X - right(b), x - left(b)?
        g = ScalarType(((left(X) - left(b)) * exp(linear * t) + left(b)).leftBound(),
                        ((right(X) - right(b)) * exp(linear * t) + right(b)).rightBound());
    }

    static void inflate(VectorType &v) {
        for(size_t i = 0; i < v.dimension(); ++i) {
            inflate(v[i]);
        }
    }
    static void inflate(ScalarType &x) {
        // TODO performance, numerics?
        x += ScalarType(-INFL_CONST, INFL_CONST) * (left(x) - right(x)) / 2;
    }

    static ScalarType diameter(const ScalarType &i) {
        return 2 * (right(i) - mid(i));
    }

public:
    template <typename Solver>
    static VectorType enclosure(Solver &ds,
                          ScalarType currentTime, 
                          VectorType &X,
                          ScalarType t) {
        
        using capd::min;
        using capd::max;
        std::cout.precision(16);
        
        auto vf = ds.getVectorField();
        // TODO t = from some step control policy
        // auto t  = ds.getSettedStep();

        // TODO better way
        VectorType linear_ = vf.getLinearDiagonal();
        VectorType linear(X.dimension());
        for(size_t i = 0; i < linear.dimension(); ++i) {
            linear[i] = linear_[i];
        }

        size_t firstDiss = 0;
        while(firstDiss < X.dimension() && linear[firstDiss] > -0.01 ) { ++firstDiss; }
        VectorType Z = X;
        VectorType vfXt = vf(currentTime, X) * ScalarType(0, t.rightBound());
        for(size_t i = 0; i < firstDiss; ++i) {
            Z[i] += vfXt[i];
            inflate(Z[i]); 
        }

        bool validated = false;
        const int MAX_TRIES = 10;
        int trial = 0;
        while(!validated && ++trial < MAX_TRIES) {
            validated = true;

            // first order nondissipative

            VectorType vfZt = vf(currentTime, Z) * ScalarType(0, t.rightBound());
            for(size_t i = 0; i < firstDiss; ++i) {
                ScalarType temp = X[i] + vfZt[i];
                if(right(temp) > right(Z[i])) { // How should it be used - when right, when rightBound?
                    validated = false;
                    Z[i].setRightBound(temp.rightBound());
                    inflate(Z[i]);
                }
                if(left(temp) < left(Z[i])) { 
                    validated = false;
                    Z[i].setLeftBound(temp.leftBound());
                    inflate(Z[i]);
                }
            }

            // dissipative enclosure

            // TODO maybe some better way
            VectorType N;
            if(Z.dimension() == vf.getNonlinear().dimension()){
                N = vf.getNonlinear()(Z);
            } else {
                VectorType Z_(N.dimension());
                for(size_t i = 0; i < Z.dimension(); ++i) {
                    Z_[i] = Z[i];
                }
                for(size_t i = Z.dimension(); i < Z_.dimension(); ++i) {
                    Z_[i] = 0;
                }
                VectorType N_ = vf.getNonlinear()(Z);
                for(size_t i = 0; i < Z.dimension(); ++i) {
                    N[i] = N_[i];
                }
            }
            for(size_t i = firstDiss; i < X.dimension(); ++i) {
                ScalarType b, g;
                set_b_g(N[i], linear[i], X[i], t, b, g);
            
                if(left(X[i]) <= left(b)) {
                    Z[i].setLeftBound(X[i].leftBound());
                }
                else if(left(g) >= left(Z[i])){
                    Z[i].setLeftBound(g.leftBound());
                }
                else {
                    validated = false;
                    ScalarType templ = left(g) - diameter(ScalarType(leftBound(min(left(X[i]), left(g))), rightBound(max(right(X[i]), right(g)))));
                    Z[i].setLeftBound(max(left(b), templ).leftBound());
                    //ScalarType templ = min(left(X[i]), left(g));
                    // TODO correct?
                    //ScalarType sign = templ > 0 ? -1 : 1; // we subtract if positive, add if negative
                    //Z[i].setLeftBound(max(left(b), (1 + 0.1 * sign) * templ).leftBound());
                }

                if(right(X[i]) >= right(b)) {
                    Z[i].setRightBound(X[i].rightBound());
                }
                else if(right(g) <= right(Z[i])){
                    Z[i].setRightBound(g.rightBound());
                }
                else {
                    validated = false;
                    ScalarType tempr = right(g) + .1 * diameter(ScalarType(leftBound(min(left(X[i]), left(g))), rightBound(max(right(X[i]), right(g)))));
                    Z[i].setRightBound(min(right(b), tempr).rightBound());
                    // TODO - compare with Piotr's
                    //ScalarType tempr = max(right(X[i]), right(g));
                    //ScalarType sign = tempr > 0 ? 1 : -1; // we substract if negative, add if positive
                    //Z[i].setRightBound(min(right(b), (1 + .1 * sign) * tempr).rightBound());
                }

                if(!validated){
                    inflate(Z);
                }
            }
        }
        if(!validated){
            throw std::runtime_error("Could not enclose bounds");
        }
        ScalarType b_, g_;
        auto N_ = vf.getNonlinear()(Z);
        set_b_g(N_[1], linear[1], X[1], t, b_, g_);
        // TODO refinement 


        // We set bound for solution, TODO probably we recalculate here
        // TODO maybe some better way
        VectorType N;
        if(Z.dimension() == vf.getNonlinear().dimension()){
            N = vf.getNonlinear()(Z);
        } else {
            VectorType Z_(N.dimension());
            for(size_t i = 0; i < Z.dimension(); ++i) {
                Z_[i] = Z[i];
            }
            for(size_t i = Z.dimension(); i < Z_.dimension(); ++i) {
                Z_[i] = 0;
            }
            VectorType N_ = vf.getNonlinear()(Z);
            for(size_t i = 0; i < Z.dimension(); ++i) {
                N[i] = N_[i];
            }
        }

        for(size_t i = firstDiss; i < X.dimension(); ++i) {
            ScalarType b, g;
            set_b_g(N[i], linear[i], X[i], t, b, g);
            X[i] = g;
        }
        COUT(Z);
        return Z;
    }

    template <typename DS>
    inline static void computeEnclosureAndRemainder(DS& ds, const ScalarType& t, const VectorType& x, VectorType& out_enc, VectorType& out_rem){
        typedef ScalarType Scalar;
        const static Scalar I(TypeTraits<Scalar>::zero().leftBound(),TypeTraits<Scalar>::one().rightBound());

        // we compute an enclosure for \varphi([0,timestep],iv)
        // TODO now our dynsys has no current time
        // ds.setCurrentTime(t);
        out_enc = enclosure(ds.getVectorField(),t,x,ds.getStep());
        ds.computeRemainderCoefficients(t + I*ds.getStep(),out_enc);
        out_rem =  ds.getRemainderCoefficients()[ds.getOrder()+1]*power(ds.getStep(),ds.getOrder()+1);
    }

    template <typename DS> 
    inline static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType& x, capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_enc, capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_rem){
        throw std::logic_error("not implemented");
    }

    template<class DS>
    inline static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType& x, capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_enc, capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_rem){
        throw std::logic_error("not implemented");
    }

    template <typename Solver>
    static VectorType enclosure(Solver &ds,
                          ScalarType currentTime, 
                          const VectorType &X){
        auto X_ = X;
        return enclosure(ds, currentTime, X_, ds.getStep());
    }
    template <typename Solver>
    static VectorType enclosure(Solver &ds,
                          ScalarType currentTime, 
                          const VectorType &X,
                          ScalarType t){
        auto X_ = X;
        return enclosure(ds, currentTime, X_, t);
    }
};