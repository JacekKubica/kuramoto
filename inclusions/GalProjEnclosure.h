#pragma once
#include "capd/capdlib.h"


template <class Solver>
class DissipativeEnclosure {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;
    typedef Solver::VectorFieldType VectorFieldType;

    static VectorType enclosure(Solver &ds,
                          ScalarType currentTime,
                          const VectorType &X) {
        
        size_t firstDiss = 0;
        while(firstDiss < X.dimension() && vf.getLinearDiagonal()[firstDiss] > -0.01 ) { ++firstDiss; }
        Vector Z = X;
        Vector vfXt = vf(X) * Scalar(0, t.rightBound());
        for(size_t i = 0; i < firstDiss; ++i) {
            Z[i] += vfXt[i];
            inflate(Z[i]); 
        }

        bool validated = false;
        while(!validated) {
            validated = true;

            // first order nondissipative
            Vector vfZt = vf(Z) * Scalar(0, t.rightBound());
            for(size_t i = 0; i < firstDiss; ++i) {
                Scalar temp = X[i] + vfZt[i];
                if(right(temp) > right(Z[i])) { // How should it be used - when right, when rightBound?
                    validated = false;
                    Z[i].setRightBound(temp.rightBound());
                }
                if(left(temp) < left(Z[i])) { 
                    validated = false;
                    Z[i].setLeftBound(temp.leftBound());
                }
                if(!validated){
                    for(auto &z: Z) {
                        inflate(z);
                    }
                }
            }


            // dissipative enclosure
            Vector N = vf.getNonlinear()(Z);
            for(size_t i = firstDiss; i < X.dimension(); ++i) {
                Scalar b = -N[i] / vf.linear[i];
                Scalar X_ = X[i];
                Scalar g = Scalar(((left(X_) - left(b)) * exp(vf.linear[i] * t) + left(b)).leftBound(),
                                  ((right(X_) - right(b)) * exp(vf.linear[i] * t) + right(b)).rightBound());
            
                COUT(t); COUT(b); COUT(g); COUT(Z);
                if(left(X[i]) <= left(b)) {
                    Z[i].setLeftBound(X[i].leftBound());
                }
                else if(left(g) >= left(Z[i])){
                    Z[i].setLeftBound(g.leftBound());
                }
                else {
                    validated = false;
                    Scalar templ = min(left(X[i]), left(g));
                    Scalar sign = templ > 0 ? -1 : 1; // we subtract if pos, add if negative
                    Z[i].setLeftBound(max(left(b), (1 + 0.1 * sign) * templ).leftBound());
                }
                if(right(X[i]) >= right(b)) {
                    Z[i].setRightBound(X[i].rightBound());
                }
                else if(right(g) <= right(Z[i])){
                    Z[i].setRightBound(g.rightBound());
                }
                else {
                    validated = false;
                    // TODO - compare with Piotr's
                    Scalar tempr = max(right(X[i]), right(g));
                    Scalar sign = tempr > 0 ? -1 : 1;
                    Z[i].setRightBound(min(right(b), (1 + .1 * sign) * tempr).rightBound());
                }
            }
        }
        // TODO refinement 
        return Z;
    }
};