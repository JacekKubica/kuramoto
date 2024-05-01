#pragma once
#include "capd/capdlib.h"
#include "DiffInclVectorField.h"


template <class DiffInclSolver>
class DiffInclusionEnclosurePolicy {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;

private:
    static void inflate(Scalar &x) {
        // TODO performance, numerics?
        x += Scalar(-INFL_CONST, INFL_CONST) * (left(x) - right(x)) / 2;
    }
public:
    // TODO getStep or argument?
    static VectorType enclosure(DiffInclSolver &ds,
                          ScalarType currentTime,
                          const VectorType &X) {
        
        
    
    }
    // TODO getStep or argument?
    //static VectorType enclosure(DiffInclSolver &ds,
                          //ScalarType currentTime,
                          //const VectorType &X) {

        //DiffInclVectorField vf = ds.getVectorField();
        //VectorType Z = X;
        //VectorType vfXt = vf(X) * Scalar(0, t.rightBound());
        //for(size_t i = 0; i < X.dimension(); ++i) {
            //Z[i] += vfXt[i];
            //inflate(Z[i]); 
        //}
        //vf.calcPerturbations(Z);

        //bool validated = false;
        //while(!validated) {
            //validated = true;
            //// first order nondissipative
            //Vector vfZt = vf(Z) * Scalar(0, t.rightBound());
            //for(size_t i = 0; i < X.dimension(); ++i) {
                //Scalar temp = X[i] + vfZt[i];
                //if(right(temp) > right(Z[i])) { // How should it be used - when right, when rightBound?
                    //validated = false;
                    //Z[i].setRightBound(temp.rightBound());
                //}
                //if(left(temp) < left(Z[i])) { 
                    //validated = false;
                    //Z[i].setLeftBound(temp.leftBound());
                //}
                //if(!validated){
                    //for(auto &z: Z) {
                        //inflate(z);
                    //}
                //}
            //}
            //vf.calcPerturbations(Z);
        //}
        //return Z;
    //}
};