#pragma once
#include "capd/capdlib.h"


class DiffInclusionEnclosurePolicy {
public:
    typedef capd::interval ScalarType;
    typedef capd::IVector VectorType;
    typedef capd::IMatrix MatrixType;
    typedef capd::IMap MapType;

    template <class DiffInclSolver>
    static VectorType enclosure(DiffInclSolver &ds,
                          ScalarType currentTime,
                          const VectorType &X) {
        
       return capd::dynsys::FirstOrderEnclosure::enclosure(ds.getVectorField(), currentTime, X, ds.getStep()); 
    
    }
};