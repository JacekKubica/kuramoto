#pragma once
#include "capd/capdlib.h"
#include "GalProjSolver.h"

class GalProjSet {
public:
    typedef capd::C0Rect2Set MainModesSetType;
    typedef capd::IVector VectorType;
    typedef VectorType::ScalarType ScalarType;

    GalProjSet(const VectorType &fullVector, const int projectionDimension) 
        : fullVector(fullVector),
          mainModesSet(VectorType()) // due to no default constructor in set
    {
        VectorType temp(projectionDimension);
        for(int i = 0; i < projectionDimension; ++i) {
            temp[i] = fullVector[i];
        }
        mainModesSet = MainModesSetType(temp);
    }
    
    ScalarType getCurrentTime() {
        return mainModesSet.getCurrentTime();
    }

    // getCurrentSet

    void move(GalProjSolver &solver) {
        // TODO save all dissipative to fullVector
        solver.calculateFullEnclosure(mainModesSet.getCurrentTime(), fullVector); // this sets perturbations, moves dissipative modes
        mainModesSet.move(solver);
        auto mainModesVector = VectorType(mainModesSet);
        for(size_t i = 0; i < mainModesVector.dimension(); ++i) {
            // TODO intersect on dissipative
            fullVector[i] = mainModesVector[i];
        }
    }

public:
    VectorType fullVector;
    MainModesSetType mainModesSet;
};