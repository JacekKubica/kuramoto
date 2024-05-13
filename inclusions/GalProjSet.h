#pragma once
#include "capd/capdlib.h"
#include "GalProjSolver.h"

class GalProjSet {
public:
    typedef capd::C0Rect2Set MainModesSetType;
    typedef capd::IVector VectorType;

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

    // getCurrentSet

    void move(GalProjSolver &solver) {
        solver.calculateFullEnclosure(mainModesSet.getCurrentTime(), fullVector);
        mainModesSet.move(solver);
        // TODO moving dissipative modes
        // TODO a good way to intersect?
    }

private:
    VectorType fullVector;
    MainModesSetType mainModesSet;
};