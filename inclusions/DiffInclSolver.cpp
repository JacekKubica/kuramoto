#pragma once
#include "DiffInclSolver.h"


void MDiffInclSolver::encloseC0Map(
    const ScalarType& t,  
    const VectorType& x0, 
    const VectorType& x,  
    VectorType& o_phi,    
    VectorType& o_rem,    
    VectorType& o_enc,    
    MatrixType& o_jacPhi  
) {
    // NOTE: We do first order enclosure here, it also sets the step.
    //       This way foe for full system with set step should also succeed
    VectorType effects = perturbationsEffects(t, x);
    selectorSolver.encloseC0Map(t, x0, x, o_phi, o_rem, o_enc, o_jacPhi);
    for(size_t i = 0; i < effects.dimension(); ++i) {
        o_rem[i] += effects[i];
    } 
}