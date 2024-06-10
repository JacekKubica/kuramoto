#pragma once
#include "DiffInclSolver.h"


template <typename VectorField,
          typename EnclosurePolicy>
void MDiffInclSolver<VectorField, EnclosurePolicy>::encloseC0Map(
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
    VectorType W1(o_enc.dimension());
    selectorSolver.encloseC0Map(t, x0, x, o_phi, o_rem, W1, o_jacPhi);
    // TODO
    VectorType throwaway(x);
    o_enc = EnclosurePolicy::enclosure(*this, t, throwaway, getSettedStep());
    o_rem += perturbationsEffects(t, x, W1, o_enc);
}

template <typename VectorField,
          typename EnclosurePolicy>
void MDiffInclSolver<VectorField, EnclosurePolicy>::encloseC0MapExternal(
    const ScalarType& t,  
    const VectorType& x0, 
    const VectorType& x,  
    VectorType& o_phi,    
    VectorType& o_rem,    
    const VectorType& enc,    
    MatrixType& o_jacPhi  
) {
    // NOTE: We do first order enclosure here, it also sets the step.
    //       This way foe for full system with set step should also succeed
    VectorType W1(enc.dimension());
    selectorSolver.encloseC0Map(t, x0, x, o_phi, o_rem, W1, o_jacPhi);
    o_rem += perturbationsEffects(t, x, W1, enc);
}