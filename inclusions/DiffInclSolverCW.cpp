#pragma once
#include "DiffInclSolverCW.h"

template <typename VectorField,
          typename EnclosurePolicy>
typename DiffInclSolverCW<VectorField, EnclosurePolicy>::VectorType 
DiffInclSolverCW<VectorField, EnclosurePolicy>::perturbationsEffects(
                                          ScalarType currentTime,
                                          const VectorType &initial,
                                          const VectorType &W1,
                                          const VectorType &W2) {
                   
    MatrixType J = Base::vf.derivative(currentTime, W2);
    VectorType delta = Base::vf.getPerturbation()(currentTime, W1); 

    VectorType C = rightVector(delta);

    size_t i, j;
    for( i=0; i < J.numberOfRows(); ++i)
      for( j = 0; j < J.numberOfColumns(); ++j)
        if(i != j)
          J[i][j] = (capd::abs(J[i][j])).right();
        else
          J[i][j] = (J[i][j]).right();
    MatrixType At = Base::getStep() * J;
    MatrixType A = MatrixType::Identity(J.numberOfRows());
    MatrixType Sum = A;
    Norm norm; // TODO should we create this norm once somewhere etc?
    ScalarType AtNorm = right(norm(At)),
               AnNorm = right(norm(A));

    ScalarType n = typename capd::TypeTraits<ScalarType>::Real(2.0);  // n = i + 2
    ScalarType q = AtNorm/n;

    while(true){
      A = A * At / n;
      Sum += A;
      AnNorm *= q;
      n += ScalarType(1.0);
      q = AtNorm / n;
      if(q < 1){
        // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
        ScalarType remainder = right(AnNorm * AtNorm / (n - AtNorm ));
        if(remainder < Base::getAbsoluteTolerance())  // TODO ok?
          break;
      }
    }
    // we recompute remainder because norm of A can be smaller than approximation : AnNorm
    ScalarType remainder = right(norm(A) * AtNorm / (n  - AtNorm ));

    for(i=0; i < J.numberOfRows(); ++i)
      for(j = 0; j < J.numberOfColumns(); ++j)
        Sum[i][j] += remainder * ScalarType(-1.0, 1.0);
    VectorType D = Base::getStep() * (Sum * C);
    VectorType result(D.dimension());

    for(i=0; i< D.dimension(); ++i)
      result[i] = ScalarType(-D[i].rightBound(), D[i].rightBound());
    return result;
}