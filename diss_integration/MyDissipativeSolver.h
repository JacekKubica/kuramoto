#include "capd/capdlib.h"
#include <vector>
#include <algorithm>


// We demand that map can be evaluated and split into linear and nonlinear part
template <class MapType>
class MyDissipativeSolver {
public:
  typedef typename MapType::ScalarType ScalarType;
  typedef typename MapType::VectorType VectorType;
  typedef typename VectorType::size_type size_type; 
  MyDissipativeSolver(const MapType &vf) : vf(vf) {}
  virtual VectorType enclosure(const ScalarType &t, const VectorType &x) const;
private:
  MapType vf;

  void inflate(ScalarType &w) const {
    ScalarType temp;
    // split gets mid into w, diameter into temp
    w.split(temp);
    temp *= INFLATION_CONSTANT;
    w += temp;
  }


  void improve_isolations(size_type firstDiss,
              VectorType &isolationBoundary) const {
    auto N = vf.evalNonlinearPlusPerturb(isolationBoundary);
    auto lambdas = vf.diagonalCoefficients();
    for (int i = firstDiss; i < isolationBoundary.dimension(); ++i) {
      isolationBoundary[i] = N[i] / (-lambdas[i]);
    }
  }

  void computeBounds(size_type firstDiss, ScalarType time,
             const VectorType &initial, const VectorType &enclosureCandidate,
             VectorType &isolationBoundary, VectorType &solutionBound) const;

  void tryDissipativeEnclosure(size_type firstDiss, const std::vector<bool> &validated, 
                 const VectorType &initial,
                 const VectorType &isolationBoundary, const VectorType &solutionBound,
                 VectorType &enclosureCandidate) const;


  void checkDissEnclosure(size_type firstDiss, 
              const VectorType &initial, const VectorType &enclosureCandidate, 
              const VectorType &isolationBoundary, const VectorType &solutionBound,
              std::vector<bool> &validated) const;
  

  void tryEncloseNondiss(size_type firstDiss, ScalarType time,
               const VectorType &initial,
               VectorType &enclosureCandidate, std::vector<bool> &validated) const;


  void refine(size_type firstDiss,
        const VectorType &initial, const VectorType &isolationBoundary, const VectorType &solutionBound,
        VectorType &enclosureCandidate) const;
  bool validate_far_tail(const VectorType &initial, const VectorType &isolationBoundary,
               const VectorType &enclosureCandidate);

  void update_tail(bool far_tail_validated, // decides course of action
          const VectorType &initial,
          VectorType &enclosureCandidate // return
          );


  const double RESERVE_CONSTANT = 0.01; // publ DG
  const double INFLATION_CONSTANT = 1.0001; // publ D2
  const double INFLATION_CONSTANT2 = 1.01 // publ also D2
  const double DISSIPATION_LIMIT = -0.01; 
  const double MAX_REFINMENT = 5;
};
