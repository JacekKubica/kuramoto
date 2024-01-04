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

    // TODO outside of class
    void computeBounds(size_type firstDiss, ScalarType time,
                       const VectorType &initial, const VectorType &enclosureCandidate,
                       VectorType &isolationBoundary, VectorType &solutionBound) const { // results 

        auto N = vf.getNonlinear()(enclosureCandidate);
        auto lambdas = vf.diagonalCoefficients();
        for (int i = firstDiss; i < initial.dimension(); ++i) {
            isolationBoundary[i] = N[i] / (-lambdas[i]);
            solutionBound[i] = 
                ScalarType(
                    ((initial[i].leftBound() - isolationBoundary[i].leftBound()) 
                     * exp(lambdas[i] * time) + isolationBoundary[i].leftBound()).leftBound(),
                    ((initial[i].rightBound() - isolationBoundary[i].rightBound()) 
                     * exp(lambdas[i] * time) + isolationBoundary[i].rightBound()).rightBound());
        }
    }

    void tryDissipativeEnclosure(size_type firstDiss, const std::vector<bool> &validated, 
                                 const VectorType &initial,
                                 const VectorType &isolationBoundary, const VectorType &solutionBound,
                                 VectorType &enclosureCandidate) const {
        for(int i = firstDiss; i < isolationBoundary.dimension(); ++i){
            if(!validated[i]){
                ScalarType w;
                intervalHull(initial[i], solutionBound[i]).split(w);
                w *= 2;
                if(isolationBoundary[i].rightBound() > enclosureCandidate[i].rightBound())
                    enclosureCandidate[i].
                        setRightBound(capd::min(isolationBoundary[i].rightBound(), 
                                          (solutionBound[i] + RESERVE_CONSTANT * w).rightBound()));
                if(isolationBoundary[i].leftBound() < enclosureCandidate[i].leftBound())
                    enclosureCandidate[i].
                        setLeftBound(capd::max(isolationBoundary[i].leftBound(),
                                     (solutionBound[i] - RESERVE_CONSTANT * w).leftBound()));
            }
            inflate(enclosureCandidate[i]);
        }
    }
    
    void checkDissEnclosure(size_type firstDiss, 
                            const VectorType &initial, const VectorType &enclosureCandidate, 
                            const VectorType &isolationBoundary, const VectorType &solutionBound,
                            std::vector<bool> &validated) const { // results
        for(int i = firstDiss; i < initial.dimension(); ++i){
            if(initial[i].rightBound() < isolationBoundary[i].rightBound() 
                && enclosureCandidate[i].rightBound() < solutionBound[i].rightBound()){
                validated[i] = false;
            }
            else if(initial[i].leftBound() > isolationBoundary[i].leftBound() 
                    && enclosureCandidate[i].leftBound() > solutionBound[i].leftBound()){
                validated[i] = false;
            }
            else{
                validated[i] = true;
            }
        }
    }
    
    void tryEncloseNondiss(size_type firstDiss, ScalarType time,
                           const VectorType &initial,
                           VectorType &enclosureCandidate, std::vector<bool> &validated) const { // results
        auto vfCandidate = vf(enclosureCandidate);
        for(int i = 0; i < firstDiss; ++i) {
            auto y = initial[i] + capd::interval(0., time.rightBound()) * vfCandidate[i];
            if(!y.subsetInterior(enclosureCandidate[i])){
                validated[i] = false;
                enclosureCandidate[i] = intervalHull(y, enclosureCandidate[i]);
                inflate(enclosureCandidate[i]);
            }
            else {
				# TODO is it okay?
				enclosureCandidate[i] = y;
                validated[i] = true;
			}
        }
    }

    void refine(size_type firstDiss,
                const VectorType &initial, const VectorType &isolationBoundary, const VectorType &solutionBound,
                VectorType &enclosureCandidate) const { // results
        for(int i = firstDiss; i < initial.dimension(); ++i){
            if(isolationBoundary[i].rightBound() <= initial[i].rightBound())
                enclosureCandidate[i].setRightBound(initial[i].rightBound());
            else if(isolationBoundary[i].rightBound() > initial[i].rightBound())
                enclosureCandidate[i].setRightBound(solutionBound[i].rightBound());
            if(isolationBoundary[i].leftBound() >= initial[i].leftBound())
                enclosureCandidate[i].setLeftBound(initial[i].leftBound());
            else if(isolationBoundary[i].leftBound() < initial[i].leftBound())
                enclosureCandidate[i].setLeftBound(solutionBound[i].leftBound());
        }
    }

    const double RESERVE_CONSTANT = 0.01; // publ DG
    const double INFLATION_CONSTANT = 1.0001; // publ D2
    const double DISSIPATION_LIMIT = -0.01; 
    const double MAX_REFINMENT = 5;
};


template <class MapType>
typename MyDissipativeSolver<MapType>::VectorType MyDissipativeSolver<MapType>::enclosure(const ScalarType &t, const VectorType &initial) const {
    auto enclosureCandidate = initial; 
    for(int i = 0; i < initial.dimension(); ++i) 
        enclosureCandidate[i] += capd::TypeTraits<ScalarType>::epsilon() * capd::interval(-1, 1); // for point data
    
    // we do first order enclosure for non-dissipative modes
    size_type firstDiss = 0;
    auto vfInitial = vf(initial);
    while(firstDiss < initial.dimension() && vf.diagonalCoefficients()[firstDiss] > DISSIPATION_LIMIT) {
        enclosureCandidate[firstDiss] += capd::interval(0., t.rightBound()) * vfInitial[firstDiss];
        inflate(enclosureCandidate[firstDiss]);
        ++firstDiss;
    }

    VectorType isolationBoundary(initial.dimension()); // publ b
    VectorType solutionBound(initial.dimension()); // publ g
    computeBounds(firstDiss, t, initial, enclosureCandidate, isolationBoundary, solutionBound);
    
    std::vector<bool> validated(initial.dimension(), false);
    tryDissipativeEnclosure(firstDiss, validated, initial, isolationBoundary, solutionBound, enclosureCandidate);
    bool allValidated = false;
    while(!allValidated) {
        tryEncloseNondiss(firstDiss, t, initial, enclosureCandidate, validated);
        computeBounds(firstDiss, t, initial, enclosureCandidate, isolationBoundary, solutionBound);
        checkDissEnclosure(firstDiss, initial, enclosureCandidate, isolationBoundary, solutionBound, validated);
        allValidated = std::all_of(validated.begin(), validated.end(), [](bool v){return v;});
        if(allValidated) {
            refine(firstDiss, initial, isolationBoundary, solutionBound, enclosureCandidate);
        } else {
            tryDissipativeEnclosure(firstDiss, validated, initial, isolationBoundary, solutionBound, enclosureCandidate);
        }
    }
    for(int i = 0; i < MAX_REFINMENT; ++i) {
		computeBounds(firstDiss, t, initial, enclosureCandidate, isolationBoundary, solutionBound);
        refine(firstDiss, initial, isolationBoundary, solutionBound, enclosureCandidate);
    }
    COUT(isolationBoundary);
    return isolationBoundary;
    return enclosureCandidate;
}