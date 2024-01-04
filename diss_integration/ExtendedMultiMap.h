#include "capd/capdlib.h"
#include "capd/diffIncl/MultiMap.h"
#include <vector>
#include <algorithm>


template<typename FMapT, typename GMapT = FMapT>
class ExtendedMultiMap : public capd::diffIncl::MultiMap<FMapT, GMapT> {
public:
    typedef typename FMapT::VectorType VectorType;
    typedef typename VectorType::size_type size_type;
    ExtendedMultiMap(FMapT &field, GMapT &perturb, FMapT &linear, FMapT &nonlinear) 
        : capd::diffIncl::MultiMap<FMapT, GMapT>(field, perturb), linear(&linear), nonlinear(&nonlinear) {
        }
    
    virtual const FMapT& getNonlinear() const { return *nonlinear; }

    VectorType diagonalCoefficients() const {
        VectorType v(linear->dimension());
        VectorType ret(linear->dimension());
        for(int i = 0; i < v.dimension(); ++i) {
            v[i] = 1; ret[i] = (*linear)(v)[i]; v[i] = 0;
        }
        return ret;
    }
    
private:
    FMapT *linear, *nonlinear;
};