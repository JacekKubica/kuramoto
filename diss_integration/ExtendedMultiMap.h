#include "capd/capdlib.h"
#include "capd/diffIncl/MultiMap.h"
#include <vector>
#include <algorithm>


template<typename FMapT, typename GMapT = FMapT, class LMapT = typename FMapT::MatrixType, class NMapT = FMapT>
class ExtendedMultiMap : public capd::diffIncl::MultiMap<FMapT, GMapT> {
public:
    typedef typename FMapT::VectorType VectorType;
    typedef typename VectorType::ScalarType ScalarType;
    typedef typename VectorType::size_type size_type;
    ExtendedMultiMap(FMapT &field,
                     GMapT &perturb,
                     LMapT &linear,
                     NMapT &nonlinear) 
        : capd::diffIncl::MultiMap<FMapT, GMapT>(field, perturb),
          linear(&linear),
          nonlinear(&nonlinear) 
        {}
    
    virtual const NMapT& getNonlinear() const { return *nonlinear; }

    VectorType diagonalCoefficients() const {
        auto dim = linear->dimension().first;
        VectorType v(dim);
        VectorType ret(dim);
        for(int i = 0; i < v.dimension(); ++i) {
            v[i] = 1; ret[i] = ((*linear) * v)[i]; v[i] = 0;
        }
        return ret;
    }

    VectorType evalLinearNondiagonal(const VectorType &X) const {
        auto dim = linear->dimension().first;
        auto linear_cpy = *linear;
        for(int i = 0; i < X.dimension(); ++i) {
            linear_cpy[i][i] = 0;
        }
        return linear_cpy * X;
    }
    
private:
    LMapT *linear; 
    NMapT *nonlinear;
};