#ifndef _CAPD_DIFFINCL_DIFFINCLUSIONDISS_H_
#define _CAPD_DIFFINCL_DIFFINCLUSIONDISS_H_

#include "DiffInclusionCW.h"
#include <vector>
#include <algorithm>
#include <iomanip>

namespace capd{
namespace diffIncl{

template<typename MapT, typename DissOp, int m, typename DynSysT = capd::dynsys::OdeSolver< typename MapT::MapType>>
class DiffInclusionDiss : public capd::diffIncl::DiffInclusionCW<MapT, DynSysT>{
public:
    typedef capd::diffIncl::DiffInclusionCW<MapT, DynSysT> BaseClass;
    typedef typename BaseClass::MultiMapType   MultiMapType;
    typedef typename BaseClass::MapType        MapType;
    typedef typename BaseClass::FunctionType   FunctionType;
    typedef typename BaseClass::MatrixType     MatrixType;
    typedef typename BaseClass::VectorType     VectorType;
    typedef typename BaseClass::ScalarType     ScalarType;
    typedef typename BaseClass::NormType       NormType;
    typedef typename MatrixType::size_type     size_type;
    using BaseClass::getStep;

    DiffInclusionDiss(
               MultiMapType& A_diffIncl,
               size_type A_order,
               NormType const & A_norm,
               DissOp & diss_op,
               ScalarType const & expErrorTolerance = capd::TypeTraits<ScalarType>::epsilon()
    )
      : BaseClass(A_diffIncl, A_order, A_norm, expErrorTolerance),
        diss_op(diss_op),
        z(m), b(m), g(m), kvalidated {std::vector<bool>(m, false)}{}

    void inflate(ScalarType &w){
        ScalarType temp;
        w.split(temp);
        temp *= D2;
        w += temp;
    }

    auto w(int i, const VectorType &x){
        ScalarType temp;
        intervalHull(x[i], g[i]).split(temp);
        return 2 * temp.rightBound();
    }

    void recompute_b_and_g(const VectorType &x){
        auto N = diss_op.N(z); // poly bounds, will be replaced by this thing below if we are dealing with proj
        for (int i = 0; i < m; ++i) {
            b[i] = N[i] / (-diss_op.lambda(i));
            g[i] = ScalarType(((x[i].leftBound() - b[i].leftBound()) * exp(diss_op.lambda(i) * getStep()) + b[i].leftBound()).leftBound(),
                              ((x[i].rightBound() - b[i].rightBound()) * exp(diss_op.lambda(i) * getStep()) + b[i].rightBound()).rightBound());
        }
    }

    VectorType diffInclusionEnclosure(const ScalarType &time, const VectorType &x) override {
        auto vf = BaseClass::getVectorField();
        // auto vf = BaseClass::getVectorField().getVectorField(); <- it was vf without perturbation
        z = x;
        for(int i = 0; i < m; ++i) z[i] += capd::TypeTraits<ScalarType>::epsilon() * capd::interval(-1, 1); // for point data
        int first_diss = 0;
        auto vfx = vf(x);
        while(diss_op.lambda(first_diss) > -0.01){
            z[first_diss] += capd::interval(0., getStep().rightBound()) * vfx[first_diss];
            inflate(z[first_diss]);
            first_diss += 1;
        }
        recompute_b_and_g(x); // HERE only x because z is class variable!
        for(int i = first_diss; i < m; ++i){
            if(z[i].leftBound() > b[i].leftBound())
                z[i].setLeftBound(max(b[i].leftBound(), g[i].leftBound() - DG * w(i, x)));
            if(z[i].rightBound() < b[i].rightBound())
                z[i].setRightBound(min(b[i].rightBound(), g[i].rightBound() + DG * w(i, x)));
        }
        bool validated = false;
        int steps = 0;

        while(!validated){ //TODO some exception when we cannot validate
            ++steps;
            auto vfz = vf(z);
            for(int i = 0; i < first_diss; ++i){
                auto yi = x[i] + capd::interval(0., getStep().rightBound()) * vfz[i];
                if(!yi.subsetInterior(z[i])){
                    kvalidated[i] = false;
                    z[i] = intervalHull(yi, z[i]);
                    inflate(z[i]);
                }
                else
                    kvalidated[i] = true;
            }
            recompute_b_and_g(x);
            for(int i = first_diss; i < m; ++i){
                if(x[i].rightBound() < b[i].rightBound() && z[i].rightBound() < g[i].rightBound()){
                    kvalidated[i] = false;
                }
                else if(x[i].leftBound() > b[i].leftBound() && z[i].leftBound() > g[i].leftBound()){
                    kvalidated[i] = false;
                }
                else{
                    kvalidated[i] = true;
                }
            }

            validated = std::all_of(kvalidated.begin(), kvalidated.end(), [](bool v){return v;});

            if(validated){
                // TODO refine nondiss
                for(int i = first_diss; i < m; ++i){
                    if(b[i].rightBound() <= x[i].rightBound())
                        z[i].setRightBound(x[i].rightBound());
                    else if(b[i].rightBound() > x[i].rightBound())
                        z[i].setRightBound(g[i].rightBound());
                    if(b[i].leftBound() >= x[i].leftBound())
                        z[i].setLeftBound(x[i].leftBound());
                    else if(b[i].leftBound() < x[i].leftBound())
                        z[i].setLeftBound(g[i].leftBound());
                }
            }
            else{
                for(int i = first_diss; i < m; ++i){
                    if(!kvalidated[i]){
                        if(b[i].rightBound() > z[i].rightBound())
                            z[i].setRightBound(min(b[i].rightBound(), g[i].rightBound() + DG * w(i, x)));
                        if(b[i].leftBound() < z[i].leftBound())
                            z[i].setLeftBound(max(b[i].leftBound(), g[i].leftBound() - DG * w(i, x)));
                    }
                    inflate(z[i]);
                }
            }
        }

        const int max_refinement = 5;
        for(int iter = 0; iter < max_refinement; ++iter){
            auto vfz = vf(z);
            for(int i = 0; i < first_diss; ++i){
                z[i] = x[i] + capd::interval(0., getStep().rightBound());
            }
            recompute_b_and_g(x);
            for(int i = first_diss; i < m; ++i){
                if(b[i].rightBound() <= x[i].rightBound())
                    z[i].setRightBound(x[i].rightBound());
                else if(b[i].rightBound() > x[i].rightBound())
                    z[i].setRightBound(g[i].rightBound());
                if(b[i].leftBound() >= x[i].leftBound())
                    z[i].setLeftBound(x[i].leftBound());
                else if(b[i].leftBound() < x[i].leftBound())
                    z[i].setLeftBound(g[i].leftBound());
            }
        }

        // std::cout << "===================" << '\n';
        // std::cout << "WE DO DISS ENCLOSURE" << '\n';
        // try{
        //     auto bd = BaseClass::diffInclusionEnclosure(time, x);
        //     std::cout << std::setprecision(9);
        //     std::cout << bd << '\n';
        //     std::cout << "***************" << '\n';
        //     std::cout << z << '\n';
        // } catch(...){
        //     std::cout << "we cannot do diff incl enclosure" << '\n';
        // }
        return z;
    }

private:

    DissOp &diss_op;
    VectorType z, b, g;
    std::vector<bool> kvalidated;
    const double DG = 0.01, D2 = 1.0001;
};
}
}

#endif
