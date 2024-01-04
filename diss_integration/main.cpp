#include "MyDissipativeSolver.h"
#include "KSExtendedMultimapContainer.h"
#include "PolynomialBound.h"
using namespace capd;

void coordChange(const ExtendedMultiMap<IMap> &emm) {\
	
}

int main() {
    constexpr int M = 16, m = 8;;
    // const capd::interval niu = capd::interval(99) / 100;
    const capd::interval niu = capd::interval(75) / 100;
    capd::IVector N(M);
    capd::IVector N_(m);
    // N[0] = N_[0] = capd::interval(-0.0482822, 0.0482822);
    N[0] = N_[0] = 0.075 * capd::interval(-1, 1);
    capd::pdes::PolynomialBound<capd::interval, unsigned, M> pb(M, 2.31307e+08 * capd::interval(-1, 1), 20, N.begin());
    KSExtendedMultiMapContainer<m, M> f(niu, pb);
    MyDissipativeSolver<KSExtendedMultiMapContainer<m, M>> mds(f);
    auto x = mds.enclosure(1e-4, N_);
    COUT(x);

}