#include "DiffInclSolverCW.h"
#include "capd/poincare/TimeMap.h"
using namespace capd;

void node_full(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i){
        int k = i + 1;
        out[i] = k * k * (1 - param[0] * k * k) * in[i];
        for(int j = 0; j < i; ++j){
            out[i] -= k * in[j] * in[i - 1 - j];
        }
        for(int j = 0; j + k < dimIn; ++j){
            out[i] += 2 * k * in[j] * in[j + k];
        }
        out[i] += param[i + 1];
    }
}

void node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    for(int i = 0; i < dimOut; ++i) {
        out[i] += param[i];
    }
}

const int m = 8, M = 16;
void recalc_perturb(IMap &selector, IMap &perturb, IVector enclosure) {
    IVector out(M); for(int i = 0; i < M; ++i) out[i] = 0;
    for(int i = m; i < M; ++i) { 
        int k = i + 1;
        for(int j = 0; j < i; ++j){
            out[i] -= k * enclosure[j] * enclosure[i - 1 - j];
        }
        for(int j = 0; j + k < M; ++j){
            out[i] += 2 * k * enclosure[j] * enclosure[j + k];
        }
        selector.setParameter(i, -out[i]);
        perturb.setParameter(i, out[i]);
    }
}

int main(int argc, char const *argv[])
{
    IMap ks(node_full, M, M, M + 1);
    ks.setParameter(0, .1);
    IMap perturb(node_perturb, M, M, M);
    for(int i = 0; i < M; ++i) {
        ks.setParameter(i + 1, 0);
        perturb.setParameter(i, 0);
    } 

    DiffInclSolverCW solver({ks, perturb});
    //capd::dynsys::OdeSolver<capd::IMap,
                            //capd::dynsys::ILastTermsStepControl,
                            //capd::dynsys::FirstOrderEnclosure> solver(ks, 20);
    IVector c(M); for(size_t i = 0; i < M; i++) { c[i] = interval(-1, 1); }
    c[0] = 1;
    
    C0Rect2Set s(c);
    
    IOdeSolver solver2(ks, 20);
    auto s2 = s;

    // TODO
    //poincare::TimeMap<DiffInclSolverCW> timeMap(solver);
    //auto res = timeMap(.1, s);

    const int STEPS = 10;
    for(int i = 0; i < STEPS; ++i) {
        s.move(solver);
        solver2.setStep(solver.getStep());
        s2.move(solver2);
    }


    COUT(IVector(s));
    COUT(solver.getStep());
    COUT(IVector(s2))
    COUT(solver2.getStep());
    
    return 0;
}
