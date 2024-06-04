#include "DiffInclSolverCW.h"
#include "capd/poincare/TimeMap.h"
using namespace capd;

capd::interval eps = interval(-.01, .01);

double factorials[100];

interval getSinTaylorCoeff(int i, const interval &time, const interval &h) {
    int sign = (i % 4 == 2 || i % 4 == 3) ? -1 : 1; 
    return power(h, i) * (i % 2 ? sign * cos(time) : sign * sin(time)) / factorials[i];
}


void fullMapNode(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    out[0] = in[1];
    out[1] = -in[0] + param[0] * sin(time);
}
void partMapNode(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    out[0] = in[1];
    out[1] = -in[0];
    capd::autodiff::Node cursteppow(1);
    capd::autodiff::Node step = time - param[noParam - 1];
    for(int i = 0; i < noParam - 1; ++i) {
        out[1] += cursteppow * param[i];
        cursteppow *= step;
    }
}
void node_perturb(capd::autodiff::Node& time, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node param[], int noParam){
    out[0] = capd::autodiff::Node(0);
    out[1] = param[0];
}

int main(int argc, char const *argv[])
{ 
    factorials[0] = 1;
    for(int i = 1; i < 100; ++i) {
        factorials[i] = i * factorials[i - 1];
    }

    const int M = 2;
    IMap fullMap(fullMapNode, M, M, 1);
    fullMap.setParameter(0, eps);
    const int MAX_DEG = 10; 
    IMap partMap(partMapNode, M, M, MAX_DEG + 1);
    IMap perturb(node_perturb, M, M, M);

    IMultiMap mm(partMap, perturb);
    CWDiffInclSolver cwDiffInclSolver(mm, 20, IMaxNorm());


    for(int i = 0; i < MAX_DEG + 1; ++i) {
        partMap.setParameter(i, 0);
    }
    for(int i = 0; i < M; ++i) {
        perturb.setParameter(i, 0);
    } 

    DiffInclSolverCW<> solver({partMap, perturb});
    //capd::dynsys::OdeSolver<capd::IMap,
                            //capd::dynsys::ILastTermsStepControl,
                            //capd::dynsys::FirstOrderEnclosure> solver(ks, 20);
    IVector c(M); for(size_t i = 0; i < M; i++) { c[i] = 0; }
    c[0] = 1;

    C0Rect2Set s(c);
    InclRect2Set cws(c);

    IOdeSolver solver2(fullMap, 20);
    solver2.setStep(1./64);
    C0Rect2Set s2(c);

    const int STEPS = 1;
    const int DEG = 0;
    for(int i = 0; i < STEPS; ++i) {
        s2.move(solver2);
        
        solver.setStep(leftBound(solver2.getStep()));        
        auto time = solver2.getCurrentTime();

        auto step = leftBound(solver2.getStep());
        auto &selector = solver.getVectorField().getSelector();
        selector.setParameter(0, 0);
        selector.setParameter(MAX_DEG, s.getCurrentTime());
        for(int i = 0; i < DEG; ++i) {
            selector.setParameter(i, eps * getSinTaylorCoeff(i, left(time), 1)); // 1 because we are gonna multiply by h^i in map
        }
        auto val = eps * getSinTaylorCoeff(DEG, time, interval(0, step));

        selector.setParameter(0, selector.getParameter(0) + val.mid());
        solver.getVectorField().getPerturbation().setParameter(0, val - val.mid());
        s.move(solver);
    }
    COUT("=");
    COUT(DEG);
    COUT(IVector(s));
    COUT(IVector(s2));
    COUT(maxDiam(IVector(s)));
    COUT(maxDiam(IVector(s2)));
    //COUT(IVector(cws))
    
    return 0;
}
