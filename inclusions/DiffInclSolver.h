#ifndef __DIFF_INCL_SOLVER_H__
#define __DIFF_INCL_SOLVER_H__
#include "capd/capdlib.h"
#include "capd/dynsys/OdeSolver.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "DiffInclVectorField.h"
#include "DiffInclEnclosurePolicy.h"

// implementuj po DynSys, rzuc wyjatek na zbednych metodach

// EnclosurePolicy, StepControlPolicy
typedef capd::IMatrix MatrixTypeOfDiffInclSolver;
class MDiffInclSolver : public capd::dynsys::DynSys<MatrixTypeOfDiffInclSolver> {
public:
    typedef MatrixTypeOfDiffInclSolver MatrixType;
    typedef capd::interval ScalarType;
    typedef ScalarType::BoundType BoundType;
    typedef capd::IVector VectorType;
    typedef capd::IMap MapType;
    typedef DiffInclusionEnclosurePolicy<MDiffInclSolver> EnclosurePolicy;

    // Legacy reasons
    virtual VectorType Phi(const ScalarType& t, const VectorType &iv) override {
        throw std::logic_error("Not implemented");
        return VectorType();
    }
    virtual MatrixType JacPhi(const ScalarType& t, const VectorType &iv) override {
        throw std::logic_error("Not implemented");
        return MatrixType();
    }
    virtual VectorType Remainder(const ScalarType& t, const VectorType &iv, VectorType &out_enc) override {
        throw std::logic_error("Not implemented");
        return VectorType();
    }
    virtual VectorType enclosure(const ScalarType& t, const VectorType& x) override {
        throw std::logic_error("Not implemented");
        return VectorType();
    }

    void encloseC0Map(
        const ScalarType& t,  //< @param[in] current time of ODE
        const VectorType& x0, //< @param[in] an internal point of x, usually center of x
        const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
        VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
        VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
        VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
        MatrixType& o_jacPhi  //< @param[out] bound for derivative Dphi(x)
    ) override; // phi(t + h, x) \subset o_phi + o_jacPhi * (x - x0) + o_rem, o_enc - enclosure, h - time step
     // o_phi dla selektora, o_jacPhi dla selektora, o_rem dla selektora + wpływ inkluzji
    /* data */
    
    void encloseC0MapExternal(
        const ScalarType& t,  
        const VectorType& x0, 
        const VectorType& x,  
        VectorType& o_phi,    
        VectorType& o_rem,    
        const VectorType& enc,    
        MatrixType& o_jacPhi  
    ); 

    virtual VectorType perturbationsEffects(ScalarType currentTime,
                                            const VectorType &initial,
                                            const VectorType &W1,
                                            const VectorType &W2) = 0;

    void setStep(BoundType h) {
        selectorSolver.setStep(h);
    }
    ScalarType getStep() const {
        return selectorSolver.getStep();
    }
    ScalarType getCurrentTime() const {
        return selectorSolver.getCurrentTime();
    }
    MDiffInclSolver(const DiffInclVectorField &vf) : vf(vf),
                                                    selectorSolver(this->vf.getSelector(), 20) {}

    DiffInclVectorField& getVectorField() { return vf; }

    void setAbsoluteTolerance(BoundType tol) {
        selectorSolver.setAbsoluteTolerance(tol);
    }

    void setRelativeTolerance(BoundType tol) {
        selectorSolver.setRelativeTolerance(tol);
    }

    BoundType getAbsoluteTolerance() {
        return selectorSolver.getAbsoluteTolerance();
    }

    BoundType getRelativeTolerance() {
        return selectorSolver.getRelativeTolerance();
    }
protected:
    DiffInclVectorField vf;
    capd::dynsys::OdeSolver<MapType,
                            capd::dynsys::ILastTermsStepControl,
                            capd::dynsys::FirstOrderEnclosure> selectorSolver;
};

#include "DiffInclSolver.cpp"
#endif