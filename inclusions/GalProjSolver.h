#pragma once
#include "DiffInclSolverCW.h"
#include "GalProjVectorField.h"
#include "GalProjEnclosure.h"


typedef capd::IMatrix MatrixTypeOfDiffInclSolver;
class GalProjSolver : public capd::dynsys::DynSys<MatrixTypeOfDiffInclSolver> {
public:
    typedef MatrixTypeOfDiffInclSolver MatrixType;
    typedef capd::interval ScalarType;
    typedef ScalarType::BoundType BoundType;
    typedef capd::IVector VectorType;
    typedef capd::IMap MapType;
    typedef DissipativeEnclosure EnclosurePolicy;
    typedef GalProjVectorField VectorFieldType;
    typedef DiffInclSolverCW<GalProjVectorField, EnclosurePolicy> InclusionSolverType;

    GalProjSolver(VectorFieldType vf) : baseSolver(vf) {}

    void encloseC0Map(
        const ScalarType& t,  //< @param[in] current time of ODE
        const VectorType& x0, //< @param[in] an internal point of x, usually center of x
        const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
        VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
        VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
        VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
        MatrixType& o_jacPhi  //< @param[out] bound for derivative Dphi(x)
    ) override {
        if(!full_enclosure_ptr) {
            throw std::logic_error("Full enclosure was not calculated yet");
        }
        for(size_t i = 0; i < o_enc.dimension(); ++i) {
            o_enc[i] = this->full_enclosure[i];
        }
        //baseSolver.encloseC0Map(t, x0, x, o_phi, o_rem, o_enc, o_jacPhi);
        baseSolver.encloseC0MapExternal(t, x0, x, o_phi, o_rem, o_enc, o_jacPhi);
    }

    void calculateFullEnclosure(const ScalarType &t, VectorType &o_x) {
        this->full_enclosure = EnclosurePolicy::enclosure(*this, t, o_x);
        this->full_enclosure_ptr = &this->full_enclosure;
        //this->baseSolver.getVectorField().setPerturbation(this->full_enclosure);
        this->baseSolver.getVectorField().setPerturbation(this->full_enclosure);
    }

    void setStep(BoundType h) {
        baseSolver.setStep(h);
        settedStep = h;
    }

    ScalarType getSettedStep() {
        return settedStep;
    }

    ScalarType getStep() {
        return baseSolver.getStep();
    }

    VectorFieldType& getVectorField() {
        return baseSolver.getVectorField();
    }

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
private:
    InclusionSolverType baseSolver;
    VectorType full_enclosure;
    VectorType *full_enclosure_ptr = nullptr;
    ScalarType settedStep;
};
