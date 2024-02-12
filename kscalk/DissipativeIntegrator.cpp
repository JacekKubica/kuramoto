#include <iomanip>

template<int m, int M, class Scalar, class Exponent, class VectorType, class SetType, class ParamType, class DissipativeOperator, class Norm>
DissipativeIntegrator<m, M, Scalar, Exponent, VectorType, SetType, ParamType, DissipativeOperator, Norm>::DissipativeIntegrator(DissipativeIntegrator<m, M, Scalar, Exponent, VectorType, SetType, ParamType, DissipativeOperator, Norm>::BoundType h, const VectorType &x, const capd::pdes::PolynomialBound<Scalar, Exponent, M> &T0, const ParamType &p)
    : x(SetType(x)),
      tv{h, x, T0, p} {}

template<int m, int M, class Scalar, class Exponent, class VectorType, class SetType, class ParamType, class DissipativeOperator, class Norm>
bool DissipativeIntegrator<m, M, Scalar, Exponent, VectorType, SetType, ParamType, DissipativeOperator, Norm>::incl_enclosure(){
    try{
        tv.get_W2() = tv.get_diff_inclusion_cw_solver().diffInclusionEnclosure(tv.get_h(), x); // TODO set_W2?
    } catch(capd::dynsys::SolverException<VectorType> &e){
        return false;
    }
    return true;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class SetType, class ParamType, class DissipativeOperator, class Norm>
bool DissipativeIntegrator<m, M, Scalar, Exponent, VectorType, SetType, ParamType, DissipativeOperator, Norm>::enclosure_with_tail(){
    const int max_iter = M / 2, maxdcount = 4;
    bool validated = false;
    auto W2_initial = tv.get_W2();
    auto T_initial = tv.get_T();
    for(int dcount = 0; !validated && dcount < maxdcount; ++dcount){
        tv.get_W2() = W2_initial;
        tv.get_T() = T_initial;

        if(dcount) { tv.decpower(dcount); }
        if(tv.validate_far_tail(true)){ // TODO should we update here or not? true yes false nope
            tv.validate_tail(true);
            for(int i = 0; !validated && i < max_iter; ++i){
                if(incl_enclosure()){
                    validated = tv.validate_tail(true);
                }
            }
        }
    }

    const int max_refinement = 1;
    for(int i = 0; validated && i < max_refinement; ++i){
        incl_enclosure();
        tv.validate_tail(true);
    }
    return validated;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class SetType, class ParamType, class DissipativeOperator, class Norm>
bool DissipativeIntegrator<m, M, Scalar, Exponent, VectorType, SetType, ParamType, DissipativeOperator, Norm>::take_step(){
    bool validated = enclosure_with_tail();
    if(!validated) {
        tv.what_failed_in_validation();
        return validated;
    }
    auto &solver = tv.get_diff_inclusion_cw_solver();
    auto mm = solver.getVectorField();
    auto vf = mm.getVectorField();
    auto perturb = mm.getPerturbations();

    x.move(solver);

    // preparing for the new step
    auto x_vector = static_cast<VectorType>(x);
    auto W2 = x_vector + Scalar(0, tv.get_h()) * tv.get_diss_op().P(x_vector);
    auto Th = tv.get_Th();
    tv.reset(Th, x_vector, W2);
    return validated;
}
