
template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
Scalar TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::eval_g(int i) const {
    return Scalar(((T0[i].leftBound() - b[i].leftBound()) * exp(diss_op.lambda(i) * h) + b[i].leftBound()).leftBound(),
                  ((T0[i].rightBound() - b[i].rightBound()) * exp(diss_op.lambda(i) * h) + b[i].rightBound()).rightBound());
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
void TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::set_N_b_and_g(){
    N = diss_op.N(W2);
    for (int i = m; i < M; ++i) {
        b[i] = N[i] / (-diss_op.lambda(i));
    }
    b.C() = N.C() / diss_op.V();
    b.exponent() = N.exponent() + diss_op.p();
}


// TODO T in diss_op should be somehow made const
template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::TailValidator(const typename Scalar::BoundType h, const VectorType &x, const capd::pdes::PolynomialBound<Scalar, Exponent, M> &T0, const ParamType &p) :
    h(h), T0(T0), diss_op{p, T}, W2(x + Scalar(0, h) * diss_op.P(x)),
    validated(false), far_tail_validated(false), kvalidated(std::vector<bool>(M, false)),
    diff_inclusion_cw_solver{diss_op.perturbated(), 4, capd::IMaxNorm(), diss_op} { //TODO norm generic

    // just dummy set to set step
    capd::C0Rect2Set x_set(x);
    diff_inclusion_cw_solver.setStep(h);
    x_set.move(diff_inclusion_cw_solver.getDynamicalSystem());
    //end setting step

    T = diss_op.guess_far_tail(x, T0); // important to note - in diss_op constructor only ref to T, which until now is empty, is saved!
    diff_inclusion_cw_solver.setStep(h);
    diss_op.update_perturbation(W2);
    set_N_b_and_g();
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
void TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::inflate(Scalar &X, typename Scalar::BoundType c){
    Scalar temp;
    X.split(temp);
    temp *= c;
    X += temp;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
bool TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::check_far_tail_inclusion(const capd::pdes::PolynomialBound<Scalar, Exponent, M> &pb1, const capd::pdes::PolynomialBound<Scalar, Exponent, M> &pb2){
    if(pb2.C() == 0 && pb1.C() != 0)
        return false;
    else if(pb1.C() != 0){
        if(pb1.exponent() < pb2.exponent() || pb1[M].rightBound() > pb2[M].rightBound())
            return false;
    }
    return true;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
bool TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::validate_far_tail(bool update) {
    if(!check_far_tail_inclusion(T0, T)) { return false; }
    if(b.exponent() > T0.exponent() && T0.C() != 0){
        int L = static_cast<int>(power(b.C() / T0.C(), Scalar(1) / (b.exponent() - T0.exponent())).rightBound());
        for(int i = M; i <= L; ++i){
            if(T[i].rightBound() < eval_g(i).rightBound()){
                if(update){
                    update_tail(true);
                }
                return false;
            }
        }
        return true;
    }
    else{
        // std::cout << "*******************" << '\n';
        // std::cout << T << '\n';
        // std::cout << "*******************" << '\n';
        // std::cout << T0 << '\n';
        // std::cout << "===================" << '\n';
        // std::cout << b << '\n';
        if(b.exponent() == T0.exponent() || T0.C() == 0){
            if(T0.C() >= b.C())
                return true;
            else if(T.exponent() <= b.exponent() && T[M].rightBound() >= b[M].rightBound()){
                return true;
            }
        }
        if(b.exponent() < T0.exponent() && T0.C() != 0){
            int L = static_cast<int>(power(T0.C() / b.C(), 1 / (T0.exponent() - b.exponent())).rightBound());
            auto p = std::max(L, M);
            if(T.exponent() <= b.exponent() && T[p].rightBound() >= b[p].rightBound())
                return true;
        }
        if(update){
            update_tail(true);
            return validate_far_tail();
        }
    }
    return false;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
void TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::update_tail(bool only_far_tail){

    if(far_tail_validated && check_far_tail_inclusion(b, T0) && T0.C() != 0){
        T.C() = T0.C() * take_power(M + 1, T.exponent() - T0.exponent());
    }
    if(!only_far_tail && validated){
        for(int i = m; i < M; ++i){
            if(b[i].rightBound() <= T0[i].rightBound())
                T[i].setRightBound(T0[i].rightBound());
            else if(b[i].rightBound() > T0[i].rightBound())
                T[i].setRightBound(eval_g(i).rightBound());
            if(b[i].leftBound() >= T0[i].leftBound())
                T[i].setLeftBound(T0[i].leftBound());
            else if(b[i].leftBound() < T0[i].leftBound())
                T[i].setLeftBound(eval_g(i).leftBound());
        }
    }
    else{
        if(!only_far_tail){
            for(int i = m; i < M; ++i){
                if(!kvalidated[i]){
                    if(b[i].rightBound() > T[i].rightBound())
                        T[i].setRightBound((1 - DG) * eval_g(i).rightBound() + DG * b[i].rightBound());
                    if(b[i].leftBound() < T[i].leftBound())
                        T[i].setLeftBound((1 - DG) * eval_g(i).leftBound() + DG * b[i].leftBound());
                }
                inflate(T[i], D2);
            }
        }
        if(!far_tail_validated){
            T.C() = capd::max(D2 * b.C() * take_power(M + 1, T.exponent() - b.exponent()),
                              T0.C() * take_power(M + 1, T.exponent() - T0.exponent()));
        }
    }

    // since T changed, we need to recompute those
    diss_op.update_perturbation(W2);
    set_N_b_and_g();
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
bool TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::validate_tail(bool gen_new){
    far_tail_validated = validate_far_tail();

    for(int i = m; i < M; ++i){
        if(T0[i].rightBound() < b[i].rightBound() && T[i].rightBound() < eval_g(i).rightBound()){
            kvalidated[i] = false;
        }
        else if(T0[i].leftBound() > b[i].leftBound() && T[i].leftBound() > eval_g(i).leftBound()){
            kvalidated[i] = false;
        }
        else{
            kvalidated[i] = true;
        }
    }

    validated = far_tail_validated && std::all_of(kvalidated.begin() + m, kvalidated.end(), [](bool v){return v;});

    if(gen_new) update_tail();
    return validated;
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
void TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::reset(const capd::pdes::PolynomialBound<Scalar, Exponent, M> &new_T0, const VectorType &x, const VectorType &new_W2) {
    validated = false; far_tail_validated = false; kvalidated = std::vector<bool>(M, false);
    T0 = new_T0;
    T = diss_op.guess_far_tail(x, T0);
    W2 = new_W2;
    diss_op.update_perturbation(W2);
    set_N_b_and_g();
}

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
capd::pdes::PolynomialBound<Scalar, Exponent, M> TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::get_Th() const {
    // TODO compute intersection of PolynomialBounds here
    diss_op.verify_step(h, b.exponent(), T0.exponent());
    capd::pdes::PolynomialBound<Scalar, Exponent, M> Th(M, T0.C() * diss_op.E(b.exponent() - T0.exponent(), h) + b.C(), b.exponent());
    Th.C() = capd::min(Th.C(), T.C() * take_power(M + 1, Th.exponent() - T.exponent())); // should it be here?
    for(int i = m; i < M; ++i) Th[i] = eval_g(i);
    for(int i = m; i < M; ++i){
        intersection(Th[i], T[i], Th[i]);
    }
    return Th;
}


// TESTING FUNCTIONS

template<int m, int M, class Scalar, class Exponent, class VectorType, class ParamType, class DissipativeOperator, class Norm>
void TailValidator<m, M, Scalar, Exponent, VectorType, ParamType, DissipativeOperator, Norm>::what_failed_in_validation(){
    std::cout << "WHAT FAILED" << '\n';
    std::cout << "with T0: " << T0 << '\n';
    std::cout << "and T: " << T << '\n';
    std::cout << "and b: " << b << '\n';
    std::cout << "and g:" << '\n';
    for (size_t i = m; i < M; i++) {
        std::cout << eval_g(i) << ", ";
    }
    std::cout << '\n';
    if(!check_far_tail_inclusion(T0, T)) { std::cout << "FAR TAIL INCLUSION" << '\n'; }
    if(b.exponent() > T0.exponent() && T0.C() != 0){
        int L = static_cast<int>(power(b.C() / T0.C(), Scalar(1) / (b.exponent() - T0.exponent())).rightBound());
        for(int i = M; i < L; ++i){
            if(T[i].rightBound() < eval_g(i).rightBound())
                std::cout << "T[i].rightBound() < eval_g(i).rightBound()" << '\n';
        }
    }
    else if(b.exponent() == T0.exponent() || T0.C() == 0){
        if(T0.C() >= b.C());
        else if(T.exponent() <= b.exponent() && T[M].rightBound() >= b[M].rightBound());
        else std::cout << "in b.exponent() == T0.exponent() || T0.C() == 0" << '\n';
    }
    else if(b.exponent() < T0.exponent() && T0.C() != 0){
        int L = static_cast<int>(power(T0.C() / b.C(), 1 / (T0.exponent() - b.exponent())).rightBound());
        auto p = std::max(L, M);
        if(T.exponent() <= b.exponent() && T[p].rightBound() >= b[p].rightBound());
        else std::cout << "in b.exponent() < T0.exponent() && T0.C() != 0" << '\n';
    }

    for(int i = m; i < M; ++i){
        if(T0[i].rightBound() < b[i].rightBound() && T[i].rightBound() < eval_g(i).rightBound()){
            std::cout << "for " << i << '\n';
            std::cout << "T0[i].rightBound() < b[i].rightBound() && T[i].rightBound() < eval_g(i).rightBound()" << '\n';
        }
        if(T0[i].leftBound() > b[i].leftBound() && T[i].leftBound() > eval_g(i).leftBound()){
            std::cout << "for " << i << '\n';
            std::cout << "T0[i].leftBound() > b[i].leftBound() && T[i].leftBound() > eval_g(i).leftBound()" << '\n';
        }

    }

}
