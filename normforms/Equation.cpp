#include "Equation.h"

Equation Equation::apply(int variable, const Equation &change) const {
    if(this->no_of_variables == 0){
        return *this;
    }
    if(this->no_of_variables != change.no_of_variables){
        COUT(this->no_of_variables);
        COUT(*this);
        COUT(change.no_of_variables);
        COUT(change);
        throw std::logic_error("wrong number of variables");
    }
    Equation ret;
    for(auto term: terms){
        
        int p = term.first[variable];
        if(p != 0){
            auto monomial_without_variable = term.first;
            monomial_without_variable[variable] = 0;
            ret = ret + change.take_power(p) * Equation(monomial_without_variable, term.second);  
        }
        else {
            ret = ret + Equation(term.first, term.second);
        }
    }
    return ret;
}


// function used to calculate bar{m}
// Equation Equation::apply_rest(int variable, const Equation &change) const {
//     return this->apply(variable, change) - *this;
// }

Equation Equation::nonlinear_part() const {
    auto ret = *this;
    for(auto it = ret.begin(); it != ret.end();){
        int count = 0;
        for(int i = 0; i < it->first.no_of_variables(); ++i){
            count += it->first[i];
        }
        if(count == 1){
            ret.terms.erase(it++);
        } else{
            ++it;
        }
    }
    return ret;
}

Equation operator+(const Equation &e1, const Equation &e2){
    if(e2.no_of_variables == 0)
        return e1;
    if(e1.no_of_variables == 0)
        return e2;
    if(e1.no_of_variables != e2.no_of_variables){
        throw std::logic_error("wrong number of variables eq +");
    }
    if(e1.no_of_variables < e2.no_of_variables)
        return e2 + e1;
    Equation ret = e1;
    for(auto const &term : e2.terms){
        map_add<Monomial, capd::interval>(ret.terms, term.first, term.second);
    }
    return ret;
}

Equation operator*(const Equation &e1, const Equation &e2){
    if(e2.no_of_variables == 0)
        return e2;
    if(e1.no_of_variables == 0)
        return e1;
    if(e1.no_of_variables != e2.no_of_variables){
        throw std::logic_error("wrong number of variables eq *");
    }
    Equation ret;
    for(auto const &term1 : e1.terms){
        for(auto const &term2 : e2.terms){
            ret = ret + Equation(term1.first * term2.first, term1.second * term2.second);
        }
    }
    return ret;
}

Equation operator*(const capd::interval &coeff, const Equation &e){
    Equation ret = e;
    for(auto &term : ret.terms){
        term.second *= coeff;
    }
    return ret;
}

Equation operator*(const Equation &e, const capd::interval &coeff){
    return coeff * e;
}

Equation Equation::take_power(int p) const {
    Equation ret = *this;
    for(int i = 0; i < p - 1; ++i){
        ret = ret * (*this);
    }
    return ret;
}

std::ostream& operator<<(std::ostream &o, const Equation &e){
    // for(const auto &term: e.terms) o << term << term == e.terms.end() ? " + " : "";
    for(const auto &term: e.terms) o << term.second << " " << term.first << " + ";
    return o;
}

Bound Equation::operator()(const Set &s) const {
    auto ret = Bound(0, 0, s.get_constant());
    if(no_of_variables == 0) return ret;
    for(const auto &term: terms) {
        Bound temp = Bound(0, 0, s.get_constant());
        for(int i = 0; i < term.first.no_of_variables(); ++i) {
            if(term.first[i] != 0){
                if(temp.coeff == 0) temp.coeff = 1;
                temp = temp * power(s[i], term.first[i]);
            }
        }
        // we need to be pessimistic in eval,
        // it is easier to change here than
        // changin transformations (and keeping transformations exact
        // may be useful anyway)
        ret = ret + capd::abs(term.second) * temp;
    }
    return ret;
}

Equation Equation::hsum() const {
    Equation ret;
    for(const auto &term: terms) {
        for(int i = 0; i < term.first.no_of_variables(); ++i){
            if(term.first[i] != 0) {
                auto temp = term.first;
                int pow = temp[i];
                temp[i] -= 1;
                ret = ret + Equation(temp, capd::abs(term.second) * term.first[i]);
            }
        }
    }
    return ret;
}

Equation Equation::partial(int i) const {
    Equation ret;
    for(const auto &term: terms) {
        if(term.first.no_of_variables() > i and term.first[i] != 0) {
            auto temp = term.first;
            int pow = temp[i];
            temp[i] -= 1;
            ret = ret + Equation(temp, term.second * term.first[i]);
        }
    }
    return ret;
}

std::vector<Equation> Equation::partials() const {
    std::vector<Equation> ret;
    if(terms.size() == 0) return ret;
    // TODO this no_of_variables should be stored in Equation and checked when adding monomial
    for(int i = 0; i < no_of_variables; ++i){
        ret.push_back(partial(i));
    }
    return ret;
}
