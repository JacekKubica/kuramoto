#ifndef _EQUATION_H_
#define _EQUATION_H_
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include "capd/capdlib.h"
#include "utils.h"
#include "Monomial.h"
#include <exception>
#include "Set.h"


class Equation {
public:
    Equation() {}
    Equation(const Monomial &m, capd::interval coeff=1) : no_of_variables(m.no_of_variables()) { terms[m] = coeff; }
    Equation(std::initializer_list<int> l, capd::interval coeff=1) : no_of_variables(l.size()) { terms[Monomial(l)] = coeff; }

    Equation apply(int variable, const Equation &change) const;
    Equation apply_rest(int variable, const Equation &change) const;
    Equation nonlinear_part() const;

    friend bool operator==(const Equation &e1, const Equation &e2){
        return e1.terms == e2.terms && e1.no_of_variables == e2.no_of_variables;
    }

    friend Equation operator+(const Equation &e1, const Equation &e2);
    friend Equation operator*(const Equation &e1, const Equation &e2);
    friend Equation operator*(const capd::interval &coeff, const Equation &e);
    Equation take_power(int p) const;
    friend Equation operator-(const Equation &e1, const Equation &e2){
        return e1 + (-1) * e2;
    }
    friend Equation operator*(const Equation &e, const capd::interval &coeff);
    friend std::ostream& operator<<(std::ostream &o, const Equation &e);

    // TODO 
    Equation operator+=(const Equation &e) const { throw std::logic_error("+= not implemented yet"); return *this + e; }

    Bound operator()(const Set &s) const;
    Equation hsum() const;
    Equation partial(int i) const;
    std::vector<Equation> partials() const;
    int get_no_of_variables() const {
        return no_of_variables;
    }

    void print_lowest_order_terms(const Set &s) const {
        int d = this->operator()(s).degree;
        for(const auto &term: terms) {
            Equation temp(term.first);
            if(temp(s).degree == d){
                COUT(temp);
                COUT(temp(s))
            }
        }
    }

    const std::map<Monomial, capd::interval>& get_terms() const{
        return terms;
    }

    
    auto begin() { return terms.begin(); }
    auto end() { return terms.end(); }
    auto cbegin() const { return terms.begin(); }
    auto cend() const { return terms.end(); }
    auto begin() const { return terms.begin(); }
    auto end() const { return terms.end(); }
    capd::interval get_coeff(const Monomial &m) const { return terms.at(m); }
    void set_coeff(const Monomial &m, capd::interval c) { terms.at(m) = c; }
    void remove(const Monomial &m) { terms.erase(m); }
    const capd::interval& operator[](const Monomial &m) const {
        return terms.at(m);
    }
    
private:
    std::map<Monomial, capd::interval> terms;
    int no_of_variables = 0;
};


#endif
