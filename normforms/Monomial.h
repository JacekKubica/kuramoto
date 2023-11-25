#ifndef _MONOMIAL_H_
#define _MONOMIAL_H_
#include <vector>
#include <initializer_list>
#include <map>
#include <iostream>
#include <numeric>
#include "capd/capdlib.h"
#include "utils.h"


// represents x0^power[0] x1^power[1] ...
class Monomial {
public:
    Monomial(std::vector<int> powers) : powers(powers) {}
    Monomial(std::initializer_list<int> powers) : powers(powers) {}

    friend std::ostream& operator<<(std::ostream &o, const Monomial &m);

    //change of form xi |-> c * xk^ak xl^al...
    Monomial apply(int variable, const Monomial &change) const;

    // so it can be a key in map
    friend bool operator<(const Monomial &m1, const Monomial &m2);
    friend bool operator==(const Monomial &m1, const Monomial &m2){
        return !((m1 < m2) || (m2 < m1));
    }

    friend Monomial operator*(const Monomial &m1, const Monomial &m2);

    int operator[](int index) const {
        return powers[index];
    }
    int& operator[](int index) {
        return powers[index];
    }
    int no_of_variables() const {
        return powers.size();
    }
private:
    std::vector<int> powers;
    friend bool smaller_degree(const Monomial &m1, const Monomial &m2);
    friend bool less_variables(const Monomial &m1, const Monomial &m2);
};

#endif
