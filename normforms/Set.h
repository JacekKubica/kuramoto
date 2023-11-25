#ifndef _SET_H_
#define _SET_H_
#include "capd/capdlib.h"
#include <iostream>
#include <exception>

struct Constant {
    Constant(capd::interval val) : val(val) {}
    capd::interval val;
};

struct Bound {
    Bound(capd::interval coeff, int degree, const Constant *c) : coeff(coeff), degree(degree), c(c) {}
    capd::interval eval() const { return coeff * power(c->val, degree); }
    friend Bound operator+(const Bound &b1, const Bound &b2);
    friend Bound operator*(const capd::interval &i, const Bound &b);
    friend Bound operator-(const Bound &b1, const Bound &b2);
    friend Bound operator*(const Bound &b1, const Bound &b2);
    friend bool operator==(const Bound &b1, const Bound &b2){
        return b1.coeff == b2.coeff && b1.degree == b2.degree && b1.c == b2.c;
    }
    friend Bound power(const Bound &b, int p);
    friend std::ostream& operator<<(std::ostream &o, const Bound &b){
        o << b.coeff << " * " << b.c->val << "^" << b.degree << " = " << b.eval();
        return o;
    }
    const Constant *c;
    capd::interval coeff;
    int degree;
};

class Set {
public:
    Set(std::vector<Bound> values, const Bound &C) : values(values), C(C) {}

    Bound operator[](int i) const {
        if(i < values.size())
            return values[i];
        if(i == values.size())
            return C;
        throw std::logic_error("out of bounds");
    }

    Bound& operator[](int i) {
        if (i >= values.size())
            throw std::logic_error("out of bounds - cannot modify");
        return values[i];
    }

    int main_modes() const {
        return values.size();
    }

    const Constant* get_constant() const {
        return C.c;
    }

    Bound get_C() const {
        return C;
    }

    friend std::ostream& operator<<(std::ostream &o, const Set &s){
        for(auto const v: s.values){
            o << v << " ";
        }
        o << s.C;
        return o;
    }
private:
    std::vector<Bound> values;
    Bound C;
};

#endif
