#ifndef _SERIES_H_
#define _SERIES_H_
#include "Equation.h"
#include "capd/capdlib.h"
#include <iostream>


class Series {
public:
    Series(int no_of_variables) : value(Equation()), hsum(Equation()), partials(std::vector<Equation>(no_of_variables)), no_of_variables(no_of_variables) {}

    Series(const Equation &value, const Equation& hsum,
           const std::vector<Equation> &partials)
        : value(value),
          hsum(hsum),
          partials(partials),
          no_of_variables(partials.size()) {}

    Series(const Equation &equation)
        : value(equation),
          hsum(equation.hsum()),
          partials(equation.partials()),
          no_of_variables(equation.get_no_of_variables()) {
          }

    static Series make_geometric(const Equation &e, int n=2);

    Series& operator*=(const Series &that);
    friend Series operator*(Series s1, const Series &s2) {
        s1 *= s2;
        return s1;
    }
    friend Series operator*(capd::interval a, Series s) {
        s.value = s.value * a;
        s.hsum = s.hsum * a;
        for(int i = 0; i < s.partials.size(); ++i) s.partials[i] = s.partials[i] * i;
        return s;
    }
    
    Series& operator+=(const Series &that);

    friend Series operator+(Series s1, const Series &s2) {
        s1 += s2;
        return s1;
    }

    friend Series operator-(const Series &s1, const Series &s2) {
        return s1 + capd::interval(-1) * s2;
    }
    void change_variables(int k, const Equation &t);

    // copy set
    Bound eval(const Set &s) const;
    Bound eval_hsum(const Set &s) const;
    Bound eval_partial(int i, const Set &s) const {
        return partials[i](s);
    }
    

    Equation get_value() const {
        return value;
    }
    
    Equation get_hsum() const {
        return hsum;
    }

    friend std::ostream& operator<<(std::ostream &o, const Series &s){
        o << "val fun: " << s.value; 
        o << "hsum fun: " << s.hsum;
        for(int i = 0; i < s.partials.size(); ++i) {
            o << "partial" << i << " " << s.partials[i];
        }
        return o;
    }

    void print(const Set &s) const {
        std::cout << "on set " << s << '\n';
        std::cout << "val fun: " << value << " eval: " << value(s) << '\n';
        std::cout << "hsum fun: " << hsum << " eval: " << hsum(s) << '\n';
        for(int i = 0; i < partials.size(); ++i) {
            std::cout << "partial" << i << " " << partials[i] << " eval: " << partials[i](s) << '\n';
        }
    }

    Equation get_partial(int i) const {
        return partials[i];
    }
private:
    std::vector<Equation> partials;
    Equation value, hsum;
    int no_of_variables;
};

#endif
