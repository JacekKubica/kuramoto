#ifndef _EQUATIONS_H_
#define _EQUATIONS_H_
#include <vector>
#include <initializer_list>
#include <iostream>
#include "Equation.h"
#include "Set.h"
#include "Series.h"
#include <algorithm>


class Equations {
public:
    auto begin() { return equations.begin(); }
    auto end() { return equations.end(); }
    auto cbegin() const { return equations.begin(); }
    auto cend() const { return equations.end(); }
    auto begin() const { return equations.begin(); }
    auto end() const { return equations.end(); }
    const std::vector<Series> get_series() const {
        return series;
    }

    Equations(std::vector<Equation> equations, int no_of_variables)
        : equations(equations), 
          series(std::vector<Series>(no_of_variables, {no_of_variables})),
          vsums_bounds(std::vector<Equation>(no_of_variables)), 
          no_of_variables(no_of_variables) {}

    Equations(const std::vector<Equation> &equations, const std::vector<Series> &series, const std::vector<Equation> &vsums_bounds, int no_of_variables)
        : equations(equations), series(series), vsums_bounds(vsums_bounds), no_of_variables(no_of_variables) {
    }

    Equations apply(int variable, const Equation &change, bool remove=false, int gn=2) const;

    Equations nonlinear_part() const{
        auto ret = *this;
        for(int i = 0; i < no_of_variables; ++i){
            ret[i] = ret[i].nonlinear_part();
        }
        return ret;
    }

    int get_no_of_variables(){
        return no_of_variables;
    }

    const Equation& operator[](int index) const {
        return equations[index];
    }
    Equation& operator[](int index) {
        return equations[index];
    }

    std::vector<Series> get_all_series() {
        return series;
    }

    friend std::ostream& operator<<(std::ostream &o, const Equations &eqs){
        int i = 0;
        for(auto const &e: eqs.equations){
            o << e << "\n"; 
            // o << ", series: " << eqs.series[i] << "\n";
            ++i;
            o << "\n";
        }
        return o;
    }

    bool check_geom_args(capd::interval max, const Set &s){
        for(auto &ga: geom_args){
            COUT(ga(s).eval());
            if(ga(s).eval() > max) return false;
        }
        return true;
    }

    std::vector<Bound> vals(const Set &s) const{
        std::vector<int> v(no_of_variables); std::iota(v.begin(), v.end(), 0); std::vector<Bound> ret; ret.reserve(no_of_variables); 
        std::transform(v.begin(), v.end(), ret.begin(),
                       [&s, this](int i){ return (*this)[i](s) + series[i].eval(s); });

        return ret;
    }

    Bound hsum(int i, const Set &s) const  {
        return equations[i].hsum()(s) + series[i].eval_hsum(s);
    }

    std::vector<Bound> hsums(const Set &s) const{
        std::vector<int> v(no_of_variables); std::iota(v.begin(), v.end(), 0); std::vector<Bound> ret; ret.reserve(no_of_variables);
        std::transform(v.begin(), v.end(), ret.begin(),
                       [&s, this](int i){ return this->hsum(i, s); });
        return ret;
    }


    Bound vsum(int k, const Set &s) const {
        Bound b(0, 0, s.get_constant());
        for(int i = 0 ; i < no_of_variables; ++i) {
            b = b + equations[i].partial(k)(s) + series[i].eval_partial(k, s);
        }
        b = b + vsums_bounds[k](s);
        return b;
    }

    std::vector<Bound> vsums(const Set &s) const{
        std::vector<int> v(no_of_variables); std::iota(v.begin(), v.end(), 0); std::vector<Bound> ret; ret.reserve(no_of_variables);
        std::transform(v.begin(), v.end(), ret.begin(),
                       [&s, this](int i){ return this->vsum(i, s); });
        return ret;
    }

    std::vector<Bound> deriv_bounds(const Set &s) const{
        std::vector<int> v(no_of_variables); std::iota(v.begin(), v.end(), 0); std::vector<Bound> ret; ret.reserve(no_of_variables);
        std::transform(v.begin(), v.end(), ret.begin(),
                       [&s, this](int i){ return this->hsum(i, s) + this->vsum(i, s); });
        return ret;
    }


    
private:
    std::vector<Equation> equations;
    std::vector<Equation> geom_args;
    std::vector<Equation> vsums_bounds;
    std::vector<Series> series;
    int no_of_variables;
    Equation value_after;
};


#endif
