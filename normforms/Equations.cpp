#include "Equations.h"

Equations Equations::apply(int variable, const Equation &change, bool remove, int gn) const {
    auto ret = *this;
    // ret[variable] = Equation();
    for(auto &equation: ret) {
        equation = equation.apply(variable, change);
    }
    for(auto &equation: ret.geom_args) {
        equation = equation.apply(variable, change);
    }
    for(auto &equation: ret.vsums_bounds) {
        equation = equation.apply(variable, change);
    }
    for(auto &one_series: ret.series) {
        one_series.change_variables(variable, change);
    }

    for(int i = 0; i < change.get_no_of_variables(); ++i) {
        if(i == variable) continue;
        auto temp = change.partial(i);
        ret[variable] = ret[variable] - temp * ret.equations[i]; // here ret because we need to change variables in f_i first THAT WAS WRONG ->  should be original equations used here 
        ret.series[variable] = ret.series[variable] + Series(temp) * ret.series[i]; // was minus here, probably wrong, because we always bound series from the above
    }

    // std::vector<int> one(no_of_variables);
    // one[0] = 1;
    // auto e_one = Equation(Monomial(one));
    auto without_linear = change.nonlinear_part();
    auto partial_variable = without_linear.partial(variable);
    auto geom = Series::make_geometric(partial_variable, gn);
    ret.geom_args.push_back(partial_variable);
    // ??? ret.series[variable] = ret.series[variable] + ret.series[variable] * Series(-1 * partial_variable) + geom;
    ret.series[variable] = ret.series[variable] + ret.series[variable] * Series(partial_variable) + ret.series[variable] * geom;
    ret.series[variable] = ret.series[variable] + geom * ret[variable];
    
    auto temp = partial_variable;
    for(int i = 0; i < gn - 1; ++i){
        ret[variable] = ret[variable] - ret[variable] * temp;
        temp = -1 * partial_variable * temp;
    }


    if(remove){
        for(auto const &term: without_linear){
            ret[variable].remove(term.first);
        }
    }
    else {
        for(auto const &term: without_linear){
            COUT(ret[variable][term.first]);
        }
    }

    return ret;
}
