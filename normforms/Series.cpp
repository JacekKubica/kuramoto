#include "Series.h"

Series& Series::operator*=(const Series &that) {
    if(this->partials.size() == 0){
        return *this;
    }
    if(that.partials.size() == 0){
        *this = that;
        return *this;
    }
    if(partials.size() != that.partials.size()){
        throw std::logic_error("wrong number of variables series *");
    }
    for(int i = 0; i < partials.size(); ++i){
        partials[i] = that.partials[i] * value + partials[i] * that.value;
    }
    hsum = value * that.hsum + that.value * hsum;
    value = value * that.value;
    return *this;
}

Series& Series::operator+=(const Series &that){
    if(that.partials.size() == 0){
        return *this;
    }
    if(this->partials.size() == 0){
        *this = that;
        return *this;
    }
    value = value + that.value;
    hsum = hsum + that.hsum;
    for(int i = 0; i < partials.size(); ++i) partials[i] = partials[i] + that.partials[i];
    return *this;
}

void Series::change_variables(int k, const Equation &t) {
    hsum = hsum.apply(k, t) - partials[k].apply(k, t) + partials[k].apply(k, t) * t.hsum();
    for(int i = 0; i < partials.size(); ++i) {
        if(i == k) partials[k] = partials[k].apply(k, t) * t.partial(k);
        else partials[i] = partials[i].apply(k, t) + partials[k].apply(k, t) * t.partial(i);
    }
    value = value.apply(k, t);
}

Bound Series::eval(const Set &s) const {
    return value(s);
}

Bound Series::eval_hsum(const Set &s) const {
    return hsum(s);
}

Series Series::make_geometric(const Equation &e, int n){
    Equation enm1 = e;
    // TODO fast exp
    for(int i = 0; i < n - 2; ++i){
        enm1 = enm1 * e;
    }
    auto g_prime_bound = 4 * (n - 1) * enm1 + 4 * n * enm1 * e;
    auto g_val = 2 * enm1 * e;

    std::vector<Equation> partials;
    for(int i = 0; i < e.get_no_of_variables(); ++i) {
        partials.push_back(e.partial(i) * g_prime_bound);
    }
    return Series(2 * e * e, e.hsum() * g_prime_bound, partials);
}
