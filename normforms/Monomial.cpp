#include "Monomial.h"
#include <exception>
#include <algorithm>

bool smaller_degree(const Monomial &m1, const Monomial &m2){
    auto lambda = [](int x, int y) { return x + y; };
    int s1 = std::accumulate(m1.powers.begin(), m1.powers.end(), 0, lambda);
    int s2 = std::accumulate(m2.powers.begin(), m2.powers.end(), 0, lambda);
    if(s1 < s2)
        return true;
    return false;
}

bool less_variables(const Monomial &m1, const Monomial &m2){
    auto lambda = [](int x, int y) { return x + int(y > 0); };
    int s1 = std::accumulate(m1.powers.begin(), m1.powers.end(), 0, lambda);
    int s2 = std::accumulate(m2.powers.begin(), m2.powers.end(), 0, lambda);
    if(s1 < s2)
        return true;
    return false;
}

bool operator<(const Monomial &m1, const Monomial &m2){
    if(smaller_degree(m1, m2)) return true;
    if(smaller_degree(m2, m1)) return false;
    if(less_variables(m1, m2)) return true;
    if(less_variables(m2, m1)) return false;

    for(int i = 0; i < m1.no_of_variables(); ++i){
        if(m1[i] == 0 && m2[i] != 0) return false;
        if(m2[i] == 0 && m1[i] != 0) return true;
        if(m1[i] < m2[i]) return true;
        if(m1[i] > m2[i]) return false;
    }
    return false;
}

Monomial operator*(const Monomial &m1, const Monomial &m2){
    if(m1.no_of_variables() != m2.no_of_variables()){
        throw std::logic_error("wrong number of variables");
    }
    if(m1.no_of_variables() < m2.no_of_variables())
        return m2 * m1;
    auto ret = m1;
    for(int i = 0; i < m2.no_of_variables(); ++i){
        ret[i] += m2[i];
    }
    return ret;
}

std::ostream& operator<<(std::ostream &o, const Monomial &m){
    bool empty = true;
    for(int i = 0; i < m.no_of_variables(); ++i){
        if(m.powers[i] != 0){ o << "x_" << i << "^" << m.powers[i] << " "; empty = false; }
    }
    if(empty) o << "1";
    return o;
}


//change of form xi |-> xk^ak xl^al...f_1^b_1...f_n^b_n
Monomial Monomial::apply(int variable, const Monomial &change) const {
    if(this->no_of_variables() != change.no_of_variables())
        throw std::logic_error("change of variables needs to describe all variables");
    Monomial temp = *this;
    int power = temp[variable];
    temp[variable] = 0;
    for(int i = 0; i < change.no_of_variables(); ++i){
        if(change[i] != 0){
            temp[i] += power * change[i];
        }
    }
    return temp;
}
