#include "Set.h"

Bound operator+(const Bound &b1, const Bound &b2){
    if(b1.coeff == 0) return b2; if(b2.coeff == 0) return b1;

    if(b1.c != b2.c) throw std::logic_error("Adding bounds with different constants!");
    if(b1.degree > b2.degree) return b2 + b1;
    return Bound(b1.coeff + b2.coeff * power(b1.c->val, b2.degree - b1.degree), b1.degree, b1.c);
}

Bound operator*(const capd::interval &i, const Bound &b){
    auto ret = b;
    ret.coeff *= i;
    return ret;
}

Bound operator-(const Bound &b1, const Bound &b2){
    return b1 + (-1) * b2;
}

Bound operator*(const Bound &b1, const Bound &b2){
    if(b1.c != b2.c) throw std::logic_error("Multiplying bounds with different constants!");
    return Bound(b1.coeff * b2.coeff, b1.degree + b2.degree, b1.c);
}

Bound power(const Bound &b, int p) {
    return Bound(power(b.coeff, p), p * b.degree, b.c);
}
