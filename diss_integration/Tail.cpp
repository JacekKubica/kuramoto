#include "Tail.h"


template <int m, int M, class Scalar, class Exponent>
bool Tail<m, M, Scalar, Exponent>::check_far_tail_inclusion(const Tail<m, M, Scalar, Exponent> &rhs) const{
    if(rhs.C == 0 && this->C != 0)
        return false;
    else if(this->C != 0){
        if(this->s < rhs.s || !(rhs.C / (Scalar(M + 1) ^ rhs.s) > rhs.C / (Scalar(M + 1) ^ this->s)))
            return false;
    }
    return true;
}

