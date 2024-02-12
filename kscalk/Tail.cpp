#include "Tail.h"


bool Tail::check_far_tail_inclusion(const Tail &rhs){
    if(rhs.C == 0 && this->C != 0)
        return false;
    else if(this->C != 0){
        if(this->s < rhs.s || !(rhs.C / (Scalar(M + 1) ^ rhs.s) > rhs.C / (Scalar(M + 1) ^ this->s)))
            return false;
    }
    return true;
}
