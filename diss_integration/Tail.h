#pragma once
#include "capd/capdlib.h"


// maybe Exponent is not really needed, int should be enough
template <int m, int M, typename Scalar, typename Exponent>
class Tail{
public:

private:
    Scalar C;
    Exponent s;
    std::array<M - m> near_tail;

    bool check_far_tail_inclusion(const Tail&) const;
};


#include "Tail.cpp"
