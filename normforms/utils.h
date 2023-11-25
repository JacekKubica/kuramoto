#ifndef _UTILS_H_
#define _UTILS_H_
#include <map>
#include "capd/capdlib.h"
#include <iostream>



template <typename T, typename Arithmetic>
void map_add(std::map<T, Arithmetic> &m, const T& key, const Arithmetic& value){
    if(m.count(key)){
        m[key] += value;
        return;
    }
    m[key] = value;
}

#endif
