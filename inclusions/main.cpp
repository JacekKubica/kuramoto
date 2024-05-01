#include "polybd_enclosure.h"

using namespace std;
using namespace capd;
#include <iostream>

int main() {
    cout.precision(12);
    IMap vectorField("var:x,y,z;fun:x,y,-z+z*z;");
    IMap nonlinear("var:x, y, z;fun:0,0,z*z;");
    //vectorField.setParameter("a",interval(57.)/10.); // a=5.7
    //vectorField.setParameter("b",interval(2.)/10.); // b=0.2
    IVector c(3);
    c[0] = 0.;
    c[1] = 0.;
    c[2] = 1;
    c += interval(-1,1)*1e-1; 
    VectorField vf {vectorField, nonlinear, {1, 1, 1}};
    auto Z = enclose(0.01, c, vf);
    COUT(Z);
    vf.linear[2] = -1;
    Z = enclose(0.01, c, vf);
    COUT(Z);
}
// TODO g definition error in my summary!