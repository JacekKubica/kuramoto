#include "capd/diffIncl/InclRect2Set.h"
#include "typedefs.h"
#include "polybd_enclosure.h"

typedef capd::diffIncl::InclRect2Set<capd::IMatrix> Set;

struct C0DubletonSetPolyBds;

struct C0DubletonSetPolyBds {
    Set set;
    PolyBdConstants pbc;  
};

void move(PolyBdSolver&, Set&, Set&); // TODO stub
void move(PolyBdSolver &solver, C0DubletonSetPolyBds &polyset, C0DubletonSetPolyBds &polyset_out) {
    Vector finitePart = Vector(polyset.set);
    PolyBd enclosure = enclose(solver.getStep(),
                               PolyBd{finitePart, polyset.pbc.s, polyset.pbc.C},
                               solver.polybd_vf,
                               polyset_out.pbc);
    move(solver, polyset.set, polyset_out.set); 
}