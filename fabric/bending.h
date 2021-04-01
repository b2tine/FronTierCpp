#ifndef BENDING_H
#define BENDING_H

#include "fabric.h"

#include <iomanip>
#include <set>

struct BOND_BENDER
{
    double kb[3];
    double gradprev_kb[3];
    double grad_kb[3];
    double gradnext_kb[3];
};


void addStringBenders(INTERFACE* intfc);

void computeBendingForce(INTERFACE* intfc, double bends, double bendd);



#endif
