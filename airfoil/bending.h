#ifndef BENDING_H
#define BENDING_H

#include "airfoil.h"

#include <iomanip>
#include <set>

struct BOND_BENDER
{
    double bends;
    double kb[3];
    double gradprev_kb[3][3];
    double grad_kb[3][3];
    double gradnext_kb[3][3];
};


void addStringBenders(Front* front);
void computeStringBendingForce(INTERFACE* intfc);
void computeSurfBendingForce(INTERFACE* intfc, double bends, double bendd);



#endif
