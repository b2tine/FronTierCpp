#ifndef IMPLICIT_H
#define IMPLICIT_H

#include "weno.h"
#include "solver.h"


void WenoFlux(int mesh_size, double *u_old,
        double *cflux, double dx, double dt);

void implicitSolver(int mesh_size, double *u_old,
        double *u_new, double *source, double dx, double dt);



#endif
