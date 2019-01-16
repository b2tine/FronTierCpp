#include "spring_solver.h"

void SpringSolver::doSolveDirect(double t) {
    const size_t size = pts.size();
    
    printf("Starting origami transformation:\n");
    for (int i = 0; i < size; i++) 
         computeAccel(pts[i]);
}
