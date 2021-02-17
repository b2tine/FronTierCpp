#ifndef IFTURB_H
#define IFTURB_H

#include "iFluid.h"

class TurbSolver
{
    public:

        TurbSolver(Front* ft);
        virtual ~TurbSolver() = default;

        virtual void initMesh();		    // set up the cartesian grid
        virtual void initFieldVariables();  // allocate field variables memory (previously setDomain())
        virtual void setComponent();	    // initialize cell components	
   
   
    protected:

        Front* front;
        
        // Topological grid variables (Dual to Comp Grid)
        RECT_GRID* top_grid;
        COMPONENT* top_comp;
        int* top_gmax;
        double* top_L;
        double* top_U;
        int* lbuf;
        int* ubuf;

	    std::vector<IF_RECTANGLE> cells; //Centered at nodes of top_grid

        // Parallel partitioning
        int NLblocks;
        int ilower;
        int iupper;
        int* n_dist;
    
        double dt; //time step
};


class TurbSolver2d : public TurbSolver
{
    protected:

        // Parallel solver index mapping
        int** ij_to_I;
        int** I_to_ij;

};



class TurbSolver3d : public TurbSolver
{
    protected:

        // Parallel solver index mapping
        int*** ijk_to_I;
        int** I_to_ijk;

};




#endif
