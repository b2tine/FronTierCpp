#include "cgal_intfc.h"
#include "airfoil.h"



extern void CGAL_MakeDiskGapBandSurf(
        Front* front,
        double* center,
        double radius,
        double height,
        int idir,
        COMPONENT neg_comp,
        COMPONENT pos_comp,
        int w_type,
        int refinement_level,
        SURFACE** surf)
{
    CGAL_MakeCylindricalShellSurf(front,center,radius,height,idir,
            neg_comp,pos_comp,w_type,refinement_level,surf);
}
