#include "Tearing.h"

void FabricTearer::collectFabricEdges(const INTERFACE* intfc)
{
    TRI* tri;
    SURFACE** s;

    intfc_surface_loop(intfc,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY)
            continue;

        surf_tri_loop(*s,tri)
        {
            //TODO:
            //visit edges of triangle, and push onto the
            //edges vector making sure not to double count.
        }
    }
}
