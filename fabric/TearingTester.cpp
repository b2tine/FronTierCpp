#include "Tearing.h"

void FabricTearer::printEdgeTensions()
{
    for (int i = 0; i < edges.size(); ++i)
    {
        printf("edges[%d]->tension = %g\n",
                i, edges[i]->tension);
    }
}

void FabricTearer::setEdgeTension(int index, double Tension)
{
    edges[index]->tension = Tension;
}
