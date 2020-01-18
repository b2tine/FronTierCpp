#include "Tearing.h"

////////////////////////////////////
///////     FabricTearer     ///////
////////////////////////////////////

/*
void FabricTearer::setEdgeTension(int index, double Tension)
{
    edges[index]->tension = Tension;
}
*/

void FabricTearer::tearFabricTest()
{
    checkForTearingEventsTest();
    processTearingEvents();
    printEdges();
}

void FabricTearer::checkForTearingEventsTest()
{
    tear_idx.clear();
    for (int i = 0; i < edges.size(); ++i)
    {
        FabricEdge* e = edges[i];
        if (e->checkForTearTest(i))
            tear_idx.push_back(i);
    }
}

////////////////////////////////////
///////     FabricEdge       ///////
////////////////////////////////////

bool FabricEdge::checkForTearTest(int i)
{
    computeTension();
    
    //artificially set tension for testing tearing
    if (i == 1000)
        tension = 5000.0;
    
    return tension > tear_threshold;
}
