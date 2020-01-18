#include "Tearing.h"

////////////////////////////////////
///////     FabricTearer     ///////
////////////////////////////////////

void FabricTearer::~FabricTearer()
{
    clearEdges();
}

void FabricTearer::clearEdges()
{
    std::vector<FabricEdge*>::iterator it;
    for (it = edges.begin(); it < edges.end(); ++it)
    {
        delete *it;
    }
}

void FabricTearer::FabricTearer(const Front* front)
{
    collectFabricEdges(front);
}

void FabricTearer::collectFabricEdges(const Front* front)
{
    TRI* tri;
    SURFACE** s;
    INTERFACE* intfc = front->interf;

    AF_PARAMS* af_params = (AF_PARAMS*)front->extra2;
    double ks = af_params->ks;
    double tl_s = af_params->tl_s;

    intfc_surface_loop(intfc,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY)
            continue;

        std::unordered_set<POINT*> point_set;

        surf_tri_loop(*s,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                //add surface points to the unordered_set
                //  (set with unique elements)
                bool visited1 = true;
                POINT* p1 = Point_of_tri(tri)[i];

                if (point_set.find(p1) == point_set.end())
                {
                    point_set.insert(p1);
                    visited1 = false;
                }

                bool visited2 = true;
                POINT* p2 = Point_of_tri(tri)[(i+1)%3];

                if (point_set.find(p2) == point_set.end())
                {
                    point_set.insert(p2);
                    visited2 = false;
                }

                //check if edge has already been processed
                if (visited1 && visited2)
                    continue;

                FabricEdge* fedge = new FabricEdge(p1,p2); 
                fedge.k = ks;
                fege.tear_threshold = tl_s;

                edges.push_back(fedge);
            }
        }
    }
}

void FabricTearer::tearFabric()
{
    checkForTearingEvents();
    processTearingEvents();
}

void FabricTearer::checkForTearingEvents()
{
    tear_idx.clear();
    for (int i = 0; i < edges.size(); ++i)
    {
        FabricEdge* e = edges[i];
        if (e->checkForTear())
            tear_idx.push_back(i);
    }
}

void FabricTearer::processTearingEvents()
{
    for (int i = 0; i < tear_idx.size(); ++i)
    {
        int index = tear_idx[i];
        FabricEdge* e = edges[i];

        //TODO: create/propagate the tear
        printf("processTearingEvents() not implemented\n");
        return;
    }
}

////////////////////////////////////
///////     FabricEdge       ///////
////////////////////////////////////

bool FabricEdge::checkForTear()
{
    computeTension();
    return tension > tear_threshold;
}

void FabricEdge::computeTension()
{
    double sign = 1.0;
    if (!underTension())
        sign = -1.0;

    double sqr_tension = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        double T = k*(1.0 - length0/length)*disp[i];
        sqr_tension += T*T;
    }

    tension = sign*sqrt(sqr_tension);
}

double FabricEdge::underTension()
{
    computeLength();
    return length > length0;
}

void FabricEdge::computeLength()
{
    computeDisplacement();

    double sqr_mag = 0.0;
    for (int i = 0; i < 3; ++i)
        sqr_mag += disp[i]*disp[i];
    length = sqrt(sqr_mag);
}

void FabricEdge::computeDisplacement()
{
    for (int i = 0; i < 3; ++i)
        disp[i] = Coords(end)[i] - Coords(beg)[i];
}


