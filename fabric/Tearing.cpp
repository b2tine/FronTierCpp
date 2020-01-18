#include "Tearing.h"
//#include "fabric.h"

////////////////////////////////////
///////     FabricTearer     ///////
////////////////////////////////////

FabricTearer::~FabricTearer()
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

void FabricTearer::collectFabricEdges(const INTERFACE* intfc)
{
    clearEdges();

    TRI* tri;
    SURFACE** s;

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
                edges.push_back(fedge);
            }
        }
    }
}

void FabricTearer::setSpringData(
        double spring_constant,
        double max_tension)
{
    for (int i = 0; i < edges.size(); ++i)
    {
        edges[i]->setSpringConstant(spring_constant);
        edges[i]->setTearingThreshold(max_tension);
    }
}

std::vector<std::pair<long int, long int>>
FabricTearer::recordGindexPairs()
{
    std::vector<std::pair<long int, long int>> gpairs;
    for (int i = 0; i < edges.size(); ++i)
    {
        long int gidx_beg = Gindex(edges[i]->beg);
        long int gidx_end = Gindex(edges[i]->end);
        gpairs.push_back(std::make_pair(gidx_beg,gidx_end));
    }
    return gpairs;
}

void FabricTearer::readGindexPairs(
        POINT** gpoints,
        const std::vector<std::pair<long int, long int>>& gpairs)
{
    printf("gpairs.size() = %d\n",gpairs.size());
    
    clearEdges();
    for (int i = 0; i < gpairs.size(); ++i)
    {
        long int gidx_beg = gpairs[i].first;
        long int gidx_end = gpairs[i].second;
        edges.push_back(
                new FabricEdge(gpoints[gidx_beg],gpoints[gidx_end]));
    }
}

std::vector<double> FabricTearer::recordRestingEdgeLengths()
{
    std::vector<double> restlengths;
    for (int i = 0; i < edges.size(); ++i)
    {
        double length = edges[i]->getLength();
        edges[i]->setRestLength(length);
        restlengths.push_back(length);
    }
    return restlengths;
}

void FabricTearer::readRestingEdgeLengths(
        const std::vector<double>& restlengths)
{
    printf("edges.size() = %d\n",edges.size());
    printf("restlengths.size() = %d\n",restlengths.size());
    
    for (int i = 0; i < restlengths.size(); ++i)
        edges[i]->setRestLength(restlengths[i]);
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
        long int itear = tear_idx[i];
        FabricEdge* e = edges[itear];

        //TODO: create/propagate the tear
        printf("tear_idx[%d] = %ld\n",i,itear);
    }

    printf("processTearingEvents() not implemented\n");
}

void FabricTearer::printEdges()
{
    for (int i = 0; i < edges.size(); ++i)
    {
        printf("i = %d   ",i);
        edges[i]->print();
    }
}

////////////////////////////////////
///////     FabricEdge       ///////
////////////////////////////////////

void FabricEdge::setRestLength(double l)
{
    length0 = l;
}

void FabricEdge::setSpringConstant(double k)
{
    ks = k;
}

void FabricEdge::setTearingThreshold(double T)
{
    tear_threshold = T;
}

double FabricEdge::getLength()
{
    computeLength();
    return length;
}

double FabricEdge::getTension()
{
    computeTension();
    return tension;
}

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
        double T = ks*(1.0 - length0/length)*disp[i];
        sqr_tension += T*T;
    }

    tension = sign*sqrt(sqr_tension);
}

bool FabricEdge::underTension()
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

void FabricEdge::print()
{
    printf("tension = %g   tear_threshold = %g   \
            length = %g   length0 = %g\n",
            tension,tear_threshold,length,length0);
}


