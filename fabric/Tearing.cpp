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
FabricTearer::recordGindexPointPairs() const
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

//std::unordered_map<long int> FabricTearer::recordGindexWeakPoints() const
std::vector<long int> FabricTearer::recordGindexWeakPoints() const
{
    return std::vector<long int>(weakpt_idx.cbegin(),weakpt_idx.cend());
}

void FabricTearer::readGindexPointPairs(
        POINT** gpoints,
        const std::vector<std::pair<long int, long int>>& gindex_pairs)
{
    printf("gindex_pairs.size() = %d\n",gindex_pairs.size());
    
    clearEdges();
    for (int i = 0; i < gindex_pairs.size(); ++i)
    {
        long int gidx_beg = gindex_pairs[i].first;
        long int gidx_end = gindex_pairs[i].second;
        edges.push_back(
                new FabricEdge(gpoints[gidx_beg],gpoints[gidx_end]));
    }
}

void FabricTearer::readRestingEdgeLengths(
        const std::vector<double>& restlengths)
{
    printf("edges.size() = %d\n",edges.size());
    printf("restlengths.size() = %d\n",restlengths.size());
    
    for (int i = 0; i < restlengths.size(); ++i)
        edges[i]->setRestLength(restlengths[i]);
}

void FabricTearer::readGindexWeakPoints(
        const std::vector<long int>& gindex_weakpts)
{
    weakpt_idx.clear();
    weakpt_idx.insert(gindex_weakpts.begin(),gindex_weakpts.end());
}

//void FabricTearer::tearFabric(Front* front)
void FabricTearer::tearFabric()
{
    checkForTearingEvents();
    processTearingEvents();
        //FT_SetGlobalIndex(front);
}

void FabricTearer::checkForTearingEvents()
{
    tear_idx.clear();
    for (int i = 0; i < edges.size(); ++i)
    {
        FabricEdge* e = edges[i];

        if (isWeakPoint(e->beg) || isWeakPoint(e->end))
            e->setWeakPointFlag(true);

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
        
            //printf("tear_idx[%d] = %ld\n",i,itear);

        //TODO: need to be careful when/how updates are done
        //      to global indices and the vectors of objects
        //      that use them.
        if (!isWeakPoint(e->beg) && !isWeakPoint(e->end))
            createNewTear(e);
        else
            propagateTear(e);
    }
}

void FabricTearer::createNewTear(FabricEdge* e)
{
    //TODO: 
    //      1. split the 2 incident triangles of the rupturing edge
    //         into 2 triangles each using midpoint of the rupturing
    //         edge and the nonadjacent triangle vertices
    //
    //      2. move the 2 newly created points off the line of the
    //         ruptured edge, shortening the 4 new edges.
    //         Or else these will rupture immediately on next check.
    //
    //      3. update the global point and tri indices, the restlength
    //         vector, and the weakpt index to account for the new
    //         points/edges/tris
    
        //insert_point_in_tri_side();

        //weakpt_idx.insert(Gindex(e->beg));
        //weakpt_idx.insert(Gindex(e->end));

    printf("FabricTearer::createNewTear() not implemented yet\n");
    clean_up(0);
}

void FabricTearer::propagateTear(FabricEdge* e)
{
    //TODO: implementation
    //
    printf("FabricTearer::propagateTear() not implemented yet\n");
    clean_up(0);
};
bool FabricTearer::isWeakPoint(POINT* p)
{
    long int gindex = Gindex(p);
    return weakpt_idx.find(gindex) != weakpt_idx.end();
}

void FabricTearer::printEdges() const
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

FabricEdge::FabricEdge(POINT* p1, POINT* p2)
    : beg{p1}, end{p2}
{
    computeTension();
}

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

double FabricEdge::getLength() const
{
    return length;
}

double FabricEdge::getTension() const
{
    return tension;
}

bool FabricEdge::checkForTear()
{
    double coeff = 1.0;
    if (has_weakpt)
        coeff = weakpt_factor;
    return tension > coeff*tear_threshold;
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

bool FabricEdge::underTension() const
{
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

void FabricEdge::setWeakPointFlag(bool flag)
{
    has_weakpt = flag;
}

bool FabricEdge::hasWeakPoint() const
{
    return has_weakpt;
}

void FabricEdge::print() const
{
    printf("tension = %g   tear_threshold = %g   \
            length = %g   length0 = %g   has_weakpt = %s\n",
            tension,tear_threshold,length,length0,
            has_weakpt ? "YES" : "NO");
}


