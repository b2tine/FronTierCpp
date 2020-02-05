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
        //std::unordered_set<TRI*> tri_set;

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

                //TODO: not sure if this is sufficient to ensure
                //      edge has been visited
                //
                //check if edge has already been processed
                if (visited1 && visited2)
                    continue;

                long int left_tri_gindex = Gindex(tri);
                long int right_tri_gindex = -1;

                TRI* ntri = nullptr;
                for (int i = 0 ; i < 3; ++i)
                {
                    ntri = Tri_neighbor(tri)[i].tri;

                    if (ntri != nullptr)
                    {
                        if (Point_on_tri(ntri,p1) && Point_on_tri(ntri,p2))
                        {
                            right_tri_gindex = Gindex(ntri);
                            break;
                        }
                    }

                    if (i == 2 && right_tri_gindex == -1)
                    {
                        ntri = nullptr;
                        //printf("can't find neighboring tri to %ld\n",
                          //      left_tri_gindex);
                        //clean_up(ERROR);
                    }
                }

                FabricEdge* fedge = new FabricEdge(p1,p2,tri,ntri); 
                edges.push_back(fedge);
            }
        }

        /*
        BOND* b;
        CURVE** c;
        BOND_TRI** btris;

        surf_pos_curve_loop(*s,c)
        {
            curve_bond_loop(*c,b)
            {
                bond_btri_loop(b,btris)
                {

                }
            }
        }
        */
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

//TODO: join with other record functions and return
//      an object containing the necessary vectors
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

std::vector<std::pair<long int, long int>>
FabricTearer::recordGindexTriPairs() const
{
    std::vector<std::pair<long int, long int>> gpairs;
    for (int i = 0; i < edges.size(); ++i)
    {
        FabricEdge* e  = edges[i];
       
        long int gidx_left = -1;
        long int gidx_right = -1;
        
        if (e->left_tri != nullptr)
            gidx_left = Gindex(e->left_tri);
        if (e->right_tri != nullptr)
            gidx_right = Gindex(e->right_tri);
        
        if (gidx_left == -1 && gidx_right == -1)
        {
            printf("edges[%d] has no incident tris\n",i);
            clean_up(ERROR);
        }

        gpairs.push_back(std::make_pair(gidx_left,gidx_right));
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

std::vector<long int> FabricTearer::recordGindexWeakPoints() const
{
    return std::vector<long int>(weakpt_idx.cbegin(),weakpt_idx.cend());
}

void FabricTearer::readEdgeData(
        POINT** gpoints,
        const std::vector<std::pair<long int, long int>>& gindex_ppairs,
        TRI** gtris,
        const std::vector<std::pair<long int, long int>>& gindex_tpairs)
{
        //printf("edges.size() = %d\n",edges.size());
        //printf("POINT gindex_pairs.size() = %d\n",gindex_pairs.size());
    
    clearEdges();
    for (int i = 0; i < gindex_ppairs.size(); ++i)
    {
        long int gidx_beg = gindex_ppairs[i].first;
        long int gidx_end = gindex_ppairs[i].second;
        
        long int gidx_left = gindex_tpairs[i].first;
        if (gidx_left != -1)
            edges[i]->left_tri = gtris[gidx_left];

        long int gidx_right = gindex_tpairs[i].second;
        if (gidx_right != -1)
            edges[i]->right_tri = gtris[gidx_right];
        
        edges.push_back(
                new FabricEdge(gpoints[gidx_beg],gpoints[gidx_end],
                                gtris[gidx_left],gtris[gidx_right]));
    }
}

/*
void FabricTearer::readGindexTriPairs(
        TRI** gtris,
        const std::vector<std::pair<long int, long int>>& gindex_pairs)
{
        //printf("edges.size() = %d\n",edges.size());
        //printf("TRI gindex_pairs.size() = %d\n",gindex_pairs.size());

    for (int i = 0; i < gindex_pairs.size(); ++i)
    {
        long int gidx_left = gindex_pairs[i].first;
        if (gidx_left != -1)
            edges[i]->left_tri = gtris[gidx_left];

        long int gidx_right = gindex_pairs[i].second;
        if (gidx_right != -1)
            edges[i]->right_tri = gtris[gidx_right];
    }
}

void FabricTearer::readRestingEdgeLengths(
        const std::vector<double>& restlengths)
{
        //printf("edges.size() = %d\n",edges.size());
        //printf("restlengths.size() = %d\n",restlengths.size());
    
    for (int i = 0; i < restlengths.size(); ++i)
        edges[i]->setRestLength(restlengths[i]);
}
*/

//TODO: join with readEdges() into single function
void FabricTearer::readGindexWeakPoints(
        const std::vector<long int>& gindex_weakpts)
{
    weakpt_idx.clear();
    weakpt_idx.insert(gindex_weakpts.begin(),gindex_weakpts.end());
    for (int i = 0; i < edges.size(); ++i)
    {
        FabricEdge* e = edges[i];
        if (isWeakPoint(e->beg) || isWeakPoint(e->end))
            e->setWeakPointFlag(true);
    }
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
        if (e->checkForTear())
            tear_idx.push_back(i);
    }
}

void FabricTearer::processTearingEvents()
{
    printf("Entering processTearingEvents()\n");
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
    
    TRI* left_tri = e->left_tri;
    TRI* right_tri = e->right_tri;

    TRI** tri_list;
    POINT *p1,*p2,*p3,*p4,*newp;
    SURFACE *surf;
    double coords[MAXD];
    int i,j,side_left,side_right;
    int num_tris;
    INTERFACE *save_intfc = current_interface();

    //IMPLEMENTATION HERE

    //These are needed for functions to be called
    surf = left_tri->surf;
    INTERFACE* intfc = surf->interface;
    set_current_interface(intfc);
    printf("intfc = %p  surf = %p\n",intfc,surf);

    // Find the side of the left_tri neighboring the right_tri
    for (i = 0; i < 3; ++i)
    {
        if (Tri_on_side(left_tri,i) == right_tri)
        {
            side_left = i;
            break;
        }
    }
    // Find the side of the right_tri neighboring the left_tri
    for (i = 0; i < 3; ++i)
    {
        if (Tri_on_side(right_tri,i) == left_tri)
        {
            side_right = i;
            break;
        }
    }
    // Check consistency
    printf("side_left  = %d\n",side_left);
    printf("side_right = %d\n",side_right);

    // Get the two points of the side two tris neighboring each other
    p1 = Point_of_tri(left_tri)[side_left];
    p2 = Point_of_tri(right_tri)[side_right];
    // Should be consistent
    printf("p1 = %p  rp1 = %p\n",p1,Point_of_tri(right_tri)[(side_right+1)%3]);
    printf("p2 = %p  rp2 = %p\n",p2,Point_of_tri(left_tri)[(side_left+1)%3]);

    // Two opposite points
    p3 = Point_of_tri(left_tri)[(side_left+2)%3];
    p4 = Point_of_tri(right_tri)[(side_right+2)%3];
    // The two points the two tris share
    printf("p1 = %p\n",p1);
    printf("p2 = %p\n",p2);
    // The two points the two tris not share
    printf("p3 = %p\n",p3);
    printf("p4 = %p\n",p4);
    printf("\n");

    // Compute coordinates of mid-point
    for (i = 0; i < 3; ++i)
    {
        coords[i] = 0.5*(Coords(p1)[i] + Coords(p2)[i]);
    }

    printf("midpoint coords = %f %f %f\n",coords[0],coords[1],coords[2]);
    printf("\n");

    // Compute direction vectors to p3 and p4 from midpoint for
    // displacing the new points
    double vec_mp2p3[3];
    double vec_mp2p4[3];
    for (int i = 0; i < 3; ++i)
    {
        vec_mp2p3[i] = Coords(p3)[i] - coords[i];
        vec_mp2p4[i] = Coords(p4)[i] - coords[i];
    }
    double mag_mp2p3 = Mag3d(vec_mp2p3);
    double mag_mp2p4 = Mag3d(vec_mp2p4);
    for (int i = 0; i < 3; ++i)
    {
        vec_mp2p3[i] /= mag_mp2p3;
        vec_mp2p4[i] /= mag_mp2p4;
    }


    // Create a point at the middle coordinate and insert in side of left_tri
    newp = Point(coords);
    printf("newp = %p\n",newp);
    insert_point_in_tri_side(newp,side_left,left_tri,surf);

    // Get the four tris after splitting
    num_tris = set_tri_list_around_point(newp,left_tri,&tri_list,intfc);
    printf("Around newp: num_tris = %d\n",num_tris);

    // Get the side of each tri for tearing
    int n,side[4];
    for (i = 0; i < num_tris; ++i)
    {
        side[i] = -1;
        for (j = 0; j < 3; ++j)
        {
            if ((Point_of_tri(tri_list[i])[j] == newp &&
                 (Point_of_tri(tri_list[i])[(j+1)%3] == p3 ||
                  Point_of_tri(tri_list[i])[(j+1)%3] == p4))
                ||
                ((Point_of_tri(tri_list[i])[j] == p3 ||
                  Point_of_tri(tri_list[i])[j] == p4) &&
                 Point_of_tri(tri_list[i])[(j+1)%3] == newp))
                side[i] = j;

        }
        printf("Tearing side of tri_list[%d]: %d\n",i,side[i]);
    }
    for (i = 0; i < num_tris; ++i)
        printf("Gindex %d: %d\n",i,Gindex(tri_list[i]));

    // Divide the four tris into two sets, each on one tearing side
    TRI *tris1[2],*tris2[2];
    n = 0;
    for (i = 0; i < num_tris; ++i)
        for (j = 0; j < 3; ++j)
            if (Point_of_tri(tri_list[i])[j] == p1)
                tris1[n++] = tri_list[i];
    printf("n1 = %d\n",n);
    n = 0;
    for (i = 0; i < num_tris; ++i)
        for (j = 0; j < 3; ++j)
            if (Point_of_tri(tri_list[i])[j] == p2)
                tris2[n++] = tri_list[i];
    printf("n2 = %d\n",n);
    printf("tris1 = %p %p\n",tris1[0],tris1[1]);
    printf("tris2 = %p %p\n",tris2[0],tris2[1]);

    // Duplicate the inserted point newp to split
    POINT *newp1 = Point(Coords(newp));
    // Reassign the newp on side 2 to newp1
    for (i = 0; i < 2; ++i)
        for (j = 0; j < 3; ++j)
            if (Point_of_tri(tris2[i])[j] == newp)
                    Point_of_tri(tris2[i])[j] = newp1;

    Gindex(newp) = intfc->max_point_gindex;
    Gindex(newp1) = intfc->max_point_gindex + 1;
    intfc->max_point_gindex += 2;
    printf("Gindex(newp) = %d\n",Gindex(newp));
    printf("Gindex(newp1) = %d\n",Gindex(newp1));

    // Open the tearing side of each tri
    printf("p3 = %p\n",p3);
    printf("p4 = %p\n",p4);
    printf("newp = %p\n",newp);
    printf("newp1 = %p\n",newp1);
    for (i = 0; i < num_tris; ++i)
    {
        Tri_on_side(tri_list[i],side[i]) = NULL;
        // Check: ps and pe must one as newp/newp1 and the other as p3/p4
        printf("ps = %p  pe = %p\n",Point_of_tri(tri_list[i])[side[i]],
                Point_of_tri(tri_list[i])[(side[i]+1)%3]);
    }

    // I already checked
    printf("\nnewp = %p\n",newp);
    printf("tris1[0] pts: %p %p %p\n",Point_of_tri(tris1[0])[0],
                Point_of_tri(tris1[0])[1],Point_of_tri(tris1[0])[2]);
    printf("tris1[0] sides: %p %p %p\n",Tri_on_side(tris1[0],0),
                Tri_on_side(tris1[0],1),Tri_on_side(tris1[0],2));
    printf("tris1[1] pts: %p %p %p\n",Point_of_tri(tris1[1])[0],
                Point_of_tri(tris1[1])[1],Point_of_tri(tris1[1])[2]);
    printf("tris1[1] sides: %p %p %p\n",Tri_on_side(tris1[1],0),
                Tri_on_side(tris1[1],1),Tri_on_side(tris1[1],2));

    printf("\nnewp1 = %p\n",newp1);
    printf("tris2[0] pts: %p %p %p\n",Point_of_tri(tris2[0])[0],
                Point_of_tri(tris2[0])[1],Point_of_tri(tris2[0])[2]);
    printf("tris2[0] sides: %p %p %p\n",Tri_on_side(tris2[0],0),
                Tri_on_side(tris2[0],1),Tri_on_side(tris2[0],2));
    printf("tris2[1] pts: %p %p %p\n",Point_of_tri(tris2[1])[0],
                Point_of_tri(tris2[1])[1],Point_of_tri(tris2[1])[2]);
    printf("tris2[1] sides: %p %p %p\n",Tri_on_side(tris2[1],0),
                Tri_on_side(tris2[1],1),Tri_on_side(tris2[1],2));


        //printf("\nInterface before install tearing curve:\n");
    // The original interface
    //print_interface(intfc);
    FT_InstallSurfEdge(surf,MONO_COMP_HSBDRY);
    // The interface after inserting tearing curve
        //printf("Interface after install tearing curve:\n");
    //print_interface(intfc);

    // Separate new points
    for (int i = 0; i < 3; ++i)
    {
        //TODO: wrong direction vector
        //Coords(newp)[i] += 0.001*vec_mp2p3[i];
        //Coords(newp1)[i] += 0.001*vec_mp2p4[i];
    }


    // Make sure you restore the global interface
    set_current_interface(save_intfc);

    //clean_up(0);
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


FabricEdge::FabricEdge(POINT* p1, POINT* p2, TRI* tl, TRI* tr)
    : beg{p1}, end{p2}, left_tri{tl}, right_tri{tr}
{
    for (int i = 0; i < 3; ++i)
    {
        if (beg == Point_of_tri(left_tri)[i])
            length0 = left_tri->side_length0[i];
    } 
    assert(length0 > 0.0);
    
    computeLength();
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
    double coeff = 1.0;//temp value
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


