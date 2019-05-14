#include "BVH.h"


//TODO: I Don't like this, but it's too early
//      to tell how/where these constants should
//      be set, and need it for testing in the
//      short term
double BVH::proximityPad = 1.0e-03;


const bool BVH::isEmpty() const noexcept
{
    return (!this->root) ? true : false;
}

BVH_Node* const BVH::getRoot() const noexcept
{
    return root;
}

BVH_Node* BVH::createLeafNode(Hse* h)
{
    assert(h);
    BVH_Node* node = new LeafNode(h);
    node->expandBV(proximityPad);
    return node;
}

BVH_Node* BVH::createInternalNode(BVH_Node* lc, BVH_Node* rc)
{
    assert(lc && rc);
    return new InternalNode(lc,rc);
}

/*
BVH::BVH(Front* front, bool draw)
    : outdir{std::string(OutName(front))},
    drawbool{draw}
{
    constructLeafNodes(front->interf);
    buildHeirarchy();
}
*/

BVH::BVH(Front* front)
    : outdir{std::string(OutName(front))}
{
    constructLeafNodes(front->interf);
}

void BVH::constructLeafNodes(std::vector<Hse*> hseList)
{
    assert(leaves.empty());
    std::vector<Hse*>::iterator it = hseList.begin();
    for( it; it != hseList.end(); ++it )
    {
        BVH_Node* leaf = BVH::createLeafNode(*it);
        leaves.push_back(leaf);
    }
}

void BVH::buildHeirarchy()
{
    assert( !leaves.empty() );

    initChildren();
    while( children.size() != 1 )
    {
        //alternate sorting direction at each level
        if( sort_iter % 2 == 0 )
        {
            std::reverse(children.begin(),children.end());
        }
        drawHeirarchyLevel();
        constructParentNodes();
    }
   
    Point_Node_Vector().swap(children);
    assert( root != nullptr );
}

//TODO: Add debugging printouts to determine what
//      the wave/hsbdry types for each element.
//      When we figure out all the correct types,
//      should put all the type enums into scoped
//      enum classes. Otherwise it is impossible to
//      distinguish the two (is there a difference?).
void BVH::constructLeafNodes(INTERFACE* intfc)
{
	SURFACE** surfaces = intfc->surfaces;
    processSurfaces(surfaces);

    CURVE** curves = intfc->curves;
    processCurves(curves);
}

//TODO: Do we need a Hse factory?
void BVH::processSurfaces(SURFACE** surf)
{
    for( surf; surf && *surf; ++surf )
	{
	    if (is_bdry(*surf)) continue;
	    
        TRI* tri;
        surf_tri_loop(*surf,tri)
	    {
            if (wave_type(*surf) == MOVABLE_BODY_BOUNDARY || 
                    wave_type(*surf) == NEUMANN_BOUNDARY)
            {
                Hse* h = new HsTri(tri,HseTag::RIGIDBODY); 
                leaves.push_back(BVH::createLeafNode(h));
                num_tris++;
            }
            else if( wave_type(*surf) == ELASTIC_BOUNDARY )
            {
                Hse* h = new HsTri(tri,HseTag::FABRIC); 
                leaves.push_back(BVH::createLeafNode(h));
                num_tris++;
            }
	    }
	}
}

void BVH::processCurves(CURVE** curve)
{
	for( curve; curve && *curve; ++curve ) 
	{
        //TODO: What about ELASTIC_STRING?
	    if (hsbdry_type(*curve) != STRING_HSBDRY
                && hsbdry_type(*curve) != MONO_COMP_HSBDRY )
            continue; 

        BOND* bond;
	    curve_bond_loop(*curve,bond)
	    {                        
            if( hsbdry_type(*curve) == STRING_HSBDRY )
            {
                Hse* b = new HsBond(bond,HseTag::STRING);
                leaves.push_back(BVH::createLeafNode(b));
		        num_bonds++;
            }
            else if( hsbdry_type(*curve) == MONO_COMP_HSBDRY )
            {
                //TODO: Is this getting counted twice,
                //      first time being with the triangles?
                Hse* b = new HsBond(bond,HseTag::FABRIC);
                leaves.push_back(BVH::createLeafNode(b));
		        num_bonds++;
            }
	    }
	}
}

void BVH::initChildren()
{
    sort_iter = 0;
    children.reserve(leaves.size());
    children = getLeafSortingData();
    sortChildren();
}

const Point_Node_Vector BVH::getLeafSortingData() const
{
    Point_Node_Vector leafdata;

    std::vector<BVH_Node*>::const_iterator it;
    for( it = leaves.cbegin(); it != leaves.cend(); ++it )
    {
        auto node = *it;
        Point_with_Node bvctr_node_pair(node->getBV().Centroid(),node);
        leafdata.push_back(bvctr_node_pair);
    }
    return leafdata;
}

/*
const Point_Node_Vector BVH::getSortedLeafData() const
{
    Point_Node_Vector leafdata(getLeafSortingData());
    CGAL::hilbert_sort(leafdata.begin(),leafdata.end(),hst);
    return leafdata;
}
*/

void BVH::sortChildren()
{
    assert( !children.empty() );
        
    sort_iter++;
    if( children.size() == 1 )
    {
        root = children[0].second;
        drawHeirarchyLevel();
    }
    else
    {
        CGAL::hilbert_sort(children.begin(),children.end(),hst);
    }
}

void BVH::constructParentNodes()
{
    Point_Node_Vector parents;
    parents.reserve(children.size()/2 + 1);
    
    //greedily pair off sorted children
    for( int i = 0; i < children.size()-1; i += 2 )
    {
        auto lc = children[i].second;
        auto rc = children[i+1].second;
        auto p = BVH::createInternalNode(lc,rc);
        Point_with_Node bvctr_node_pair(p->getBV().Centroid(),p);
        parents.push_back(bvctr_node_pair);
    }

    //if an odd number of children, move the unpaired
    //node up to parent level unchanged
    if( children.size() % 2 != 0 )
    {
        auto oc = children[children.size()-1].second;
        Point_with_Node bvctr_node_pair(oc->getBV().Centroid(),oc);
        parents.push_back(bvctr_node_pair);
    }

    std::swap(parents,children);
    children.shrink_to_fit();
    sortChildren();
}

void BVH::drawHeirarchyLevel() const
{
    if( drawbool == true )
        writeHilbertCurveFiles(sort_iter);
}

void BVH::setDrawBool(bool draw)
{
    drawbool = draw;
}

void BVH::setDrawDirectory(std::string dir)
{
    outdir = dir;
}

