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
    : outdir{std::string(OutName(front))},
    timestep{&front->step}
{
    drawdir = outdir + "/BVH/";
    constructLeafNodes(front->interf);
    DrawUnlock();
    buildHeirarchy();
}

//TODO: For testing only right now, but may be a
//      good starting point for rebuilding/refitting.
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

//TODO: Add debugging printouts to determine what
//      the wave/hsbdry types for each element.
//      When we figure out all the correct types,
//      should put all the type enums into scoped
//      enum classes. Otherwise it is impossible to
//      distinguish the two (is there a difference?).
void BVH::constructLeafNodes(INTERFACE* intfc)
{
    assert(leaves.empty());
	SURFACE** surfaces = intfc->surfaces;
    processSurfaces(surfaces);

    CURVE** curves = intfc->curves;
    processCurves(curves);
    assert(!leaves.empty());
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

const Point_Node_Vector BVH::getLeafSortingData() const
{
    Point_Node_Vector leafdata;
    leafdata.reserve(leaves.size());

    std::vector<BVH_Node*>::const_iterator it;
    for( it = leaves.cbegin(); it != leaves.cend(); ++it )
    {
        auto node = *it;
        Point_with_Node bvctr_node_pair(node->getBV().Centroid(),node);
        leafdata.push_back(bvctr_node_pair);
    }
    return leafdata;
}

void BVH::updateHeirarchy()
{
    //TODO: Refit bounding volumes etc.
    std::vector<BVH_Node*>::iterator it;
    for( it = leaves.begin(); it != leaves.end(); ++it )
    {
        auto node = *it;
        node->refitBV();
    }    
    initChildren();
    sortChildren();
}

void BVH::buildHeirarchy()
{
    initChildren();
    while( children.size() != 1 )
    {
        //alternate sorting direction at each level
        if( sort_iter % 2 == 1 )
        {
            std::reverse(children.begin(),children.end());
        }
        drawHeirarchyLevel();
        constructParentNodes();
        sortChildren();
        sort_iter++;
    }
    constructRootNode();
    Point_Node_Vector().swap(children);
}

void BVH::initChildren()
{
    sort_iter = 0;
    children = getLeafSortingData();
    sortChildren();
}

void BVH::sortChildren()
{
    //sort_iter++;
    if( children.size() == 1 )
    {
        return;
        //root = children[0].second;
        //drawHeirarchyLevel();
    }
    //CGAL::hilbert_sort(children.begin(),children.end(),hst);
    CGAL::spatial_sort(children.begin(),children.end(),
            hst,CGAL::Hilbert_sort_median_policy());
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
        auto orphan = children[children.size()-1].second;
        Point_with_Node bvctr_node_pair(orphan->getBV().Centroid(),orphan);
        parents.push_back(bvctr_node_pair);
    }

    std::swap(parents,children);
    children.shrink_to_fit();
}


void BVH::constructRootNode()
{
    root = children[0].second;
    assert( root != nullptr );
    drawHeirarchyLevel();
    isDrawTime = false;
}
