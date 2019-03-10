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

BVH_Node* BVH::getRoot() const noexcept
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


BVH::BVH(Front* front)
{
    constructLeafNodes(front->interf);
    buildHeirarchy();
}


void BVH::constructLeafNodes(INTERFACE* intfc)
{
	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;
	
    int n_tri = 0;
    int n_bond = 0;

    intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    
        surf_tri_loop(*s,tri)
	    {
            //TODO: what other wave_types/boundaries need
            //      to be considered?
            if (wave_type(*s) == MOVABLE_BODY_BOUNDARY || 
                    wave_type(*s) == NEUMANN_BOUNDARY)
            {
                Hse* h = new HsTri(tri,HseTag::RIGIDBODY); 
                leaves.push_back(BVH::createLeafNode(h));
            }
            else
            {
                //TODO: what is the wave_type for fabric?
                //      MONO_COMP_HSBDRY? ELASTIC?
                Hse* h = new HsTri(tri,HseTag::FABRIC); 
                leaves.push_back(BVH::createLeafNode(h));
            }
            n_tri++;
	    }
	}

	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 
	    curve_bond_loop(*c,b)
	    {
            leaves.push_back(BVH::createLeafNode(
                        new HsBond(b,HseTag::STRING)));
		    n_bond++;
	    }
	}

	if( debugging("BVH") )
    {
	    printf("%d num tris, %d num bonds\n",n_tri,n_bond);
	    printf("%lu total elements\n",leaves.size());
	}

}

void BVH::buildHeirarchy()
{
    assert( !leaves.empty() );

    initChildren();
    while( children.size() > 2 )
    {
        constructParentNodes();
    }
    
    constructRootNode();    
    Point_Node_Vector().swap(children);
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

const Point_Node_Vector BVH::getSortedLeafData() const
{
    Point_Node_Vector leafdata(getLeafSortingData());
    CGAL::hilbert_sort(leafdata.begin(),leafdata.end(),hst);
    return leafdata;
}

void BVH::sortChildren()
{
    assert( !children.empty() );
    CGAL::hilbert_sort(children.begin(),children.end(),hst);
    sort_iter++;
    //TODO: add option and function to write the children to
    //      output file at each level of the heirarchy.
}

void BVH::constructParentNodes()
{
    //alternate sorting direction at each level
    if( sort_iter % 2 == 0 )
    {
        std::reverse(children.begin(),children.end());
    }

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

    //if odd number of children, the unpaired one gets bumped up to parents 
    if( children.size() % 2 != 0 )
    {
        auto oc = children[children.size()-1].second;
        Point_with_Node bvctr_node_pair(oc->getBV().Centroid(),oc);
        parents.push_back(bvctr_node_pair);
    }

    std::swap(parents,children);
    Point_Node_Vector().swap(parents);
    children.shrink_to_fit();
    sortChildren();
}

void BVH::constructRootNode()
{
    assert( children.size() == 2 );
    auto lc = children[0].second;
    auto rc = children[1].second;
    root = BVH::createInternalNode(lc,rc);
}


//Testing function
void BVH::buildTester(std::vector<Hse*> hseList)
{
    assert(children.empty());
    std::vector<Hse*>::iterator it = hseList.begin();
    for( it; it != hseList.end(); ++it )
    {
        BVH_Node* leaf = BVH::createLeafNode(*it);
        Point_with_Node ctr_bv_pair(leaf->getBV().Centroid(),leaf);
        children.push_back(ctr_bv_pair);
    }

    while( children.size() > 2 )
    {
        constructParentNodes();
    }
    constructRootNode();    
}


