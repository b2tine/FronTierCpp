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

const std::weak_ptr<BVH_Node> BVH::getRoot() const noexcept
{
    return std::weak_ptr<BVH_Node>(root);
}


//TODO: Is std::move() Hindering RVO in these factory functions?
//      Or are the RVO conditions not met?
//      i.e. does it work for polymorphic types?

std::shared_ptr<BVH_Node>
BVH::createLeafNode(Hse* h)
{
    auto node = std::make_shared<LeafNode>(h);
    node->expandBV(proximityPad);
    return std::move(node);
}


std::shared_ptr<BVH_Node>
BVH::createInternalNode(std::shared_ptr<BVH_Node> lc,
        std::shared_ptr<BVH_Node> rc)
{
    auto node = std::make_shared<InternalNode>(lc,rc);
    node->setChildren(lc,rc);
    return std::move(node);
}


BVH::BVH(const Front* const front)
{
    constructLeafNodes(front->interf);
    buildHeirarchy();
}


void BVH::constructLeafNodes(const INTERFACE* const intfc)
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

                leaves.push_back(BVH::createLeafNode(
                            new HsTri(tri,HseTag::RIGIDBODY)));
            }
            else
            {
                //TODO: what is the wave_type for fabric?
                //      MONO_COMP_HSBDRY? ELASTIC?
                leaves.push_back(BVH::createLeafNode(
                            new HsTri(tri,HseTag::FABRIC)));
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
    assert(!leaves.empty());

    initChildren();
    while( children.size() > 2 )
    {
        constructParentNodes();
    }
    constructRootNode();    
}

void BVH::initChildren()
{
    sort_iter = 0;
    Point_Node_Vector().swap(children);
    children.reserve(leaves.size());
    children = getLeafSortingData();
    sortChildren();
}

Point_Node_Vector BVH::getLeafSortingData() const
{
    Point_Node_Vector leafdata;

    std::vector<std::shared_ptr<BVH_Node>>::const_iterator it;
    for( it = leaves.cbegin(); it != leaves.cend(); ++it )
    {
        auto node = *it;
        Point_with_Node bvctr_node_pair(node->getBV().Centroid(),node);
        leafdata.push_back(bvctr_node_pair);
    }
    return leafdata;
}

void BVH::sortChildren()
{
    assert(!children.empty());
    CGAL::hilbert_sort(children.begin(),children.end(),hst);
    sort_iter++;
}



/*
void BVH::clearVectors()
{
    std::vector<std::shared_ptr<BVH_Node>>().swap(leaves);
    Point_Node_Vector().swap(children);

    if(num_leaves != 0) 
    {
        leaves.reserve(num_leaves);
        children.reserve(num_leaves);
    }
}
*/


void BVH::constructParentNodes()
{
    assert(!children.empty());

    //alternate sweep direction at each level
    if( sort_iter % 2 == 0 )
    {
        std::reverse(children.begin(),children.end());
    }

    Point_Node_Vector parents;
    
    //greedily pair off sorted children
    for( int i = 0; i < children.size()-1; i += 2 )
    {
        auto lc = children[i].second;
        auto rc = children[i+1].second;
        auto p = BVH::createInternalNode(lc,rc);
        Point_with_Node bvctr_node_pair(p->getBV().Centroid(),p);
        parents.push_back(bvctr_node_pair);
    }

    //if odd number of leafnodes on the first pass
    if( children.size() % 2 != 0 )
    {
        auto oc = children[children.size()-1].second;
        auto p = BVH::createInternalNode(oc,oc);
        Point_with_Node bvctr_node_pair(p->getBV().Centroid(),p);
        parents.push_back(bvctr_node_pair);
    }

    std::swap(parents,children);
    sortChildren();

    Point_Node_Vector().swap(parents);
}

void BVH::constructRootNode()
{
    assert(children.size() == 2);
    auto lc = children[0].second;
    auto rc = children[1].second;
    root = BVH::createInternalNode(lc,rc);
    assert(root);
    sort_iter = 0;
}

const Point_Node_Vector BVH::getSortedLeafData() const
{
    Point_Node_Vector leafdata(getLeafSortingData());
    CGAL::hilbert_sort(leafdata.begin(),leafdata.end(),hst);
    return leafdata;
}

//Testing function
void BVH::buildTester(std::vector<Hse*> hseList)
{
    assert(children.empty());
    std::vector<Hse*>::iterator it = hseList.begin();
    for( it; it != hseList.end(); ++it )
    {
        auto leaf = BVH::createLeafNode(*it);
        Point_with_Node ctr_bv_pair(leaf->getBV().Centroid(),leaf);
        children.push_back(ctr_bv_pair);
    }

    while( children.size() > 2 )
    {
        constructParentNodes();
    }
    constructRootNode();    
}



