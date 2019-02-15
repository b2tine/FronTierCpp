#include "BVH.h"

std::shared_ptr<LeafNode>
BVH::createLeafNode(Hse* h)
{
    return std::make_shared<LeafNode>(h);
}


std::shared_ptr<InternalNode>
BVH::createInternalNode(std::shared_ptr<BVH_Node> lc,
        std::shared_ptr<BVH_Node> rc)
{
    auto node = std::make_shared<InternalNode>(lc,rc);
    node->setChildren(lc,rc);
    return std::move(node);
}


void BVH::sortChildNodes()
{
    assert(!children.empty());
    CGAL::hilbert_sort(children.begin(),children.end(),hst);
    sort_iter++;
}


BVH::BVH(const Front* const front)
{
    constructLeafNodes(front->interf);
    while( children.size() > 2 )
    {
        constructParentNodes();
    }

    constructRootNode();    
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


//Does this need to be called more than once if
//all hypersurf elements remain in the leaves vector?
void BVH::constructLeafNodes(const INTERFACE* const intfc)
{
	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;
	
	//clearVectors();

    int n_tri = 0;
    int n_bond = 0;

    intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    //unsort_surface_point(*s);
	    surf_tri_loop(*s,tri)
	    {
            if (wave_type(*s) == MOVABLE_BODY_BOUNDARY || 
                    wave_type(*s) == NEUMANN_BOUNDARY)
            {

                leaves.push_back(BVH::createLeafNode(
                            new HsTri(tri,HseTag::RIGIDBODY)));
            }
            else
            {
                leaves.push_back(BVH::createLeafNode(
                            new HsTri(tri,HseTag::FABRIC)));
            }
            
            auto node = leaves.back();
            Point_with_Node ctr_bv_pair(node->getBV().Centroid(),node);
            children.push_back(ctr_bv_pair);
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

            auto node = leaves.back();
            Point_with_Node ctr_bv_pair(node->getBV().Centroid(),node);
            children.push_back(ctr_bv_pair);
		    n_bond++;
	    }
	}

	//makeSet(hseList);
	//createImpZoneForRG(intfc);
	//setDomainBoundary(intfc->table->rect_grid.L,intfc->table->rect_grid.U);
	if( debugging("BVH") )
    {
	    printf("%d num of tris, %d num of bonds\n",n_tri,n_bond);
	    printf("%lu number of elements is assembled\n",leaves.size());
	    //printf("%lu number of elements is assembled\n",hseList.size());
	}

    sort_iter = 0;
    sortChildNodes();
}

/*
void BVH::clearVectors()
{
    std::vector<std::shared_ptr<BVH_Node>>().swap(leaves);
    Point_Node_Vector().swap(children);

    if(numLeaves != 0) 
    {
        leaves.reserve(numLeaves);
        children.reserve(numLeaves);
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
        Point_with_Node ctr_bv_pair(p->getBV().Centroid(),p);
        parents.push_back(ctr_bv_pair);
    }

    //if odd number of leafnodes on the first pass
    if( children.size() % 2 != 0 )
    {
        auto oc = children[children.size()-1].second;
        auto p = BVH::createInternalNode(oc,oc);
        Point_with_Node ctr_bv_pair(p->getBV().Centroid(),p);
        parents.push_back(ctr_bv_pair);
    }

    std::swap(parents,children);
    sortChildNodes();
    //Point_Node_Vector().swap(parents);
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


const std::weak_ptr<InternalNode> BVH::getRoot() const
{
    return std::weak_ptr<InternalNode>(root);
}




