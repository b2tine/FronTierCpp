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
}


void BVH::constructBVH(const Front* const front)
{
    constructLeafNodes(front->interf);
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
                //hseList.push_back(new HsTri(tri,HseTag::RIGIDBODY));
                //hseList.push_back(new CD_TRI(tri, "tris_rigid"));
            }
            else
            {
                leaves.push_back(BVH::createLeafNode(
                            new HsTri(tri,HseTag::FABRIC)));
                //hseList.push_back(new HsTri(tri,HseTag::FABRIC));
                //hseList.push_back(new CD_TRI(tri, "tris"));
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
            //hseList.push_back(new HsBond(b,HseTag::STRING));
            //hseList.push_back(new CD_BOND(b,m_dim, "string_bond"));
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

    //sortNodes(children,hst);
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

    //if odd number of leafnodes
    if( children.size() % 2 != 0 )
    {
        auto oc = children[children.size()-1].second;
        auto p = BVH::createInternalNode(oc,oc);
        Point_with_Node ctr_bv_pair(p->getBV().Centroid(),p);
        parents.push_back(ctr_bv_pair);
    }

    //sortNodes(parents,hst);
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
}

/*
void sortNodes(Point_Node_Vector& pn_vec,
        const BV_HilbertSortingTraits& traits)
{
    CGAL::hilbert_sort(pn_vec.begin(),pn_vec.end(),traits);
}
*/




