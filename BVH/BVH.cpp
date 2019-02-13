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


//void BVH::assembleHseListFromInterface(const INTERFACE* const intfc)
void BVH::constructLeafNodes(const INTERFACE* const intfc)
{
	SURFACE** s;
	CURVE** c;
	TRI *tri;
	BOND *b;
	
    int n_tri = 0;
    int n_bond = 0;

	clearLeafNodes();
	//clearHseList();
    //if(lastHseCount != 0) 
    if(lastCountLeaves != 0) 
    {
        leaves.reserve(lastCountLeaves);
        //hseList.reserve(lastHseCount);
    }
	
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
            ctrVec.push_back(node->getBV().Centroid());
            bvMap[ctrVec[n_tri++]] = node;
            //n_tri++;
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
            ctrVec.push_back(node->getBV().Centroid());
            bvMap[ctrVec[n_tri + n_bond++]] = node;
		    //n_bond++;
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

}


//void BVH::clearHseList()
void BVH::clearLeafNodes()
{
    int lastCountLeaves = leaves.size();
    ctrVec.clear();
    bvMap.clear();
    leaves.clear();
    /*
    int lastHseCount = hseList.size();
    for( int i = 0; i < hseList.size(); ++i)
    {
        delete hseList[i];
    }
    hseList.clear();
    */
}

void BVH::sortNodes()
{
    BV_HilbertSortingTraits hst;
    CGAL::hilbert_sort(ctrVec.begin(),ctrVec.end(),hst);
}



