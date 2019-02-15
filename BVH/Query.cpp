#include "BVH.h"

using NodePair = std::pair<std::shared_ptr<BVH_Node>,
                            std::shared_ptr<BVH_Node>>;


const bool checkProximity(const BVH* A, const BVH* B)
{
    assert(A && B);
    auto nodeA = A->getRoot().lock();
    auto nodeB = B->getRoot().lock();
    assert(nodeA && nodeB);

    queryProximity(nodeA,nodeB);

    //force failure of checkProximityTest for now
    return false;
    
}


void queryProximity(std::shared_ptr<BVH_Node> nodeA,
        std::shared_ptr<BVH_Node> nodeB)
{
    //NodePair pair = std::make_pair(nodeA,nodeB);
    //std::stack<NodePair> stack(pair);
    std::stack<NodePair> stack(std::make_pair(nodeA,nodeB));

    while( !stack.empty() )
    {
        auto A = stack.top().first;
        auto B = stack.top().second;
        stack.pop();

        if( A->overlaps(B) )
        {
            if( A->isLeaf() && B->isLeaf() )
            {
                //TODO: implement this calculation
                //compute distance / check intersection
            }
            else
            {
                if( A->volume() < B->volume() )
                {
                    auto rc = B->getRightChild().lock();
                    stack.push(std::make_pair(A,rc));
                    auto lc = B->getLeftChild().lock();
                    stack.push(std::make_pair(A,lc));
                }
                else
                {
                    auto rc = A->getRightChild().lock();
                    stack.push(std::make_pair(rc,B));
                    auto lc = A->getLeftChild().lock();
                    stack.push(std::make_pair(lc,B));
                }
            }
        }
    }

}



