#include "BVH.h"


using NodePair = std::pair<std::shared_ptr<BVH_Node>,
                            std::shared_ptr<BVH_Node>>;


const bool checkProximity(BVH* A, BVH* B)
{
    assert(A && B);
    auto rootA = A->getRoot().lock();
    auto rootB = B->getRoot().lock();
    assert(rootA && rootB);
    
    auto proximity_stack = queryProximity(
            std::move(rootA),std::move(rootB));

    if( !proximity_stack.empty() )
        return true;
    else
        return false;
}

std::stack<NodePair> queryProximity(
        std::shared_ptr<BVH_Node>&& nodeA,
        std::shared_ptr<BVH_Node>&& nodeB)
{
    std::stack<NodePair> qstack;
    qstack.push(std::make_pair(nodeA,nodeB);
    std::stack<NodePair> proximity_stack;

    while( !qstack.empty() )
    {
        auto A = qstack.top().first;
        auto B = qstack.top().second;
        qstack.pop();

        if( A->overlaps(B) )
        {
            if( A->isLeaf() && B->isLeaf() )
            {
                proximity_stack.push(std::make_pair(A,B));
            }
            else if( A->isLeaf() )
            {
                auto rc = B->getRightChild().lock();
                qstack.push(std::make_pair(A,rc));
                auto lc = B->getLeftChild().lock();
                qstack.push(std::make_pair(A,lc));
            }
            else if( B->isLeaf() )
            {
                auto rc = A->getRightChild().lock();
                qstack.push(std::make_pair(rc,B));
                auto lc = A->getLeftChild().lock();
                qstack.push(std::make_pair(lc,B));
            }
            else
            {
                if( A->volume() < B->volume() )
                {
                    auto rc = B->getRightChild().lock();
                    qstack.push(std::make_pair(A,rc));
                    auto lc = B->getLeftChild().lock();
                    qstack.push(std::make_pair(A,lc));
                }
                else
                {
                    auto rc = A->getRightChild().lock();
                    qstack.push(std::make_pair(rc,B));
                    auto lc = A->getLeftChild().lock();
                    qstack.push(std::make_pair(lc,B));
                }
            }
        }
    }

    return proximity_stack;
}



