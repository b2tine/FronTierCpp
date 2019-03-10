#include "BVH.h"

using NodePair = std::pair<BVH_Node*,BVH_Node*>;

static std::stack<NodePair> queryProximity(
        BVH_Node* nodeA, BVH_Node* nodeB);


const bool checkProximity(BVH* A, BVH* B)
{
    assert(A && B);
    BVH_Node* rootA = A->getRoot();
    BVH_Node* rootB = B->getRoot();
    assert(rootA && rootB);
    
    auto proximity_stack = queryProximity(rootA,rootB);

    if( !proximity_stack.empty() )
        return true;
    else
        return false;
}

std::stack<NodePair> queryProximity(
        BVH_Node* nodeA,
        BVH_Node* nodeB)
{
    std::stack<NodePair> qstack;
    qstack.push(std::make_pair(nodeA,nodeB));
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
                auto rc = B->getRightChild();
                qstack.push(std::make_pair(A,rc));
                auto lc = B->getLeftChild();
                qstack.push(std::make_pair(A,lc));
            }
            else if( B->isLeaf() )
            {
                auto rc = A->getRightChild();
                qstack.push(std::make_pair(rc,B));
                auto lc = A->getLeftChild();
                qstack.push(std::make_pair(lc,B));
            }
            else
            {
                if( A->volume() < B->volume() )
                {
                    auto rc = B->getRightChild();
                    qstack.push(std::make_pair(A,rc));
                    auto lc = B->getLeftChild();
                    qstack.push(std::make_pair(A,lc));
                }
                else
                {
                    auto rc = A->getRightChild();
                    qstack.push(std::make_pair(rc,B));
                    auto lc = A->getLeftChild();
                    qstack.push(std::make_pair(lc,B));
                }
            }
        }
    }

    return proximity_stack;
}

