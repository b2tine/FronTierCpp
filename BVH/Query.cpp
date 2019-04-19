#include "BVH.h"

using NodePair = std::pair<BVH_Node*,BVH_Node*>;

static std::stack<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA, BVH_Node* const nodeB);

static void ProcessProximityCandidates(std::stack<NodePair>& candidates);


//TODO: Return type of this function is temporary for testing.
//      Need to actually perform distance computations for the
//      candidate Hse's corresponding to the nodes/BVs returned
//      in the stack by queryProximity().
const bool checkProximity(const BVH* A, const BVH* B)
{
    assert(A && B);
    BVH_Node* rootA = A->getRoot();
    BVH_Node* rootB = B->getRoot();
    assert(rootA && rootB);
    
    auto proximity_stack = GetProximityCandidates(rootA,rootB);
    ProcessProximityCandidates(proximity_stack);
    //NOTE: at this point the CollisionSolver can apply the
    //      necessary repulsion forces to the nodes remaining
    //      in the stack i.e. the ones that are actually in
    //      proximity of each other.
    //      
    //      Could have the function return the stack, but that could be
    //      an expensive copy/move.
    //
    //      
    //


    /*
    if( !proximity_stack.empty() )
        return true;
    else
        return false;
    */
}

std::stack<NodePair> GetProximityCandidates(
        BVH_Node* const nodeA,
        BVH_Node* const nodeB)
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
                if( A->hasAdjacentHse(B) )
                    continue;
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


void ProcessProximityCandidates(std::stack<NodePair> candidates)
{
    //TODO: Compute distance between each candidate pair,
    //and remove from the stack

}



