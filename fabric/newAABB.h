#ifndef AABB_H
#define AABB_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGAL_Point = Kernel::Point_3;


enum class MotionState {STATIC, MOVING};

using BV_Point = std::vector<double>;


class CD_HSE;

class AABB
{
    public:
        
        BV_Point lowerbound;
        BV_Point upperbound;

        AABB();
        AABB(const CD_HSE* const hse, double pad);
        AABB(const CD_HSE* const hse, double dt, double pad);
        
            //AABB(const CD_HSE* const hse, double TOL);
            //AABB(const CD_HSE* const hse, double TOL, double DT);
    
        AABB(const AABB& A, const AABB& B);
    
        AABB(const AABB&) = default;
        AABB& operator=(const AABB&) = default;
        AABB(AABB&&) = default;
        AABB& operator=(AABB&&) = default;
        ~AABB() = default;
    
        void encloseHSE(const CD_HSE* const hse);
        void encloseMovingHSE(const CD_HSE* const hse, double dt);
            //void encloseHSE(const CD_HSE* const hse, double pad);
            //void encloseMovingHSE(const CD_HSE* const hse, double dt, double pad);
        void expand(double pad);
    
        const double volume() const noexcept;
        const bool overlaps(const AABB& box) const;
        const bool contains(const AABB& box) const;

        const CGAL_Point centroid() const;

    //AABB merge(const AABB&) const;//replaced by ctor taking 2 AABBs
    //bool isCollid(const AABB&);//replaced by overlaps()
    //bool contain(const AABB*);//replaced by contains()

    //protected:
    private:
    
        //double tol;
        //double dt {-1.0};
        //MotionState mstate;
    
    //void updateAABBInfo(double);//replaced by encloseHSE()
                                  //and AABB_Node::update()

    //TODO: Do we really need these?
    //      If so, AABB_Node better suited to hold.
    //
    // indices will store the index of points on the  
    // corresponding triangle or bond.
    //
    //std::vector<long> indices;
    //CD_HSE* hse = nullptr;
};


/*
//dcollid.cpp
bool getProximity(const CD_HSE*,const CD_HSE*);
bool getCollision(const CD_HSE*,const CD_HSE*);
*/

#endif
