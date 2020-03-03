#include "newAABB.h"

///////////////////////////////////
//////     AABB methods     //////
/////////////////////////////////

AABB::AABB()
    : lower(3,HUGE), upper(3,-HUGE)
{}

AABB::AABB(const CD_HSE* const h, double TOL)
    : AABB()
{
    tol = TOL;
    mstate = MotionState::STATIC;

    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = h->min_static_coord(i) - tol;
        upperbound[i] = h->max_static_coord(i) + tol;
    }

    /*
    //TODO:
    for (int i = 0; i < 3; ++i)
        indices.push_back(h->Point_of_hse(i)->global_index); 
    */
}

AABB::AABB(const CD_HSE* const h, double TOL, double DT)
    : AABB()
{
    dt = DT;
    tol = TOL;
    mstate = MotionState::MOVING;

    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = h->min_moving_coord(i,dt) - tol;
        upperbound[i] = h->max_moving_coord(i,dt) + tol;
    }

    /*
    //TODO:
    for (int i = 0; i < 3; ++i)
        indices.push_back(h->Point_of_hse(i)->global_index); 
    */
}

//replaces AABB::merge() method
AABB::AABB(const AABB& A, const AABB& B)
    : AABB()
{
    assert(A.volume() > 0 && B.volume() > 0);
    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = std::min(A.lowerbound[i],B.lowerbound[i]);
        upperbound[i] = std::max(A.upperbound[i],B.upperbound[i]);
    }
}

const double AABB::volume() const noexcept
{
    double vol = 1.0;
    for (int i = 0; i < 3; ++i)
        vol *= upperbound[i] - lowerbound[i];
    return vol;
}

const bool AABB::overlaps(const AABB& BB) const
{
    for (int i = 0; i < 3; ++i)
    {
        if (BB.upperbound[i] < lowerbound[i]) return false;
        if (BB.lowerbound[i] > upperbound[i]) return false;
    }
    return true;
}

const bool AABB::contains(const AABB& BB) const
{
    for (int i = 0; i < 3; ++i)
    {
        if (BB.lowerbound[i] <= lowerbound[i]) return false;
        if (BB.upperbound[i] >= upperbound[i]) return false;
    }
    return true;
}

void AABB::encloseHSE(const CD_HSE* const h)
{
    assert(h != nullptr);
    if (mstate == MotionState::STATIC)
    {
        for (int i = 0; i < 3; ++i)
        {
            lowerbound[i] = h->min_static_coord(i) - tol;
            upperbound[i] = h->max_static_coord(i) + tol;
        }
    }
    else
    {
        assert(dt > 0);
        for (int i = 0; i < 3; ++i)
        {
            lowerbound[i] = h->min_moving_coord(i,dt) - tol;
            upperbound[i] = h->max_moving_coord(i,dt) + tol;
        }
    }
}
