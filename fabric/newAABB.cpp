#include "newAABB.h"

///////////////////////////////////
//////     AABB methods     //////
/////////////////////////////////

AABB::AABB()
    : lowerbound(3,HUGE), upperbound(3,-HUGE)
{}

AABB::AABB(const CD_HSE* const hse, double pad)
    : AABB()
{
    encloseHSE(hse);
    expand(pad);
}

/*
//AABB::AABB(const CD_HSE* const hse, double TOL)
AABB::AABB(const CD_HSE* const hse, double pad)
    : AABB()
{
    //tol = TOL;
    //mstate = MotionState::STATIC;

    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = hse->min_static_coord(i) - pad;
        upperbound[i] = hse->max_static_coord(i) + pad;
    }

    //TODO:
    //for (int i = 0; i < 3; ++i)
      //  indices.push_back(h->Point_of_hse(i)->global_index); 
}
*/

/*
TODO: do we really need this ctor?
//AABB::AABB(const CD_HSE* const hse, double TOL, double DT)
AABB::AABB(const CD_HSE* const hse, double pad, double dt)
    : AABB()
{
    //dt = DT;
    //tol = TOL;
    //mstate = MotionState::MOVING;

    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = hse->min_moving_coord(i,dt) - pad;
        upperbound[i] = hse->max_moving_coord(i,dt) + pad;
    }

    //TODO:
    //for (int i = 0; i < 3; ++i)
      //  indices.push_back(h->Point_of_hse(i)->global_index); 
}
*/

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

const bool AABB::overlaps(const AABB& box) const
{
    for (int i = 0; i < 3; ++i)
    {
        if (box.upperbound[i] < lowerbound[i]) return false;
        if (box.lowerbound[i] > upperbound[i]) return false;
    }
    return true;
}

const bool AABB::contains(const AABB& box) const
{
    for (int i = 0; i < 3; ++i)
    {
        if (box.lowerbound[i] <= lowerbound[i]) return false;
        if (box.upperbound[i] >= upperbound[i]) return false;
    }
    return true;
}

//void AABB::encloseHSE(const CD_HSE* const hse, double pad)
void AABB::encloseHSE(const CD_HSE* const hse)
{
    assert(hse);
    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = hse->min_static_coord(i);
        upperbound[i] = hse->max_static_coord(i);
        //lowerbound[i] = hse->min_static_coord(i) - pad;
        //upperbound[i] = hse->max_static_coord(i) + pad;
    }
}

//void AABB::encloseMovingHSE(const CD_HSE* const hse, double dt, double pad)
void AABB::encloseMovingHSE(const CD_HSE* const hse, double dt)
{
    assert(hse);
    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] = hse->min_moving_coord(i,dt);
        upperbound[i] = hse->max_moving_coord(i,dt);
        //lowerbound[i] = hse->min_moving_coord(i,dt) - pad;
        //upperbound[i] = hse->max_moving_coord(i,dt) + pad;
    }
}

//TODO: see header
/*
void AABB::encloseHSE(const CD_HSE* const hse, double pad)
{
    assert(hse != nullptr);
    if (mstate == MotionState::STATIC)
    {
        for (int i = 0; i < 3; ++i)
        {
            lowerbound[i] = hse->min_static_coord(i) - pad;
            upperbound[i] = hse->max_static_coord(i) + pad;
        }
    }
    else
    {
        assert(dt > 0);
        for (int i = 0; i < 3; ++i)
        {
            lowerbound[i] = hse->min_moving_coord(i,dt) - pad;
            upperbound[i] = hse->max_moving_coord(i,dt) + pad;
        }
    }
}
*/

void AABB::expand(double pad)
{
    assert(pad >= 0);
    for (int i = 0; i < 3; ++i)
    {
        lowerbound[i] -= pad;
        upperbound[i] += pad;
    }
}

