#ifndef PROXIMITY_H
#define PROXIMITY_H

#include <FronTier.h>
#include "state.h"


void updateAverageVelocityProximity(POINT** pts);

//TODO: implement this
//void updateAverageVelocityCollision(POINT** pts);


class Proximity
{
    public:

        POINT** pts {nullptr};
        double nor[3] {0.0};
        double dist {HUGE};
        
        virtual void computeImpulse() = 0;
        virtual void updateAverageVelocity() = 0;
        
        virtual ~Proximity() = default;
};

class EdgeEdgeProximity : public Proximity
{
    public:

        double a {-1.0};
        double b {-1.0};
        
        EdgeEdgeProximity(POINT** Pts, double* Nor,
                double A, double B, double Dist)
            : pts{Pts}, a{A}, b{B}, dist{Dist}
        {
            for (int i = 0; i < 3; ++i)
                nor[i] = Nor[i];
        }
        
        virtual void computeImpulse() override
        {
            EdgeToEdgeProximityImpulse(pts,nor,a,b,dist);
        }

        virtual void updateAverageVelocity() override
        {
            updateAverageVelocityProximity(pts);
        }
};

class PointTriProximity : public Proximity
{
    public:

        double w[3] {-1.0};

        PointTriProximity(POINT** Pts,
                double* Nor, double* W, double Dist)
            : pts{Pts}, dist{Dist}
        {
            for (int i = 0; i < 3; ++i)
            {
                w[i] = W[i];
                nor[i] = Nor[i];
            }
        }

        void computeImpulse() override
        {
            PointToTriProximityImpulse(pts,nor,w,dist);
        }

        virtual void updateAverageVelocity() override
        {
            updateAverageVelocityProximity(pts);
        }
};

class Collision
{
    public:

        POINT** pts {nullptr};
        double nor[3] {0.0};
        double dist {HUGE};
        double dt {-1.0};
        
        virtual void computeImpulse() = 0;
        virtual void updateAverageVelocity() = 0;
        
        virtual ~Collision() = default;
};

class EdgeEdgeCollision : public Collision
{
    public:

        double a {-1.0};
        double b {-1.0};

        EdgeEdgeCollision(POINT** Pts, double* Nor,
                double A, double B, double Dist, double Dt)
            : pts{Pts}, a{A}, b{B}, dist{Dist}, dt{DT}
        {}

        void computeImpulse() override
        {
            EdgeToEdgeCollisionImpulse(pts,nor,a,b,dist,dt);
        }

        virtual void updateAverageVelocity() override
        {

        }
};

class PointTriCollision : public Collision
{
    public:

        double w[3] {-1.0};

        PointTriCollision(POINT** Pts, double* Nor,
                double* W, double Dist, double Dt)
            : pts{Pts}, dist{Dist}, dt{Dt}
        {
            for (int i = 0; i < 3; ++i)
            {
                w[i] = W[i];
                nor[i] = Nor[i];
            }
        }

        void computeImpulse() override
        {
            PointToTriCollisionImpulse(pts,nor,w,dist,dt);
        }

        virtual void updateAverageVelocity() override
        {

        }
};



#endif
