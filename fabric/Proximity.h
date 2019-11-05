#ifndef PROXIMITY_H
#define PROXIMITY_H

#include <FronTier.h>
#include "state.h"

std::unique_ptr<Proximity> checkProximity(const CD_HSE*,const CD_HSE*,double);
std::unique_ptr<Collision> checkCollision(const CD_HSE*,const CD_HSE*,double);


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

class PointTriProximity : public Proximity
{
    public:

        double w[3] {-1.0};


        PointTriProximity(POINT** Pts, double* Nor,
                         double* W, double Dist);


        void computeImpulse() override;
        void updateAverageVelocity() override;
};

class EdgeEdgeProximity : public Proximity
{
    public:

        double a {-1.0};
        double b {-1.0};
        

        EdgeEdgeProximity(POINT** Pts, double* Nor,
                        double A, double B, double Dist);
        

        void computeImpulse() override;
        void updateAverageVelocity() override;
};

class Collision
{
    public:

        POINT** pts {nullptr};
        double nor[3] {0.0};
        double dist {HUGE};

        double dt {-1.0};
        double maxdt {-1.0};
        
        virtual void computeImpulse() = 0;
        virtual void updateState() = 0;
        virtual void restorePrevState() = 0;
        
        virtual ~Collision() = default;
};

class PointTriCollision : public Collision
{
    public:

        double w[3] {-1.0};

        
        PointTriCollision(POINT** Pts, double* Nor,
                double* W, double Dist, double Dt, double MaxDt);


        void computeImpulse() override;
        void updateState() override;
        void restorePrevState() override;
};

class EdgeEdgeCollision : public Collision
{
    public:

        double a {-1.0};
        double b {-1.0};

        EdgeEdgeCollision(POINT** Pts, double* Nor, double A,
                double B, double Dist, double Dt, double MaxDt)
            : pts{Pts}, a{A}, b{B}, dist{Dist}, dt{DT}, maxdt{MaxDt}
        {}

        void computeImpulse() override
        {
            EdgeToEdgeCollisionImpulse(pts,nor,a,b,dist,dt);
        }

        void updateState() override
        {
            UpdateState(pts,dt);
        }

        void restorePrevState() override
        {
            RestorePrevState(pts);
        }
};



#endif
