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
        virtual void computePostCollisionImpulse(double dt) = 0;
        //virtual void updatePostCollisionState(double dt) = 0;
        
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
        void computePostCollisionImpulse(double dt) override;
        //void updatePostCollisionState(double dt) override;
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
        void computePostCollisionImpulse(double dt) override;
        //void updatePostCollisionState(double dt) override;
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
        virtual void checkNewStateProximity(double tol) = 0;
        virtual void mergeImpactZones() = 0;
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
        void checkNewStateProximity(double tol) override;
        void mergeImpactZones() override;
        void restorePrevState() override;
};

class EdgeEdgeCollision : public Collision
{
    public:

        double a {-1.0};
        double b {-1.0};


        EdgeEdgeCollision(POINT** Pts, double* Nor, double A,
                double B, double Dist, double Dt, double MaxDt);


        void computeImpulse() override;
        void updateState() override;
        void checkNewStateProximity(double tol) override;
        void mergeImpactZones() override;
        void restorePrevState() override;
        }
};



#endif
