#ifndef PROXIMITY_H
#define PROXIMITY_H

#include <FronTier.h>
#include "state.h"

std::unique_ptr<Proximity> checkProximity(const CD_HSE*,const CD_HSE*,double);

//std::unique_ptr<Proximity> TriToTri(const TRI*,const TRI*,double);
//std::unique_ptr<Proximity> TriToBond(const TRI*,const BOND*,double);
//std::unique_ptr<Proximity> BondToBond(const BOND*,const BOND*,double);

//std::unique_ptr<Proximity> StaticPointToTri(POINT**,double);
//std::unique_ptr<Proximity> StaticEdgeToEdge(POINT**,double);

std::unique_ptr<Collision> checkCollision(const CD_HSE*,const CD_HSE*,double);

//std::unique_ptr<Collision> MovingTriToTri(const TRI*,const TRI*,double);
//std::unique_ptr<Collision> MovingTriToBond(const TRI*,const BOND*,double);
//std::unique_ptr<Collision> MovingBondToBond(const BOND*,const BOND*,double);

//std::unique_ptr<Collision> KineticPointToTri(POINT**,double,double,double);
//std::unique_ptr<Collision> KineticEdgeToEdge(POINT**,double,double,double);

//Moved to Impulse.h
//
//void PointToTriProximityImpulse(POINT**,double*,double*,double);
//void EdgeToEdgeProximityImpulse(POINT**,double*,double,double,double);
//void PointToTriCollisionImpulse(POINT**,double*,double*,double,double);
//void EdgeToEdgeCollisionImpulse(POINT**,double*,double,double,double,double);


/*
void UpdateAverageVelocity(POINT** pts);
//void updateAverageVelocityProximity(POINT** pts);

void UpdateState(POINT** pts, double dt);
//void updateAverageVelocityCollision(POINT** pts);
void RestorePrevState(POINT** pts);
*/

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
