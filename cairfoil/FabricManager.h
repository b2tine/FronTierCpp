#ifndef FABRIC_MANAGER_H
#define FABRIC_MANAGER_H

#include "airfoil.h"
#include "collid.h"


class FabricManager
{
public: //public temporarily

    Front* front;
    FABRIC_COLLISION_PARAMS collsn_params;
    std::unique_ptr<CollisionSolver3d> collision_solver;
        //CollisionSolver3d* collision_solver;

    int substep;    //use to track which collsn substep we are on
	
    std::vector<CD_HSE*> hseList;

public:

    //TODO: Will need to abandon this Ctor taking Front* front
    explicit FabricManager(Front* fr)
        : front{fr}, 
        collision_solver{std::unique_ptr<CollisionSolver3d>(new CollisionSolver3d)}
    {
        collision_solver->initFront(front);
    }

        /*
        FabricManager()
        {
            collision_solver = new CollisionSolver3d;
        }
        */

    ~FabricManager()
    {
        //delete collision_solver;
        clearHseList();
    }

public:

    void setCollisionTimeStep(double collsn_dt);
    void setCollisionParams(const FABRIC_COLLISION_PARAMS& params);

    void initializeSystem();

	void clearHseList();
    const std::vector<CD_HSE*>& getHseList() const;
    
    void resolveCollisionSubstep();

private:

    //TODO: Write allocation routine for this like in the "first" initialization block in afcnpy.cpp
    SPRING_VERTEX* sv;
    ELASTIC_SET* geom_set;
    GLOBAL_POINT** point_set;
    GLOBAL_POINT* point_set_store;
    GLOBAL_POINT** client_point_set_store;

    void assembleHseListFromInterface();
    void recordOriginalPositions();
};



#endif
