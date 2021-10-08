#ifndef FABRIC_MANAGER_H
#define FABRIC_MANAGER_H

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

    FabricManager(Front* fr)
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
    void setCollisionParams(FABRIC_COLLISION_PARAMS params);

    void initializeSystem();

	void clearHseList();
    const std::vector<CD_HSE*>& getHseList() const;
    
    void resolveCollisionSubstep();

private:

    void assembleHseListFromInterface();
    void recordOriginalPositions();
};



#endif
