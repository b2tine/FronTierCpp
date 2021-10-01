#ifndef FABRIC_MANAGER_H
#define FABRIC_MANAGER_H

#include "collid.h"


class FabricManager
{
public: //public temporarily

    Front* front;
    std::unique_ptr<CollisionSolver3d> collision_solver;
        //CollisionSolver3d* collision_solver;

	std::vector<CD_HSE*> hseList; //will need to move out of CollisionSolver3d

public:

    FabricManager(Front* fr)
        : front{fr}, 
        collision_solver{
            std::unique_ptr<CollisionSolver3d>(new CollisionSolver3d)
        }
    {}

        /*
        FabricManager()
        {
            collision_solver = new CollisionSolver3d;
        }
        */

    ~FabricManager()
    {
        delete collision_solver;
        clearHseList();
    }

public:

    //Temporarily pass in newfront until we have figured out exactly how we want to write ctor
    void initializeSystem(Front* newfront, double collsn_dt);
    void assembleFromInterface(INTERFACE* intfc);

	void clearHseList();
    const std::vector<CD_HSE*>& getHseList() const;
    
    void resolveCollisionSubstep();
};



#endif
