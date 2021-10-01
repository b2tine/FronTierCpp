#include "FabricManager.h"


void FabricManager::initializeSystem(Front* newfront, double collsn_dt)
{
    front = newfront; //temporary until we figure out how ctor should be written

    collision_solver->initializeSystem(newfront);//will want assembleFromInterface() to be called at FabricManager level
    collision_solver->gpoints = newfront->gpoints;
    collision_solver->gtris = newfront->gtris;
    
    //Overwrite dt set by initializeSystem() with collsn_dt
    CollisionSolver3d::setTimeStepSize(collsn_dt);

    //TODO: Does extra2 get copied when the front is copied? Looks like it ...
    AF_PARAMS* af_params = (AF_PARAMS*)newfront->extra2;
    
    collision_solver->setFabricRoundingTolerance(af_params->fabric_eps);
    collision_solver->setFabricThickness(af_params->fabric_thickness);
    collision_solver->setFabricFrictionConstant(af_params->mu_s);
    collision_solver->setFabricSpringConstant(af_params->ks);
    collision_solver->setFabricPointMass(af_params->m_s);

    collision_solver->setStringRoundingTolerance(af_params->string_eps);
    collision_solver->setStringThickness(af_params->string_thickness);
    collision_solver->setStringFrictionConstant(af_params->mu_l);
    collision_solver->setStringSpringConstant(af_params->kl);
    collision_solver->setStringPointMass(af_params->m_l);

    collision_solver->setStrainLimit(af_params->strain_limit);
    collision_solver->setStrainRateLimit(af_params->strainrate_limit);
    //TODO: add input file options for number of strain limiting iterations

    //For elastic collision impulse control
    collision_solver->setOverlapCoefficient(af_params->overlap_coefficient);

    collision_solver->setVolumeDiff(af_params->vol_diff);
    collision_solver->setRestitutionCoef(1.0);
}

//TODO: Can the action of setCollisionFreePoints() be performed in assembleFromInterface()?
void FabricManager::assembleFromInterface(INTERFACE* intfc)
{
	clearHseList();

	SURFACE** s;
	TRI *tri;
	
    intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    
        surf_tri_loop(*s,tri)
	    {
            CD_HSE_TYPE tag;

            switch (wave_type(*s))
            {
                case ELASTIC_BOUNDARY:
                    tag = CD_HSE_TYPE::FABRIC_TRI;
                    break;
                case NEUMANN_BOUNDARY:
                    tag = CD_HSE_TYPE::STATIC_RIGID_TRI;
                    break;
                case MOVABLE_BODY_BOUNDARY:
                    tag = CD_HSE_TYPE::MOVABLE_RIGID_TRI;
                    break;
                default:
                    printf("assembleFromInterface() ERROR: unknown surface type\n");
                    LOC(); clean_up(EXIT_FAILURE);
            }
            
            hseList.push_back(new CD_TRI(tri,tag));
	    }
	}


	CURVE** c;
	BOND *b;

	intfc_curve_loop(intfc,c)
	{
        if (is_bdry(*c)) continue;
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue; 

        CD_HSE_TYPE tag = CD_HSE_TYPE::STRING_BOND;
	    curve_bond_loop(*c,b)
	    {
            hseList.push_back(new CD_BOND(b,tag));
	    }
	}
}

void FabricManager::clearHseList()
{
    for (unsigned i = 0; i < hseList.size(); ++i)
    {
        delete hseList[i];
    }
    hseList.clear();
}

const std::vector<CD_HSE*>& FabricManager::getHseList() const
{
    return hseList;
}


//TODO: Figure out a way to separate collision handling and strain limiting.
//      Will likely need to move many of the CollisionSolver3d member variables
//      and functions into the FabricManager. The hseLists that are used in
//      collision handling and strain limiting will definitely need to become
//      a member of the FabricManager.
void FabricManager::resolveCollisionSubstep()
{
    collision_solver->resolveCollisionSubstep();
}
