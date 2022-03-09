#include "FabricManager.h"


void FabricManager::setCollisionTimeStep(double collsn_dt)
{
    collision_solver->setTimeStepSize(collsn_dt);
}
    
void FabricManager::setCollisionParams(const FABRIC_COLLISION_PARAMS& params)
{
    //TODO: DO WE EVEN NEED A COPY?
    std::copy(&params, &params + 1, &collsn_params);
    
    if (collsn_params.collision_off)
    {
        collision_solver->turnCollision_OFF();
    }

    collision_solver->setFabricRoundingTolerance(collsn_params.fabric_eps);
    collision_solver->setFabricThickness(collsn_params.fabric_thickness);
    collision_solver->setFabricFrictionConstant(collsn_params.mu_s);
    collision_solver->setFabricSpringConstant(collsn_params.k_s);
    collision_solver->setFabricPointMass(collsn_params.m_s);

    collision_solver->setStringRoundingTolerance(collsn_params.string_eps);
    collision_solver->setStringThickness(collsn_params.string_thickness);
    collision_solver->setStringFrictionConstant(collsn_params.mu_l);
    collision_solver->setStringSpringConstant(collsn_params.k_l);
    collision_solver->setStringPointMass(collsn_params.m_l);

    collision_solver->setStrainLimit(collsn_params.strain_limit);
    collision_solver->setCompressiveStrainLimit(collsn_params.compressive_strain_limit);
    collision_solver->setStrainRateLimit(collsn_params.strainrate_limit);
    collision_solver->setStrainVelocityTol(collsn_params.strain_vel_tol);

    collision_solver->setInelasticImpulseCoefficient(collsn_params.inelastic_impulse_coeff);
    collision_solver->setOverlapCoefficient(collsn_params.overlap_coefficient);

    collision_solver->setRestitutionCoef(collsn_params.coefRestitution);
}

void FabricManager::initializeSystem()
{
    assembleHseListFromInterface();
    recordOriginalPositions();
    collision_solver->initializeSystem(hseList);
}

//TODO: Can the action of setCollisionFreePoints() be performed in assembleFromInterface()?
void FabricManager::assembleHseListFromInterface()
{
	clearHseList();

    INTERFACE* intfc = front->interf;
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

void FabricManager::recordOriginalPositions()
{
    std::vector<CD_HSE*>::iterator it;
	for (it = hseList.begin(); it < hseList.end(); ++it)
    {
	    for (int i = 0; i < (*it)->num_pts(); ++i)
        {
            POINT* pt = (*it)->Point_of_hse(i);
            STATE* sl = (STATE*)left_state(pt); 
            
            //NOTE: sl->x_old and sl->x_prevstep for movable rigid body points
            //      are recorded in rgbody_point_propagate()
            if (isMovableRigidBody(pt)) continue;

            for (int j = 0; j < 3; ++j)
            {
                sl->x_prevstep[j] = Coords(pt)[j];
            
                if (std::isnan(sl->x_prevstep[j]))
                {
                    std::cout << "nan_x_old" << std::endl;
                    LOC(); clean_up(ERROR);
                }
            }
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

