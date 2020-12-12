#include "iFluid.h"
#include "iFturb.h"
#include "keps.h"


KE_PARAMS* Incompress_Solver_Smooth_Basis::computeMuOfKepsModel()
{
    static KE_PARAMS params;
    static KE_CARTESIAN *keps_solver;
    static bool first = true;

    if (first)
    {
        keps_solver = new KE_CARTESIAN(*front);
        keps_solver->read_params(InName(front),&params);
        keps_solver->eqn_params = &params;
        keps_solver->field = NULL;
        keps_solver->initMesh();
        keps_solver->field->vel = iFparams->field->vel;
        keps_solver->field->f_surf = iFparams->field->f_surf;
        keps_solver->eqn_params->mu = iFparams->mu2;
        keps_solver->eqn_params->rho = iFparams->rho2;
        keps_solver->setInitialCondition();
        first = false;
    }

    keps_solver->solve(front->dt);
    return &params;
}

//TODO: Implementation not correct
double Incompress_Solver_Smooth_Basis::computeMuOfBaldwinLomax(
    int *icoords,
	double dist,
	boolean first)
{
	int index;
	COMPONENT comp;
	static double udif;
	static double Fmax;
	double speed, umax = -HUGE, umin = HUGE;
	double mag_vort, wmax = -HUGE;
	double mu_t, mu_out;
	double rho = iFparams->rho2;
	double nu = iFparams->mu2/rho;
	double Fkleb, Fwake;

    //TODO: ymax and Fmax are supposed to be solved for by maximizing the function
    //      F(y) = y*|vorticity|*(1 - exp(-y_plus/A_plus))
	double ymax = iFparams->ymax;

	if (first == YES)
	{
	    first = NO;
	    switch (dim)
	    {
            case 2:
	
            for (int j = jmin; j < jmax; ++j)
	    	for (int i = imin; i < imax; ++i)
	    	{
                index = d_index2d(i,j,top_gmax);
                comp  = cell_center[index].comp;
                if (!ifluid_comp(comp)) continue;

                speed = sqrt(sqr(field->vel[0][index])
                            + sqr(field->vel[1][index]));

                mag_vort = fabs(field->vort[index]);

                umax = std::max(umax,speed);
                umin = std::min(umin,speed);
                wmax = std::max(wmax,mag_vort);
	    	}
            break;

	        case 3:	
        
            for (int k = kmin; k < kmax; ++k)
            for (int j = jmin; j < jmax; ++j)
            for (int i = imin; i < imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                comp  = cell_center[index].comp;
                if (!ifluid_comp(comp)) continue;
                
                speed = sqrt(sqr(field->vel[0][index])
                            + sqr(field->vel[1][index])
                            + sqr(field->vel[2][index]));

                mag_vort = sqrt(sqr(field->vorticity[0][index])
                            + sqr(field->vorticity[1][index])
                            + sqr(field->vorticity[2][index]));

                umax = std::max(umax,speed);
                umin = std::min(umin,speed);
                wmax = std::max(wmax,mag_vort);
            }
            
            break;
        }

	    udif = umax - umin;
	    mag_vort = wmax;
	    Fmax = ymax*mag_vort;
        //TODO: ymax and Fmax are supposed to be solved for by maximizing the function
        //          F(y) = y*|vorticity|*(1 - exp(-y_plus/A_plus))
        //      Above Fmax is not using the damping factor of F(y), 1-exp(..).
	}

	index = d_index(icoords,top_gmax,dim);

    switch (dim)
    {
        case 2:
            mag_vort = fabs(field->vort[index]);
            break;
        case 3:
            mag_vort = sqrt(sqr(field->vorticity[0][index])
                        + sqr(field->vorticity[1][index])
                        + sqr(field->vorticity[2][index]));
            break;
    }

	double l = 0.41*dist;
	double mu_in = rho*l*l*mag_vort; 

	Fwake = std::min(ymax*Fmax,0.25*ymax*sqr(udif)/Fmax);
	Fkleb = 1.0/(1.0 + 5.5*pow((dist*0.3/ymax),6));
	mu_out = rho*0.0168*1.6*Fwake*Fkleb;

    //TODO: This method may select the wrong viscosity.Need to solve for the
    //      crossover distance, y_crx, which is the smallest distance
    //      from the nearest wall in which mu_inner == mu_outer (root finding problem)
    //      and then select mu_t based on whether the distance from the point to the
    //      wall is less/greater than y_crx.
	if (mu_in < mu_out)
	    mu_t = mu_in;
	else
	    mu_t = mu_out;

	return mu_t;
}	/* end computeMuOfBaldwinLomax */

//TODO: 3d runs produce unphysical horizontal bands of eddy viscosity
//      normal to rigid body center of mass motion.
double Incompress_Solver_Smooth_Basis::computeMuofSmagorinskyModel(
                int *icoords)
{
        double delta;
        double S[MAXD][MAXD] = {{0,0,0}, {0,0,0}, {0,0,0}};
        double alpha[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};
        int index0, index[6];

        double **vel = field->vel;
        double C_s = iFparams->C_s;
        
        //TODO: Need to detect boundary's and apply boundary condition (slip, noslip etc.)???
        switch (dim)
	    {
            case 2:
            	
                delta = sqrt(top_h[0]*top_h[1]);
                
                index0 = d_index2d(icoords[0],icoords[1],top_gmax);
                index[0] = d_index2d(icoords[0]-1,icoords[1],top_gmax);
                index[1] = d_index2d(icoords[0]+1,icoords[1],top_gmax);
                index[2] = d_index2d(icoords[0],icoords[1]-1,top_gmax);
                index[3] = d_index2d(icoords[0],icoords[1]+1,top_gmax);
                break;

            case 3:
                
                delta = pow(top_h[0]*top_h[1]*top_h[2], 1.0/3.0);
                
                index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
                index[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax);
                index[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
                index[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
                index[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
                index[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
                index[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
                break;
        }

        for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
        {
            alpha[i][j] = (vel[j][index[2*i+1]] - vel[j][index[2*i]])/(2.0*top_h[i]);
        }

        double sum_S = 0.0;
        for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
        {
            S[i][j] = 0.5*(alpha[i][j] + alpha[j][i]);
            sum_S += S[i][j]*S[i][j];
        }

        double mod_S = sqrt(2.0*sum_S);
        double mu_t = field->rho[index0]*sqr(C_s*delta)*mod_S;
        
        return mu_t;
}       /* end of computeMuofSmagorinskyModel */

//Vreman 2004 paper
double Incompress_Solver_Smooth_Basis::computeMuOfVremanModel(
	int *icoords)
{
    int index = d_index(icoords,top_gmax,dim);
        //int index[6], index0;
        //double alpha[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};
    double beta[MAXD][MAXD] = {{0,0,0}, {0, 0, 0}, {0, 0, 0}};
        //double **vel = field->vel;
        
    double C_v = iFparams->C_v;

    /*
    switch (dim)
    {
        case 2:
            index0 = d_index2d(icoords[0],icoords[1],top_gmax);
            index[0] = d_index2d(icoords[0]-1,icoords[1],top_gmax);
            index[1] = d_index2d(icoords[0]+1,icoords[1],top_gmax);
            index[2] = d_index2d(icoords[0],icoords[1]-1,top_gmax);
            index[3] = d_index2d(icoords[0],icoords[1]+1,top_gmax);
            break;
        
        case 3:
            index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax); 
            index[0] = d_index3d(icoords[0]-1,icoords[1],icoords[2],top_gmax); 
            index[1] = d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            index[2] = d_index3d(icoords[0],icoords[1]-1,icoords[2],top_gmax);
            index[3] = d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            index[4] = d_index3d(icoords[0],icoords[1],icoords[2]-1,top_gmax);
            index[5] = d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            break;
    }

   	double sum_alpha = 0.0;
    for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
	{
	    alpha[i][j] = 0.5*(vel[j][index[2*i+1]] - vel[j][index[2*i]])/top_h[i];
	    sum_alpha += alpha[i][j]*alpha[i][j];
	}
    */

    //Boundary aware computation appears to be working correctly
    auto alpha = computeVelocityGradient(icoords);
   	
    double sum_alpha = 0.0;
    for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
	{
	    sum_alpha += alpha[i][j]*alpha[i][j];
    }

    for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
	{
	    beta[i][j] = 0.0;
	    for (int k = 0; k < dim; ++k)
            beta[i][j] += top_h[k]*top_h[k]*alpha[k][i]*alpha[k][j];
	}

    double B_beta = beta[0][0]*beta[1][1] - beta[0][1]*beta[0][1]
                  + beta[0][0]*beta[2][2] - beta[0][2]*beta[0][2]
                  + beta[1][1]*beta[2][2] - beta[1][2]*beta[1][2];
 
    //see Vreman's implementation, he actually uses 1.0e-12 when checking B_beta
    //  (about 10x larger than MACH_EPS)
    double nu_t;
    if (B_beta < MACH_EPS) //if (sum_alpha >= MACH_EPS)
        nu_t = 0.0;
    else
        nu_t = C_v*sqrt(B_beta/sum_alpha);
    
    double mu_t = nu_t*field->rho[index];

    if (std::isinf(mu_t) || std::isnan(mu_t))
    {
        printf("ERROR: inf/nan eddy viscosity!\n");
        printf("nu_t = %g  rho[%d] = %g\n",nu_t,index,field->rho[index]);
        LOC(); clean_up(EXIT_FAILURE);
    }

    return mu_t;
}	/* end computeMuOfVremanModel*/

std::vector<std::vector<double>>
Incompress_Solver_Smooth_Basis::computeVelocityGradient(
	int *icoords)
{
    double** vel = field->vel;

    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    boolean fr_crx_grid_seg;
    INTERFACE *grid_intfc = front->grid_intfc;
    STATE* intfc_state;
    HYPER_SURF *hs;
    double crx_coords[MAXD];

    std::vector<std::vector<double>> J(dim,std::vector<double>(dim,0.0));

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];
    if (!ifluid_comp(comp)) return J;

    double vel_nb[2];
    double d_h[2];

    for (int l = 0; l < dim; ++l)
    for (int m = 0; m < dim; ++m)
    {
        //l derivatives in m direction
        for (int nb = 0; nb < 2; ++nb)
        {
            d_h[nb] = top_h[m];

            fr_crx_grid_seg = FT_StateStructAtGridCrossing(front,
                    grid_intfc,icoords,dir[m][nb],comp,
                    (POINTER*)&intfc_state,&hs,crx_coords);

            if (!fr_crx_grid_seg)
            {
                int index_nb = next_index_in_dir(icoords,dir[m][nb],dim,top_gmax);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == ELASTIC_BOUNDARY)
            {
                //NOTE: zeroing the relative fluid velocity normal to canopy
                //      via the slip boundary condition will cause the pressure
                //      jump imposed by the ghost fluid method porosity model
                //      to vanish -- destroying the porosity of the canopy
                
                //TODO: Is this the correct approach??
                int index_nb = next_index_in_dir(icoords,dir[m][nb],dim,top_gmax);
                vel_nb[nb] = vel[l][index_nb];
            }
            else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                    wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                double v_slip[MAXD] = {0.0};
                setSlipBoundary(icoords,m,nb,comp,hs,intfc_state,field->vel,v_slip);
                vel_nb[nb] = v_slip[l];
            }
            else if (wave_type(hs) == DIRICHLET_BOUNDARY)
            {
                if (boundary_state_function_name(hs) &&
                        strcmp(boundary_state_function_name(hs),
                            "flowThroughBoundaryState") == 0)
                {
                    //OUTLET
                    //TODO: determine which is correct boundary condition...
                    vel_nb[nb] = intfc_state->vel[l];
                        //vel_nb[nb] = vel[l][index];
                }
                else
                {
                    //INLET
                    vel_nb[nb] = intfc_state->vel[l];
                }
            }
            else
            {
                printf("ERROR: Unknown Boundary Type\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
        }

        J[l][m] = (vel_nb[1] - vel_nb[0])/(d_h[1] + d_h[0]);
    }

    return J;
}

double computeWallShearStress(
        double u_tan,
        double walldist,
        double mu,
        double rho)
{
    if (u_tan < MACH_EPS) return 0.0;
    double u_friction = computeFrictionVelocity(u_tan,walldist,mu,rho);
    double tau_wall = u_friction*u_friction*rho;
    return tau_wall;
}

double computeFrictionVelocity(
        double u_tan,
        double walldist,
        double mu,
        double rho)
{
    double u_wall = computeWallVelocity(u_tan,walldist,mu,rho);
    if (u_wall < MACH_EPS) return 0.0;
    double u_friction = u_tan/u_wall;
    return u_friction;
}

//Computes the dimensionless wall velocity u^{+} = u_tan/u_friction
double computeWallVelocity(
        double u_tan,
        double walldist,
        double mu,
        double rho)
{
    SpaldingWallLaw wallfunc(u_tan,walldist,mu/rho);
    double u_wall_initialguess = u_tan;//TODO: method for better initial guess?
    double u_wall = wallfunc.solve(u_wall_initialguess);
    return u_wall;
}

