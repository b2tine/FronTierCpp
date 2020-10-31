/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/

/*******************************************************************
 * 			iFcartsn3d.cpp	
 *******************************************************************/
#include "iFluid.h"
#include <solver.h>

//--------------------------------------------------------------------------
// 		   Incompress_Solver_Smooth_3D_Cartesian	
//--------------------------------------------------------------------------

static double (*getStateVel[3])(POINTER) = {getStateXvel,getStateYvel,
                                        getStateZvel};


void Incompress_Solver_Smooth_3D_Cartesian::computeAdvection(void)
{
	int i,j,k,l,index;
	static HYPERB_SOLVER hyperb_solver(*front);
    double speed;

	static double *rho;
	if (rho == NULL)
	{
	    int size = (top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1);
	    FT_VectorMemoryAlloc((POINTER*)&rho,size,sizeof(double));
	}
	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    rho[index] = field->rho[index];
	}
	hyperb_solver.rho = rho;

	hyperb_solver.obst_comp = SOLID_COMP;
	switch (iFparams->adv_order)
	{
	case 1:
	    hyperb_solver.order = 1;
	    hyperb_solver.numericalFlux = upwind_flux;
	    break;
	case 4:
	    hyperb_solver.order = 4;
	    hyperb_solver.numericalFlux = weno5_flux;
	    break;
	default:
	    (void) printf("Advection order %d not implemented!\n",
					iFparams->adv_order);
	    clean_up(ERROR);
	}
    //TODO: probably need to save field->vel before overwriting with soln..
    //      flux soln should be sent to rhs of diffusion computation.... 
	hyperb_solver.dt = m_dt;
	hyperb_solver.var = field->vel;
	hyperb_solver.soln = field->vel;
	hyperb_solver.soln_comp1 = LIQUID_COMP1;
	hyperb_solver.soln_comp2 = LIQUID_COMP2;
	hyperb_solver.rho1 = iFparams->rho1;
	hyperb_solver.rho2 = iFparams->rho2;
	hyperb_solver.getStateVel[0] = getStateXvel;
	hyperb_solver.getStateVel[1] = getStateYvel;
	hyperb_solver.getStateVel[2] = getStateZvel;
	hyperb_solver.findStateAtCrossing = findStateAtCrossing;
	hyperb_solver.solveRungeKutta();

        max_speed = 0.0;
        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            speed = sqrt(sqr(field->vel[0][index]) +
                         sqr(field->vel[1][index]) +
                         sqr(field->vel[2][index]));
            if (max_speed < speed)
            {
                max_speed = speed;
                icrds_max[0] = i;
                icrds_max[1] = j;
                icrds_max[2] = k;
            }
            for (l = 0; l < dim; ++l)
            {
                if (vmin[l] > field->vel[l][index])
                    vmin[l] = field->vel[l][index];
                if (vmax[l] < field->vel[l][index])
                    vmax[l] = field->vel[l][index];
            }
        }
        pp_global_max(&max_speed,1);
        pp_global_min(vmin,dim);
        pp_global_max(vmax,dim);
}

void Incompress_Solver_Smooth_3D_Cartesian::computeNewVelocity(void)
{
	int i, j, k, l, index;
	COMPONENT comp;
	int icoords[MAXD];
	double **vel = field->vel;
	double *phi = field->phi;
	double mag_grad_phi,max_grad_phi,ave_grad_phi;
	double speed;
	int icrds_max[MAXD];

	double grad_phi[MAXD];
    double rho;

	max_grad_phi = ave_grad_phi = 0.0;

	for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = phi[index];
	}
	
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    comp = top_comp[index];
	    if (!ifluid_comp(comp))
	    {
		    for (l = 0; l < 3; ++l)
                vel[l][index] = 0.0;
    		continue;
	    }
	    
        rho = field->rho[index];
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
        
        //Coupling of the velocity impulse due to force of the
        //immersed porous interface accounted for by the addition
        //of GFM pressure jump source terms into grad_phi 
        computeFieldPointGradJump(icoords,phi,grad_phi);

        speed = 0.0;
	    for (l = 0; l < 3; ++l)
	    {
	    	vel[l][index] -= accum_dt/rho*grad_phi[l];
		    speed += sqr(vel[l][index]);
	    }
        speed = sqrt(speed);

	    mag_grad_phi = Mag3d(grad_phi);
	    ave_grad_phi += mag_grad_phi;

        if (mag_grad_phi > max_grad_phi)
	    {
            max_grad_phi = mag_grad_phi;	
            icrds_max[0] = i;
            icrds_max[1] = j;
            icrds_max[2] = k;
	    }

        if (speed > iFparams->ub_speed)
	    {
	    	for (l = 0; l < 3; ++l)
                vel[l][index] *= iFparams->ub_speed/speed;
            speed = iFparams->ub_speed;
	    }
	}

	FT_ParallelExchGridVectorArrayBuffer(vel,front);
	
    if (debugging("step_size"))
	{
	    (void) printf("Max gradient phi = %f  occuring at: %d %d %d\n",
			max_grad_phi,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf("Ave gradient phi = %f\n",ave_grad_phi/(imax-imin+1)
			/(jmax-jmin+1)/(kmax-kmin+1));
	}
	if (debugging("check_div"))
	{
	    checkVelocityDiv("After computeNewVelocity3d()");
	}
}	/* end computeNewVelocity3d */

//TODO: move these to iFbasic.cpp
void Incompress_Solver_Smooth_Basis::printEnstrophy()
{
    switch (dim)
    {
        case 2:
            printEnstrophy2d();
            break;
        case 3:
            printEnstrophy3d();
            break;
        default:
            printf("ERROR: Dimension must be 2 or 3\n");
            LOC(); clean_up(EXIT_FAILURE);
    }
}

void Incompress_Solver_Smooth_Basis::printEnstrophy2d()
{
    static bool first = true;
    static FILE* efile;
	static char fname[512];

    if (first)
    {
	    sprintf(fname,"%s/enstrophy.xg",OutName(front));
        efile = fopen(fname,"w");
        first = false;
    }
    else
    {
        efile = fopen(fname,"a");
    }

    int index;
    double enstrophy = 0.0;
    double* vorticity = field->vort;
    
    for (int i = imin; i < imax; ++i)
    for (int j = jmin; j < jmax; ++j)
    {
        index = d_index2d(i,j,top_gmax);
        if (!ifluid_comp(top_comp[index])) continue;
        
        enstrophy += sqr(vorticity[index]);
    }
    
    double vol_elem = top_h[0]*top_h[1];
    enstrophy *= vol_elem;

    fprintf(efile,"%g %g\n",front->time,enstrophy);
    fclose(efile);
}

void Incompress_Solver_Smooth_Basis::printEnstrophy3d()
{
    static bool first = true;
    static FILE* efile;
	static char fname[512];

    if (first)
    {
	    sprintf(fname,"%s/enstrophy.xg",OutName(front));
        efile = fopen(fname,"w");
        first = false;
    }
    else
    {
        efile = fopen(fname,"a");
    }

    int index;
    double enstrophy = 0.0;
    double** vorticity = field->vorticity;
    
    for (int i = imin; i < imax; ++i)
    for (int j = jmin; j < jmax; ++j)
    for (int k = kmin; k < kmax; ++k)
    {
        index = d_index3d(i,j,k,top_gmax);
        if (!ifluid_comp(top_comp[index])) continue;
        
        double sqrmag_vort = 0.0;
        for (int l = 0; l < dim; ++l)
            sqrmag_vort += sqr(vorticity[l][index]);

        enstrophy += sqrmag_vort;
    }
    
    double vol_elem = top_h[0]*top_h[1]*top_h[2];
    enstrophy *= vol_elem;

    fprintf(efile,"%g %g\n",front->time,enstrophy);
    fclose(efile);
}

/*
void Incompress_Solver_Smooth_3D_Cartesian::computeVorticity()
{
    int index;
    int icoords[MAXD];

    int index_xnb0, index_xnb1;
    int icnb_x0[MAXD], icnb_x1[MAXD];

    int index_ynb0, index_ynb1;
    int icnb_y0[MAXD], icnb_y1[MAXD];

    int index_znb0, index_znb1;
    int icnb_z0[MAXD], icnb_z1[MAXD];
	
    double **vel = field->vel;
	double **vorticity = field->vorticity;

	for (int k = kmin; k <= kmax; k++)
	for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
        if (!ifluid_comp(top_comp[index]))
        {
            for (int l = 0; l < 3; ++l)
                vorticity[l][index] = 0.0;
            continue;
        }

        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;

        for (int l = 0; l < dim; ++l)
        {
            icnb_x0[l] = icoords[l];
            icnb_x1[l] = icoords[l];
            icnb_y0[l] = icoords[l];
            icnb_y1[l] = icoords[l];
            icnb_z0[l] = icoords[l];
            icnb_z1[l] = icoords[l];
        }

        //cell centered derivative indices wrt x
        icnb_x0[0] = icoords[0] - 1;
        icnb_x1[0] = icoords[0] + 1;
        index_xnb0 = d_index(icnb_x0,top_gmax,dim);
        index_xnb1 = d_index(icnb_x1,top_gmax,dim);

        //cell centered derivative indices wrt y
        icnb_y0[1] = icoords[1] - 1;
        icnb_y1[1] = icoords[1] + 1;
        index_ynb0 = d_index(icnb_y0,top_gmax,dim);
        index_ynb1 = d_index(icnb_y1,top_gmax,dim);

        //cell centered derivative indices wrt z
        icnb_z0[2] = icoords[2] - 1;
        icnb_z1[2] = icoords[2] + 1;
        index_znb0 = d_index(icnb_z0,top_gmax,dim);
        index_znb1 = d_index(icnb_z1,top_gmax,dim);


        //x component vorticity
        double u2_wrty = 0.5*(vel[2][index_ynb1] - vel[2][index_ynb0])/top_h[1];
        double u1_wrtz = 0.5*(vel[1][index_znb1] - vel[1][index_znb0])/top_h[2];
        vorticity[0][index] = u2_wrty - u1_wrtz;

        //y component vorticity
        double u0_wrtz = 0.5*(vel[0][index_znb1] - vel[0][index_znb0])/top_h[2];
        double u2_wrtx = 0.5*(vel[2][index_xnb1] - vel[2][index_xnb0])/top_h[0];
        vorticity[1][index] = u0_wrtz - u2_wrtx;

        //z component vorticity
        double u1_wrtx = 0.5*(vel[1][index_xnb1] - vel[1][index_xnb0])/top_h[0];
        double u0_wrty = 0.5*(vel[0][index_ynb1] - vel[0][index_ynb0])/top_h[1];
        vorticity[2][index] = u1_wrtx - u0_wrty;
    }
	
    //TODO: where/how to do this? See cFluid copyMeshStates() ...
    //FT_ParallelExchGridVectorArrayBuffer(vorticity,front);
}
*/

void Incompress_Solver_Smooth_3D_Cartesian::computeVorticity()
{
    double** vel = field->vel;
    double** vorticity = field->vorticity;
    std::vector<double> curl_vel(3);

    int index;
    int icoords[MAXD];

    for (int k = kmin; k <= kmax; ++k)
    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        index = d_index3d(i,j,k,top_gmax);
        if (!ifluid_comp(top_comp[index])) continue;
        
        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;
        
        curl_vel = computePointVorticity(icoords,vel);

        vorticity[0][index] = curl_vel[0];
        vorticity[1][index] = curl_vel[1];
        vorticity[2][index] = curl_vel[2];
    }
	
    FT_ParallelExchGridVectorArrayBuffer(vorticity,front);
}

//TODO: Turn into global curl function
std::vector<double> Incompress_Solver_Smooth_3D_Cartesian::
    computePointVorticity(int* icoords, double** vel)
{
    HYPER_SURF *hs;
    POINTER intfc_state;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    double crx_coords[MAXD];
    int status;

    int index = d_index(icoords,top_gmax,dim);
    COMPONENT comp = top_comp[index];
    
    int icnb[MAXD], i_index_nb;
    int iicnb[MAXD], ii_index_nb;

    double vel_i[2];
    double vel_ii[2];

    std::vector<double> curl_vel(dim);
    for (int idir = 0; idir < 3; ++idir)
    {
        for (int k = 0; k < 3; ++k)
        {
            icnb[k] = icoords[k];
            iicnb[k] = icoords[k];
        }

        int i = (idir+1)%3;
        int ii = (idir+2)%3;

        for (int nb = 0; nb < 2; nb++)
        {
            icnb[i] = (nb == 0) ? icoords[i] - 1 : icoords[i] + 1;
            iicnb[ii] = (nb == 0) ? icoords[ii] - 1 : icoords[ii] + 1;
            
            i_index_nb = d_index(icnb,top_gmax,dim);
            ii_index_nb = d_index(iicnb,top_gmax,dim);

            //differentiate vel[ii=(idir+2)%3] with respect to coordinate i = (idir+1)%3;
            vel_ii[nb] = vel[ii][i_index_nb];

            //differentiate vel[i=(idir+1)%3] with respect to coordinate ii = (idir+2)%3;
            vel_i[nb] = vel[i][ii_index_nb];
        }

        double uii_wrt_xi = 0.5*(vel_ii[1] - vel_ii[0])/top_h[i];
        double ui_wrt_xii = 0.5*(vel_i[1] - vel_i[0])/top_h[ii];
        curl_vel[idir] = uii_wrt_xi - ui_wrt_xii;
    }

    return curl_vel;
}       /* end computePointVorticity */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeSourceTerm(double *coords, double *source) 
{
    if(iFparams->if_buoyancy)
    {
        int ic[MAXD],index;
        rect_in_which(coords,ic,top_grid);
        index = d_index(ic,top_gmax,dim);
        for (int i = 0; i < dim; ++i)
        {
            source[i] = field->ext_accel[i][index];
        }
    }
    else
    {
        for (int i = 0; i < dim; ++i)
            source[i] = iFparams->gravity[i];
    }
} 	/* computeSourceTerm */

#include<fstream>
void Incompress_Solver_Smooth_3D_Cartesian::solve(double dt)
{
	static boolean first = YES;
	if (first)
	{
	    accum_dt = 0.0;
	    first = NO;
	}
	m_dt = dt;

    if (debugging("step_size"))
    {
        computeMaxSpeed();
        printf("max speed entering solve(): %20.14f\n",max_speed);
    }

	start_clock("solve");
	setDomain();

    //TODO: Put into a debugging string block
    /*
    ///////////////////////////////////////
    ////     For grid debugging       ////
    //////////////////////////////////////
    printf("\ntopological grid:\n");
    print_RECT_GRID_structure(top_grid);
    printf("\ncomputational grid:\n");
    RECT_GRID *rgr = computational_grid(front->interf);
    print_RECT_GRID_structure(rgr);
    auto coords0 = cell_center[0].getCoords();
    printf("\ncell_center[0].m_coords = %f %f %f\n",
            coords0[0],coords0[1],coords0[2]);
    auto coords1 = cell_center[1].getCoords();
    printf("\ncell_center[1].m_coords = %f %f %f\n",
            coords1[0],coords1[1],coords1[2]);
    clean_up(0);
    ///////////////////////////////////////
    */

	setComponent();
	if (debugging("trace"))
	    printf("Passed setComponent()\n");

	paintAllGridPoint(TO_SOLVE);
	setGlobalIndex();
	if (debugging("trace"))
	    printf("Passed setGlobalIndex()\n");

	start_clock("setSmoothedProperties");
	setSmoothedProperties();
	stop_clock("setSmoothedProperties");
	if (debugging("trace"))
	    printf("Passed setSmoothedProperties()\n");
    
    addImmersedForce();
	
	appendOpenEndStates();
	// 1) solve for intermediate velocity
	start_clock("computeAdvection");
	computeAdvection();
	stop_clock("computeAdvection");
	if (debugging("check_div") || debugging("step_size"))
	{
	    computeMaxSpeed();
	    checkVelocityDiv("After computeAdvection()");
	    (void) printf("max_speed after computeAdvection(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
	if (debugging("sample_velocity"))
	    sampleVelocity();
	
	start_clock("computeDiffusion");
	computeDiffusion();
	stop_clock("computeDiffusion");

	if (debugging("step_size"))
	{
	    computeMaxSpeed();
	    checkVelocityDiv("After computeDiffusion()");
	    (void) printf("max_speed after computeDiffusion(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
	if (debugging("sample_velocity"))
	    sampleVelocity();

	// 2) projection step
	accum_dt += m_dt;
	if (accum_dt >= min_dt)
	{
	    start_clock("computeProjection");
	    computeProjection();
	    stop_clock("computeProjection");

	    start_clock("computePressure");
	    computePressure();
	    stop_clock("computePressure");
	    if (debugging("trace"))
		printf("min_pressure = %f  max_pressure = %f\n",
			min_pressure,max_pressure);

        //TODO: appendOpenEndStates() appears to be deprecated,
        //      and the OPEN_BOUNDARY condition replaced by the
        //      FLOW_THROUGH_BOUNDARY condition.
        
        //TODO: Is this getting taken care of??
        //      If not, need a new function to update the values
        //      of pressure and phi at the boundary.
        appendOpenEndStates(); //necessary since phi is updated
	    
        start_clock("computeNewVelocity");
	    computeNewVelocity();
	    stop_clock("computeNewVelocity");

        accum_dt = 0.0;
	}
	
    computeMaxSpeed();
    computeVorticity();

	if (debugging("step_size"))
	{
	    printf("max_speed after computeNewVelocity(): %20.14f\n",
				max_speed);
	    (void) printf("max speed occured at (%d %d %d)\n",icrds_max[0],
				icrds_max[1],icrds_max[2]);
	}
	if (debugging("sample_velocity"))
	    sampleVelocity();
	start_clock("copyMeshStates");
	copyMeshStates();
	stop_clock("copyMeshStates");

        if(iFparams->if_ref_pres == YES)
            setReferencePressure();

	setAdvectionDt();

    if (debugging("step_size"))
    {
        computeMaxSpeed();
        printf("max speed leaving solve(): %20.14f\n",max_speed);
    }

	stop_clock("solve");
}	/* end solve */

void Incompress_Solver_Smooth_3D_Cartesian::copyMeshStates()
{
	int i,j,k,d,index;
	double *pres = field->pres;
    
    //TODO: what about velocity?

	for (i = imin; i <= imax; ++i)
	for (j = jmin; j <= jmax; ++j)
	for (k = kmin; k <= kmax; ++k)
	{
	    index  = d_index3d(i,j,k,top_gmax);
        if (!ifluid_comp(top_comp[index]))
	    	pres[index] = 0.0;
	}
    //TODO: compare to cFluid G_CARTESIAN::copyMeshStates()
    FT_ParallelExchGridArrayBuffer(pres,front,NULL);
}	/* end copyMeshStates */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusion(void)
{
	return computeDiffusionCN();
	//return computeDiffusionImplicit();
}	/* end computeDiffusion */

/*
void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionCN(void)
{
    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    double **grad_phi = field->grad_phi;
    double **vel = field->vel;
	double **f_surf = field->f_surf;
    double *mu = field->mu;
    double *rho = field->rho;
	
    POINTER intfc_state;
	HYPER_SURF *hs;
    
    COMPONENT comp;
    int index,index_nb;
    int I,I_nb;
    int crx_status;
    boolean status;

	int icoords[MAXD], icnb[MAXD];
	double crx_coords[MAXD];
    double nor[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionCN()\n");

    setIndexMap();
    int size = iupper - ilower;
    
    double *x;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	double source[MAXD];
    for (int l = 0; l < dim; ++l)
        source[l] = iFparams->gravity[l];

    //TODO: Save last step solns and provide as initial guess.
    //      Pass in as argument like in ELLIPTIC_SOLVER::solve3d()?
    for (int l = 0; l < 3; ++l)
    {
        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
        solver.Reset_A();
        solver.Reset_b();
        solver.Reset_x();
        boolean use_neumann_solver = YES;

        for (int k = kmin; k <= kmax; k++)
        for (int j = jmin; j <= jmax; j++)
        for (int i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;
            
            index = d_index3d(i,j,k,top_gmax);
            comp = top_comp[index];
            
            if (!ifluid_comp(comp))
            {
                for (int m = 0; m < 3; ++m)
                    vel[m][index] = 0.0;
                continue;
            }
            
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            
            double aII = 1.0;
            double coeff_rhs = 1.0;
            double RHS = 0.0;

            double nu_index = mu[index]/rho[index];
            
            for (int idir = 0; idir < dim; ++idir)
            {
                for (int m = 0; m < dim; ++m)
                    icnb[m] = icoords[m];
                    
                double lambda = 0.5*m_dt/sqr(top_h[idir]);

                for (int nb = 0; nb < 2; ++nb)
                {
                    icnb[idir] = (nb == 0) ?
                        icoords[idir] - 1 : icoords[idir] + 1;

                    index_nb  = d_index(icnb,top_gmax,dim);
                    I_nb = ijk_to_I[icnb[0]][icnb[1]][icnb[2]];

                    double coeff_nb = 0.0;
                    double nu_halfidx = 0.5*nu_index;
                    
                    crx_status = (*findStateAtCrossing)(front,icoords,
                            dir[idir][nb],comp,&intfc_state,&hs,crx_coords);

                    if (crx_status)
                    {
                        if (wave_type(hs) == DIRICHLET_BOUNDARY)
                        {
                            //INFLOW/OUTFLOW BOUNDARY ONLY, NOT NOSLIP WALL BOUNDARY!
                            nu_halfidx += 0.5*getStateMu(intfc_state)/rho[index_nb];
                            coeff_nb = -1.0*lambda*nu_halfidx;
                            aII -= coeff_nb;
                            coeff_rhs += coeff_nb;
                            
                            double bval = getStateVel[idir](intfc_state);
                            bval += m_dt*field->grad_phi[idir][index_nb];
                            RHS -= 2.0*coeff_nb*bval;
                            use_neumann_solver = NO;
                            
                            //if (iFparams->num_scheme.projc_method == SIMPLE ||
                            //    iFparams->num_scheme.projc_method == KIM_MOIN)
                            //    bval += m_dt*field->grad_phi[idir][index_nb];
                            }//
                        }
                        else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                                 wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                        {
                            status = FT_NormalAtGridCrossing(front,icoords,
                                    dir[idir][nb],comp,nor,&hs,crx_coords);

                  //          //
                  //          //ghost point
                  //          double coords_ghost[MAXD];
                  //          getRectangleCenter(index_nb,coords_ghost);
                  //          
                  //          //Reflect the ghost point through intfc-mirror at crossing.
                  //          //first reflect across the grid line containing intfc crossing.
                  //          double coords_reflect[MAXD];
                  //          for (int m = 0; m < dim; ++m)
                  //              coords_reflect[m] = coords_ghost[m];
                  //          coords_reflect[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];
                  //          //(^should just be the coords at current index)
                  //          //

                            //Reflect the ghost point through intfc-mirror at crossing.
                            //
                            //first reflect across the grid line containing intfc crossing.
                            //Should be the coords of the current index
                            double coords_reflect[MAXD];
                            getRectangleCenter(index,coords_reflect);

                            //Reflect the displacement vector across the line
                            //containing the intfc normal vector
                            double v[MAXD];
                            double vn = 0.0;

                            for (int m = 0; m < dim; ++m)
                            {
                                v[m] =  coords_reflect[m] - crx_coords[m];
                                vn += v[m]*nor[m];
                            }

                            for (int m = 0; m < dim; ++m)
                                v[m] = 2.0*vn*nor[m] - v[m];

                            //The desired reflected point
                            for (int m = 0; m < dim; ++m)
                                coords_reflect[m] = crx_coords[m] + v[m];

                            //Interpolate the velocity at the reflected point
                            double vel_reflect[MAXD];
                            for (int m = 0; m < dim; ++m)
                            {
                                FT_IntrpStateVarAtCoords(front,comp,
                                        coords_reflect,vel[m],getStateVel[m],
                                        &vel_reflect[m],&vel[m][index]);
                            }

                            //Ghost vel has relative normal velocity component equal
                            //in magnitude to reflected point's relative normal velocity
                            //and opposite in direction.
                            vn = 0.0;
                            double vel_rel[MAXD];
                            double* vel_intfc = ((STATE*)intfc_state)->vel;
                            for (int m = 0; m < dim; ++m)
                            {
                                vel_rel[m] = vel_reflect[m] - vel_intfc[m];
                                vn += vel_rel[m]*nor[m];
                            }

                            double vel_ghost[MAXD];
                            for (int m = 0; m < dim; ++m)
                                vel_ghost[m] = vel_reflect[m] - 2.0*vn*nor[m];
                        
                            //nu_halfidx += 0.0;
                            coeff_nb = -1.0*lambda*nu_halfidx;
                                //solver.Set_A(I,I_nb,coeff_nb);
                            aII -= coeff_nb;
                            coeff_rhs += coeff_nb;
                                //RHS -= coeff_nb*vel[l][index_nb];
                            RHS -= 2.0*coeff_nb*vel_ghost[l];
                        }
                    }
                    else
                    {
                        //NO_PDE_BOUNDARY
                        nu_halfidx += 0.5*mu[index_nb]/rho[index_nb];
                        coeff_nb = -1.0*lambda*nu_halfidx;
                        solver.Set_A(I,I_nb,coeff_nb);
                        aII -= coeff_nb;
                        coeff_rhs += coeff_nb;
                        RHS -= coeff_nb*vel[l][index_nb];
                    }
                }

            }

            //TODO: Was this not working because the solver isn't in the correct loop block????
            //      Or did a bracket get deleted when commenting out everything?
            solver.Set_A(I,I,aII);

            RHS += coeff_rhs*vel[l][index];
		    RHS += m_dt*source[l];
		    RHS += m_dt*f_surf[l][index];
            solver.Set_b(I,RHS);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-10);

        use_neumann_solver = pp_min_status(use_neumann_solver);
        bool Try_GMRES = false;
	
        PetscInt num_iter;
	    double rel_residual;

	    start_clock("Before Petsc solve");
        if (use_neumann_solver)
        {
            //
            //if (skip_neumann_solver)
            //{
            //    //TODO: is this needed??
            //}
            ///

            printf("\ncomputeDiffusionCN(): Using Neumann Solver!\n");
            solver.Solve_withPureNeumann();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if (rel_residual > 1)
            {
                printf("\n The solution diverges! The residual \
                   is %g. Solve again using GMRES!\n",rel_residual);
                Try_GMRES = true;
            }
        }
        else
        {
            printf("\ncomputeDiffusionCN(): Using non-Neumann Solver!\n");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

            if (rel_residual > 1)
            {
                printf("\n The solution diverges! The residual \
                   is %g. Solve again using GMRES!\n",rel_residual);
                Try_GMRES = true;
            }
        }

        if (Try_GMRES)
        {
            solver.Reset_x();
            solver.Solve_GMRES();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);
            
            if(rel_residual > 1)
            {
                printf("\n The solution diverges using GMRES! \
                        The residual is %g. Exiting ...\n",rel_residual);
                clean_up(EXIT_FAILURE);
            }
        }

	    stop_clock("After Petsc solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            printf("L_CARTESIAN::computeDiffusionCN(): \
                    num_iter = %d, rel_residual = %g. \n",
                    num_iter,rel_residual);
        }

	    for (int k = kmin; k <= kmax; k++)
        for (int j = jmin; j <= jmax; j++)
        for (int i = imin; i <= imax; i++)
        {
            //TODO: write to another soln array instead
            //      of immediately overwriting the vel array.
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
                vel[l][index] = x[I-ilower];
            else
                vel[l][index] = 0.0;
        }

    }
	
    FT_ParallelExchGridVectorArrayBuffer(vel,front);
    FT_FreeThese(1,x);

	if (debugging("trace"))
    {
        printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::\
                computeDiffusionCN()\n");
    }
}*/       /* end computeDiffusionCN */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionCN(void)
{
    COMPONENT comp;
    int index,index_nb[6],size;
    int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
    double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;
    int status;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionCN()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
        PETSc solver;
        solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I  = ijk_to_I[i][j][k];
            if (I == -1) continue;

            index  = d_index3d(i,j,k,top_gmax);
            index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            index_nb[5] = d_index3d(i,j,k+1,top_gmax);

            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            comp = top_comp[index];

            I_nb[0] = ijk_to_I[i-1][j][k]; //west
            I_nb[1] = ijk_to_I[i+1][j][k]; //east
            I_nb[2] = ijk_to_I[i][j-1][k]; //south
            I_nb[3] = ijk_to_I[i][j+1][k]; //north
            I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            mu0 = field->mu[index];
            rho = field->rho[index];

            for (nb = 0; nb < 6; nb++)
            {
                if ((*findStateAtCrossing)(front,icoords,dir[nb],comp,
                                    &intfc_state,&hs,crx_coords))
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            U_nb[nb] = vel[l][index];
                        }
                        else
                        {
                            //INLET
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        if (!is_bdry_hs(hs))
                        {
                            //TODO: shouldn't use slip boundary until turb model is activated
                            double v_slip[MAXD] = {0.0};
                            int idir = nb/2; int nbr = nb%2; //quick hack to avoid restructuring loop while prototyping
                            setSlipBoundary(icoords,idir,nbr,comp,hs,intfc_state,field->vel,v_slip);
                            U_nb[nb] = v_slip[l];
                        }
                        else
                        {
                            //TODO: Without this rayleigh-taylor with NEUMANN boundaries
                            //      crashes for some reasone
                            U_nb[nb] = getStateVel[l](intfc_state);
                        }
                    }
                    else
                    {
                        printf("Unkown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }

                    if (wave_type(hs) == DIRICHLET_BOUNDARY || neumann_type_bdry(wave_type(hs)))
                        mu[nb] = mu0;
                    else
                        mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                
                }
                else
                {
                    U_nb[nb] = vel[l][index_nb[nb]];
                    mu[nb] = 0.5*(mu0 + field->mu[index_nb[nb]]);
                }
                        
            }
            
            coeff[0] = 0.5*m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            coeff[1] = 0.5*m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            coeff[2] = 0.5*m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            coeff[3] = 0.5*m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            coeff[4] = 0.5*m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            coeff[5] = 0.5*m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            getRectangleCenter(index, coords);
            computeSourceTerm(coords, source);

            aII = 1.0+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
            rhs = (1.0-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*(vel[l][index]);

            for(nb = 0; nb < 6; nb++)
            {
                status = (*findStateAtCrossing)(front,icoords,dir[nb],comp,
                        &intfc_state,&hs,crx_coords);
                
                if (status == NO_PDE_BOUNDARY)
                {
                    solver.Set_A(I,I_nb[nb],-coeff[nb]);
                    rhs += coeff[nb]*U_nb[nb];
                }
                else
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (boundary_state_function(hs) &&
                            strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                        {
                            //OUTLET
                            aII -= coeff[nb];
                            rhs += coeff[nb]*U_nb[nb];
                        }
                        else
                        {
                            //INLET
                            rhs += 2.0*coeff[nb]*U_nb[nb];
                        }
                    }
                    else if (neumann_type_bdry(wave_type(hs)))
                    {
                        //NEUMANN
                        rhs += 2.0*coeff[nb]*U_nb[nb];
                    }
                    else
                    {
                        printf("Unkown Boundary Type!\n");
                        LOC(); clean_up(EXIT_FAILURE);
                    }
                }
            }

            rhs += m_dt*source[l];
            rhs += m_dt*f_surf[l][index];

                //rhs -= m_dt*grad_q[l][index]/rho;
            
            solver.Set_A(I,I,aII);
            solver.Set_b(I,rhs);
        }

        solver.SetMaxIter(40000);
        solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
        solver.Solve();
        solver.GetNumIterations(&num_iter);
        solver.GetResidualNorm(&residual);

	    stop_clock("After Petsc solve");

        // get back the solution
        solver.Get_x(x);

        if (debugging("PETSc"))
        {
            (void) printf("L_CARTESIAN::"
                    "computeDiffusionCN: "
                    "num_iter = %d, residual = %g. \n",
                    num_iter,residual);
        }

	    for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            I = ijk_to_I[i][j][k];
            index = d_index3d(i,j,k,top_gmax);
            if (I >= 0)
                vel[l][index] = x[I-ilower];
            else
                vel[l][index] = 0.0;
        }

    }

	FT_ParallelExchGridVectorArrayBuffer(vel,front);
    FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionCN()\n");
}       /* end computeDiffusionCN */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionExplicit(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionExplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            for (k = kmin; k <= kmax; k++)
	    {
            	index = d_index3d(13,13,k,top_gmax);
	    }
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index  = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	mu0   = field->mu[index];
            	rho   = field->rho[index];

            	for (nb = 0; nb < 6; nb++)
            	{
                    if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			{
                    	    U_nb[nb] = vel[l][index];
			}
			else
			{
			    U_nb[nb] = getStateVel[l](intfc_state);
			}
			if (wave_type(hs) == DIRICHLET_BOUNDARY || 
			    neumann_type_bdry(wave_type(hs)))
			    mu[nb] = mu0;
			else
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
                    else
		    {
                    	U_nb[nb] = vel[l][index_nb[nb]];
			mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
            	}

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		rhs = (-coeff[0]-coeff[1]-coeff[2]-coeff[3]-coeff[4]-coeff[5])*
		      		vel[l][index];

		int num_nb = 0;
		for(nb = 0; nb < 6; nb++)
		{
		    rhs += coeff[nb]*U_nb[nb];
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
                x[I-ilower] = vel[l][index] + rhs;
            }

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionExplicit()\n");
}       /* end computeDiffusionExplicit */

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionImplicit(void)
{
        COMPONENT comp;
        int index,index_nb[6],size;
        int I,I_nb[6];
	int i,j,k,l,nb,icoords[MAXD];
	double coords[MAXD], crx_coords[MAXD];
	double coeff[6],mu[6],mu0,rho,rhs,U_nb[6];
        double *x;
	GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	POINTER intfc_state;
	HYPER_SURF *hs;
	PetscInt num_iter;
	double rel_residual;
	double aII;
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	INTERFACE *grid_intfc = front->grid_intfc;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionImplicit()\n");

        setIndexMap();

        size = iupper - ilower;
        FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

	for (l = 0; l < dim; ++l)
	{
            PETSc solver;
            solver.Create(ilower, iupper-1, 7, 7);
	    solver.Reset_A();
	    solver.Reset_b();
	    solver.Reset_x();

            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
            	I  = ijk_to_I[i][j][k];
            	if (I == -1) continue;

            	index  = d_index3d(i,j,k,top_gmax);
            	index_nb[0] = d_index3d(i-1,j,k,top_gmax);
            	index_nb[1] = d_index3d(i+1,j,k,top_gmax);
            	index_nb[2] = d_index3d(i,j-1,k,top_gmax);
            	index_nb[3] = d_index3d(i,j+1,k,top_gmax);
            	index_nb[4] = d_index3d(i,j,k-1,top_gmax);
            	index_nb[5] = d_index3d(i,j,k+1,top_gmax);

		icoords[0] = i;
		icoords[1] = j;
		icoords[2] = k;
		comp = top_comp[index];

            	I_nb[0] = ijk_to_I[i-1][j][k]; //west
            	I_nb[1] = ijk_to_I[i+1][j][k]; //east
            	I_nb[2] = ijk_to_I[i][j-1][k]; //south
            	I_nb[3] = ijk_to_I[i][j+1][k]; //north
            	I_nb[4] = ijk_to_I[i][j][k-1]; //lower
            	I_nb[5] = ijk_to_I[i][j][k+1]; //upper


            	mu0   = field->mu[index];
            	rho   = field->rho[index];

        for (nb = 0; nb < 6; nb++)
        {
            if (FT_StateStructAtGridCrossing(front,grid_intfc,icoords,
				dir[nb],comp,&intfc_state,&hs,crx_coords) &&
                wave_type(hs) != FIRST_PHYSICS_WAVE_TYPE)
		    {
                if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                                boundary_state_function(hs) &&
                                strcmp(boundary_state_function_name(hs),
                                "flowThroughBoundaryState") == 0)
                {
                                U_nb[nb] = vel[l][index];
                }
                else
                {
                    U_nb[nb] = getStateVel[l](intfc_state);
                }

                if (wave_type(hs) == DIRICHLET_BOUNDARY || neumann_type_bdry(wave_type(hs)))
                    mu[nb] = mu0;
                else
                    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    
            }
            else
		    {
                U_nb[nb] = vel[l][index_nb[nb]];
			    mu[nb] = 1.0/2*(mu0 + field->mu[index_nb[nb]]);
		    }
        }

            	coeff[0] = m_dt/rho*mu[0]/(top_h[0]*top_h[0]);
            	coeff[1] = m_dt/rho*mu[1]/(top_h[0]*top_h[0]);
            	coeff[2] = m_dt/rho*mu[2]/(top_h[1]*top_h[1]);
            	coeff[3] = m_dt/rho*mu[3]/(top_h[1]*top_h[1]);
            	coeff[4] = m_dt/rho*mu[4]/(top_h[2]*top_h[2]);
            	coeff[5] = m_dt/rho*mu[5]/(top_h[2]*top_h[2]);

            	getRectangleCenter(index, coords);
            	computeSourceTerm(coords, source);

		aII = 1+coeff[0]+coeff[1]+coeff[2]+coeff[3]+coeff[4]+coeff[5];
		rhs = vel[l][index];

		for(nb = 0; nb < 6; nb++)
		{
	            if (!(*findStateAtCrossing)(front,icoords,dir[nb],comp,
			        &intfc_state,&hs,crx_coords))
			solver.Set_A(I,I_nb[nb],-coeff[nb]);
		    else
		    {
			if (wave_type(hs) == DIRICHLET_BOUNDARY &&
                    	    boundary_state_function(hs) &&
                    	    strcmp(boundary_state_function_name(hs),
                    	    "flowThroughBoundaryState") == 0)
			    aII -= coeff[nb];
			else
			    rhs += coeff[nb]*U_nb[nb];
		    }
		}
		rhs += m_dt*source[l];
		rhs += m_dt*f_surf[l][index];
		//rhs -= m_dt*grad_q[l][index]/rho;
            	solver.Set_A(I,I,aII);

		solver.Set_b(I, rhs);
            }

            solver.SetMaxIter(40000);
            solver.SetTol(1e-10);
                //solver.SetTol(1e-14);

	    start_clock("Befor Petsc solve");
            solver.Solve();
            solver.GetNumIterations(&num_iter);
            solver.GetFinalRelativeResidualNorm(&rel_residual);

	    stop_clock("After Petsc solve");

            // get back the solution
            solver.Get_x(x);

            if (debugging("PETSc"))
                (void) printf("L_CARTESIAN::"
			"computeDiffusionImplicit: "
                        "num_iter = %d, rel_residual = %g. \n",
                        num_iter,rel_residual);

	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                I = ijk_to_I[i][j][k];
                index = d_index3d(i,j,k,top_gmax);
                if (I >= 0)
                    vel[l][index] = x[I-ilower];
                else
                    vel[l][index] = 0.0;
            }
	}
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
        FT_FreeThese(1,x);

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionImplicit()\n");
}       /* end computeDiffusionImplicit */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmI(void)
{
    int index;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	int icrds_max[MAXD],icrds_min[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmI()\n");

    for (int k = 0; k <= top_gmax[2]; k++)
	for (int j = 0; j <= top_gmax[1]; j++)
    for (int i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
        pres[index] += phi[index];
	    q[index] = pres[index];
	    
        if (i < imin || i > imax ||
            j < jmin || j > jmax ||
            k < kmin || k > kmax) continue;

	    if (min_pressure > pres[index])
	    {
            min_pressure = pres[index];
            icrds_min[0] = i;
            icrds_min[1] = j;
            icrds_min[2] = k;
	    }
	    if (max_pressure < pres[index])
	    {
            max_pressure = pres[index];
            icrds_max[0] = i;
            icrds_max[1] = j;
            icrds_max[2] = k;
	    }
	}
    
    //TODO: need to scatter pres and q?

	if (debugging("step_size"))
	{
	    (void) printf(" Max pressure = %f  occuring at: %d %d %d\n",
			max_pressure,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf(" Min pressure = %f  occuring at: %d %d %d\n",
			min_pressure,icrds_min[0],icrds_min[1],icrds_min[2]);
	    (void) printf("Diff pressure = %f\n",max_pressure-min_pressure);
	}

    if (debugging("trace"))
	    (void) printf("Leaving computePressurePmI()\n");
}        /* end computePressurePmI3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmII(void)
{
    int i,j,k,index;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;
    double mu0;
	
    if (debugging("trace"))
	    (void) printf("Entering computePressurePmII()\n");
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
        mu0 = 0.5*field->mu[index];
        pres[index] = q[index] + phi[index] - accum_dt*mu0*div_U[index];
        q[index] = pres[index];

	    if (min_pressure > pres[index])
		min_pressure = pres[index];
	    if (max_pressure < pres[index])
		max_pressure = pres[index];
	}

    //TODO: need to scatter pres and q?

	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmII()\n");
}        /* end computePressurePmII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressurePmIII(void)
{
    int i,j,k,index;
    double mu0;
	double *pres = field->pres;
	double *phi = field->phi;
	double *q = field->q;
	double *div_U = field->div_U;

	if (debugging("trace"))
	    (void) printf("Entering computePressurePmIII()\n");
    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
        mu0 = 0.5*field->mu[index];
        pres[index] = phi[index] - accum_dt*mu0*div_U[index];
	    q[index] = 0.0;

	    if (min_pressure > pres[index])
		min_pressure = pres[index];
	    if (max_pressure < pres[index])
		max_pressure = pres[index];
	}

    //TODO: need to scatter pres?

	if (debugging("trace"))
	    (void) printf("Leaving computePressurePmIII()\n");
}        /* end computePressurePmIII3d */

void Incompress_Solver_Smooth_3D_Cartesian::computePressure(void)
{
    int i,j,k,index;
	double *pres = field->pres;
	min_pressure =  HUGE;
	max_pressure = -HUGE;

    for (k = 0; k <= top_gmax[2]; k++)
	for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
	{
        index = d_index3d(i,j,k,top_gmax);
	    if (min_pressure > pres[index])
            min_pressure = pres[index];
	    if (max_pressure < pres[index])
            max_pressure = pres[index];
	}
	min_pressure =  HUGE;
	max_pressure = -HUGE;
	switch (iFparams->num_scheme.projc_method)
	{
	case BELL_COLELLA:
	    computePressurePmI();
	    break;
	case KIM_MOIN:
	    computePressurePmII();
	    break;
	case SIMPLE:
	case PEROT_BOTELLA:
	    computePressurePmIII();
	    break;
	case ERROR_PROJC_SCHEME:
	default:
	    (void) printf("Unknown computePressure() scheme!\n");
	    clean_up(ERROR);
	}

    computeGradientQ();
}	/* end computePressure */

/*
void Incompress_Solver_Smooth_3D_Cartesian::computeGradientPhi()
{
	int i,j,k,l,index;
	double **grad_phi = field->grad_phi;
	int icoords[MAXD];
	double *phi = field->phi;
	double point_grad_phi[MAXD];

	for (k = 0; k < top_gmax[2]; ++k)
	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = phi[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
        computeFieldPointGrad(icoords,array,point_grad_phi);
	    for (l = 0; l < dim; ++l)
            grad_phi[l][index] = point_grad_phi[l];
	}
	FT_ParallelExchGridVectorArrayBuffer(grad_phi,front);
}*/	/* end computeGradientPhi */

void Incompress_Solver_Smooth_3D_Cartesian::computeGradientQ()
{
	int i,j,k,l,index;
	double **grad_q = field->grad_q;
	int icoords[MAXD];
	double *q = field->q;
	double point_grad_q[MAXD];

	for (k = 0; k < top_gmax[2]; ++k)
	for (j = 0; j < top_gmax[1]; ++j)
	for (i = 0; i < top_gmax[0]; ++i)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    array[index] = q[index];
	}
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
	for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
        computeFieldPointGrad(icoords,array,point_grad_q);
	    for (l = 0; l < dim; ++l)
            grad_q[l][index] = point_grad_q[l];
	}
	FT_ParallelExchGridVectorArrayBuffer(grad_q,front);
}	/* end computeGradientQ */

#define		MAX_TRI_FOR_INTEGRAL		100
void Incompress_Solver_Smooth_3D_Cartesian::surfaceTension(
	double *coords,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *force,
	double sigma)
{
	int i,j,k,num_tris;
	TRI *tri,*tri_list[MAX_TRI_FOR_INTEGRAL];
	double kappa_tmp,kappa,mag_nor,area,delta;
	double median[MAXD],nor[MAXD];
	POINT *p;

	TriAndFirstRing(hse,hs,&num_tris,tri_list);
	for (i = 0; i < num_tris; ++i)
	{
	    kappa = 0.0;
	    tri = tri_list[i];
	    for (j = 0; j < 3; ++j) median[j] = 0.0;
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tri)[j];
		for (k = 0; k < 3; ++k) 
		    median[k] += Coords(p)[k];
	    	GetFrontCurvature(p,Hyper_surf_element(tri),hs,
				&kappa_tmp,front);
		kappa += kappa_tmp;
		nor[j] = Tri_normal(tri)[j];
	    }
	    kappa /= 3.0;
	    mag_nor = mag_vector(nor,3);
	    area = 0.5*mag_nor;
	    for (j = 0; j < 3; ++j)  
	    {
		nor[j] /= mag_nor;
		median[j] /= 3.0;
	    }
	    delta = smoothedDeltaFunction(coords,median);
	    if (delta == 0.0) continue;
	    for (j = 0; j < dim; ++j) 
	    {
		force[j] += delta*sigma*area*kappa*nor[j];
	    }
	}
}	/* end surfaceTension3d */

void Incompress_Solver_Smooth_3D_Cartesian::setInitialCondition()
{
	int i;
	COMPONENT comp;
	double coords[MAXD];
	int size = (int)cell_center.size();
	double *pres = field->pres;
	double *phi = field->phi;

	FT_MakeGridIntfc(front);
	setDomain();

    m_rho[0] = iFparams->rho1;
    m_rho[1] = iFparams->rho2;
    m_mu[0] = iFparams->mu1;
    m_mu[1] = iFparams->mu2;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	m_sigma = iFparams->surf_tension;
	mu_min = rho_min = HUGE;
	for (i = 0; i < 2; ++i)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
        	mu_min = std::min(mu_min,m_mu[i]);
        	rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	// Initialize state at cell_center
        for (i = 0; i < size; i++)
        {
            getRectangleCenter(i, coords);
	    //cell_center[i].m_state.setZero();
	    comp = top_comp[i];
	    if (getInitialState != NULL)
	    {
            //TODO: see comments in getPressure() function
	    	(*getInitialState)(comp,coords,field,i,dim,iFparams);
		    pres[i] = getPressure(front,coords,NULL);
            phi[i] = getPhiFromPres(front,pres[i]);
	    }
        }

    computeGradientQ();
    copyMeshStates();
	setAdvectionDt();
}       /* end setInitialCondition */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjection(void)
{
        switch (iFparams->num_scheme.ellip_method)
        {
        case SIMPLE_ELLIP:
            computeProjectionSimple();
            return;
        case DOUBLE_ELLIP:
            computeProjectionDouble();
            return;
        default:
            printf("Elliptic Method Not Implemented\n");
            clean_up(1);
        }
}       /* end computeProjection */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionDouble(void)
{
}	/* end computeProjectionDouble */

void Incompress_Solver_Smooth_3D_Cartesian::computeProjectionSimple(void)
{
	static ELLIPTIC_SOLVER elliptic_solver(*front);
    int index;
    int i,j,k,l,icoords[MAXD];
    double **vel = field->vel;
    double *phi = field->phi;
    double *div_U = field->div_U;
    double sum_div,L1_div;
    double value;
    double min_phi,max_phi;
	int icrds_max[MAXD],icrds_min[MAXD];

	if (debugging("trace"))
	    (void) printf("Entering computeProjectionSimple()\n");

	for (l = 0; l < dim; ++l)
    {
        vmin[l] = HUGE;
        vmax[l] = -HUGE;
    }

    /* Compute velocity divergence */
    for (k = kmin; k <= kmax; k++)
    for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
    {
        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;
        index  = d_index(icoords,top_gmax,dim);
        source[index] = computeFieldPointDiv(icoords,vel);
        diff_coeff[index] = 1.0/field->rho[index];
        div_U[index] = source[index];
        source[index] /= accum_dt;

        /*Compute pressure jump due to porosity*/
        source[index] += computeFieldPointPressureJump(icoords,
                         iFparams->porous_coeff[0],
                         iFparams->porous_coeff[1]);
        
        array[index] = phi[index];

        if (debugging("check_div"))
        {
            for (l = 0; l < dim; ++l)
            {
                if (vmin[l] > field->vel[l][index])
                    vmin[l] = field->vel[l][index];
                if (vmax[l] < field->vel[l][index])
                    vmax[l] = field->vel[l][index];
            }
        }
    }

    FT_ParallelExchGridArrayBuffer(source,front,NULL);
    FT_ParallelExchGridArrayBuffer(diff_coeff,front,NULL);

    /*
    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index  = d_index3d(i,j,k,top_gmax);
        div_U[index] = source[index];
        source[index] /= accum_dt;
        icoords[0] = i; 
    icoords[1] = j;
    icoords[2] = k;
        source[index] += computeFieldPointPressureJump(icoords,
                         iFparams->porous_coeff[0],
                         iFparams->porous_coeff[1]);
        array[index] = phi[index];
    }
    */

    if(debugging("step_size"))
    {
        sum_div = 0.0;
        min_value =  HUGE;
        max_value = -HUGE;
        min_phi =  HUGE;
        max_phi = -HUGE;

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            value = fabs(div_U[index]);
            sum_div = sum_div + fabs(div_U[index]);
        
            if (min_value > div_U[index]) 
            {
                min_value = div_U[index];
                icrds_min[0] = i;
                icrds_min[1] = j;
                icrds_min[2] = k;
            }

            if (max_value < div_U[index]) 
            {
                max_value = div_U[index];
                icrds_max[0] = i;
                icrds_max[1] = j;
                icrds_max[2] = k;
            }
        }

        pp_global_min(&min_value,1);
        pp_global_max(&max_value,1);
        L1_div = sum_div/(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1);

        (void) printf("Before computeProjection:\n");
        (void) printf("Sum div(U) = %f\n",sum_div);
        (void) printf("Min div(U) = %f  ",min_value);
        (void) printf("occuring at: %d %d %d\n",icrds_min[0],icrds_min[1],
                icrds_min[2]);
        (void) printf("Max div(U) = %f  ",max_value);
        (void) printf("occuring at: %d %d %d\n",icrds_max[0],icrds_max[1],
                icrds_max[2]);
        (void) printf("L1  div(U) = %f\n",L1_div);
    }

    if (debugging("check_div"))
    {
        checkVelocityDiv("Before computeProjection()");
    }
    
    elliptic_solver.D = diff_coeff;
    elliptic_solver.source = source;
    elliptic_solver.soln = array;
    elliptic_solver.set_solver_domain();
    //TODO: Need getStatePres() for interpolation when
    //      applying boundary conditions at solid walls?
    elliptic_solver.getStateVar = getStatePhi;
    elliptic_solver.findStateAtCrossing = findStateAtCrossing;
	elliptic_solver.skip_neumann_solver = skip_neumann_solver;
	
    paintAllGridPoint(TO_SOLVE);
    setGlobalIndex();
    setIndexMap();

    elliptic_solver.ijk_to_I = ijk_to_I;
    elliptic_solver.ilower = ilower;
    elliptic_solver.iupper = iupper;
    elliptic_solver.skip_neumann_solver = skip_neumann_solver;
    elliptic_solver.solve(array);

	FT_ParallelExchGridArrayBuffer(array,front,NULL);

    /*
	if (iFparams->with_porosity)
	{
	    paintAllGridPoint(TO_SOLVE);
	    setGlobalIndex();
	    setIndexMap();
        elliptic_solver.ijk_to_I = ijk_to_I;
        elliptic_solver.ilower = ilower;
        elliptic_solver.iupper = iupper;
        elliptic_solver.skip_neumann_solver = skip_neumann_solver;
        elliptic_solver.solve(array);
	}
    else
	{	
	    int num_colors = drawColorMap();
        std::vector<int> ncell;
        ncell.resize(num_colors);
        paintAllGridPoint(NOT_SOLVED);
        for (i = 1; i < num_colors; ++i)
        {
            paintToSolveGridPoint2(i);
            setGlobalIndex();
            setIndexMap();
            elliptic_solver.ijk_to_I = ijk_to_I;
            elliptic_solver.ilower = ilower;
            elliptic_solver.iupper = iupper;

            ncell[i] = iupper - ilower;
            printf("ilower = %d  iupper = %d\n",ilower,iupper);
            
            elliptic_solver.solve(array);
            paintSolvedGridPoint();
        }
        //for (i = 1; i < num_colors; ++i)
        //{
        //    if (ncell[i] > 10) continue;
        //    setIsolatedSoln(i,array);
        //}

	}
	FT_ParallelExchGridArrayBuffer(array,front,NULL);
    */

	min_phi =  HUGE;
	max_phi = -HUGE;
    for (k = 0; k <= top_gmax[2]; k++)
    for (j = 0; j <= top_gmax[1]; j++)
    for (i = 0; i <= top_gmax[0]; i++)
    {
        index  = d_index3d(i,j,k,top_gmax);
        phi[index] = array[index];

        if (min_phi > phi[index])
        {
            min_phi = phi[index];
            icrds_min[0] = i;
            icrds_min[1] = j;
            icrds_min[2] = k;
        }
        
        if (max_phi < phi[index])
        {
            max_phi = phi[index];
            icrds_max[0] = i;
            icrds_max[1] = j;
            icrds_max[2] = k;
        }
    }
        
    if(debugging("projection"))
	{
	    (void) printf("After computeProjection:\n");
	    (void) printf("min_phi = %f  occuring at: %d %d %d\n",
			min_phi,icrds_min[0],icrds_min[1],icrds_min[2]);
	    (void) printf("max_phi = %f  occuring at: %d %d %d\n",
			max_phi,icrds_max[0],icrds_max[1],icrds_max[2]);
	    (void) printf("diff_phi = %f\n",max_phi-min_phi);
	}
	if (debugging("trace"))
	    (void) printf("Leaving computeProjectionSimple()\n");
}	/* end computeProjectionSimple */

void Incompress_Solver_Smooth_3D_Cartesian::solveTest(const char *msg)
{
	// This function is reserved for various debugging tests.
	int k,index;
	int icoords[MAXD];
	double **vel = field->vel;
	(void) printf("%s\n",msg);

    for (k = kmin; k <= kmax; k++)
    {
        icoords[0] = 30;
        icoords[1] = 30;
        icoords[2] = k;
        index  = d_index(icoords,top_gmax,dim);
        source[index] = computeFieldPointDiv(icoords,vel);
        source[index] /= m_dt;
        printf("div[%d]/dt = %14.8f\n",k,source[index]);
    }
}	/* end solveTest */

static int parab_find_state_at_crossing(Front*,int*,GRID_DIRECTION,int,
                                POINTER*,HYPER_SURF**,double*);

void Incompress_Solver_Smooth_3D_Cartesian::
	computeDiffusionParab(void)
{
	static PARABOLIC_SOLVER parab_solver(*front);
	static boolean first = YES;

        COMPONENT comp;
        int index;
	int i,j,k,l,icoords[MAXD];
	double source[MAXD];
	double **vel = field->vel;
	double **f_surf = field->f_surf;
	static double *nu;

	if (debugging("trace"))
	    (void) printf("Entering Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionParab()\n");

	if (nu == NULL)
	{
	    int size = 1;
	    for (i = 0; i < dim; ++i)  size *= (top_gmax[i] + 1);
	    FT_VectorMemoryAlloc((POINTER*)&nu,size,sizeof(double));
	}
        setIndexMap();
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index  = d_index3d(i,j,k,top_gmax);
	    nu[index] = field->mu[index]/field->rho[index];
	}
	FT_ParallelExchGridArrayBuffer(nu,front,NULL);

	parab_solver.soln_comp = LIQUID_COMP2;
	parab_solver.obst_comp = SOLID_COMP;
	parab_solver.ilower = ilower;
        parab_solver.iupper = iupper;
	parab_solver.dt = m_dt;
	parab_solver.order = 2;
	parab_solver.a = NULL;
	parab_solver.findStateAtCrossing = parab_find_state_at_crossing;
	parab_solver.first = first;
	parab_solver.set_solver_domain();
	first = NO;
	switch(dim)
        {
        case 2:
            parab_solver.ij_to_I = ij_to_I;
            break;
        case 3:
            parab_solver.ijk_to_I = ijk_to_I;
            break;
        }
	
	for (l = 0; l < dim; ++l)
	{
	    parab_solver.var = vel[l];
	    parab_solver.soln = vel[l];
	    parab_solver.getStateVarFunc = getStateVel[l];
	    parab_solver.solveIM();
	}

	if (debugging("trace"))
	    (void) printf("Leaving Incompress_Solver_Smooth_3D_Cartesian::"
			"computeDiffusionParab()\n");
}       /* end computeDiffusionParab */

//TODO: get this working
void Incompress_Solver_Smooth_3D_Cartesian::vtk_plot_scalar(
        char *outname, const char* varname)
{
        std::vector<int> ph_index;
        int i,j,k,index;
        char dirname[256],filename[256];
        FILE *outfile;
        double coord_x,coord_y,coord_z,xmin,ymin,zmin;
        COMPONENT comp;
        int pointsx,pointsy,pointsz,num_points,num_cells,num_cell_list;
        int icoords[3],p_gmax[3];

        int ii,jj,kk;
        double ih,jh,kh;

        sprintf(filename, "%s/vtk/vtk.ts%s",outname,
                right_flush(front->step,7));
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
        //cell-based liquid phase
        ph_index.clear();
        if (!create_directory(filename,NO))
        {
            printf("Cannot create directory %s\n",filename);
            clean_up(ERROR);
        }
        sprintf(filename,"%s/%s.vtk",filename,varname);
        outfile = fopen(filename,"w");
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"%s\n",varname);
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");

        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            if (ifluid_comp(cell_center[index].comp))
                ph_index.push_back(index);
        }

        pointsx = top_gmax[0] + 2;
        pointsy = top_gmax[1] + 2;
        pointsz = top_gmax[2] + 2;
        num_points = pointsx*pointsy*pointsz;

        num_cells = (int)ph_index.size();
        num_cell_list = 9*num_cells;
        fprintf(outfile,"POINTS %d double\n", num_points);

        for (k = 0; k <= top_gmax[2]; k++)
        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }
        for (i = 0; i <= top_gmax[0]; i++)
        for (j = 0; j <= top_gmax[1]; j++)
        {
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (i = 0; i <= top_gmax[0]; i++)
        {
            j = top_gmax[1];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] - top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        for (j = 0; j <= top_gmax[1]; j++)
        {
            i = top_gmax[0];
            k = top_gmax[2];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] - top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }
        for (k = 0; k <= top_gmax[2]; k++)
        {
            i = top_gmax[0];
            j = top_gmax[1];
            index = d_index3d(i,j,k,top_gmax);
            coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
            coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
            coord_z = cell_center[index].m_coords[2] - top_h[2]/2.0;
            fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);
        }

        i = top_gmax[0];
        j = top_gmax[1];
        k = top_gmax[2];
        index = d_index3d(i,j,k,top_gmax);
        coord_x = cell_center[index].m_coords[0] + top_h[0]/2.0;
        coord_y = cell_center[index].m_coords[1] + top_h[1]/2.0;
        coord_z = cell_center[index].m_coords[2] + top_h[2]/2.0;
        fprintf(outfile,"%f %f %f\n",coord_x,coord_y,coord_z);


        fprintf(outfile,"CELLS %i %i\n", num_cells,num_cell_list);
        for (i = 0; i < num_cells; i++)
        {
            int index0,index1,index2,index3,index4,index5,index6,index7;
            index = ph_index[i];
            icoords[0] = cell_center[index].icoords[0];
            icoords[1] = cell_center[index].icoords[1];
            icoords[2] = cell_center[index].icoords[2];
            index0 = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
            index1 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2],top_gmax);
            index2 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2],top_gmax);
            index3 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2],top_gmax);
            index4 =
                d_index3d(icoords[0],icoords[1],icoords[2]+1,top_gmax);
            index5 =
                d_index3d(icoords[0]+1,icoords[1],icoords[2]+1,top_gmax);
            index6 =
                d_index3d(icoords[0],icoords[1]+1,icoords[2]+1,top_gmax);
            index7 =
                d_index3d(icoords[0]+1,icoords[1]+1,icoords[2]+1,top_gmax);

            fprintf(outfile,"8 %i %i %i %i %i %i %i %i\n",
                index0,index1,index2,index3,index4,index5,index6,index7);
        }

        fprintf(outfile, "CELL_TYPES %i\n", num_cells);
        for (i = 0; i < num_cells; i++)
            fprintf(outfile,"11\n");

        fprintf(outfile, "CELL_DATA %i\n", num_cells);
        fprintf(outfile, "SCALARS %s double\n",varname);
        fprintf(outfile, "LOOKUP_TABLE default\n");
        for (i = 0; i < num_cells; i++)
        {
            index = ph_index[i];
            if(strcmp(varname,"pres") == 0)
              fprintf(outfile,"%f\n",field->pres[index]);
        }
        fclose(outfile);
}

static int parab_find_state_at_crossing(
	Front *front,
        int *icoords,
        GRID_DIRECTION dir,
        int comp,
        POINTER *state,
        HYPER_SURF **hs,
        double *crx_coords)
{
	boolean status;
        INTERFACE *grid_intfc = front->grid_intfc;

        status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
                                comp,state,hs,crx_coords);
        if (status == NO) return NO_PDE_BOUNDARY;

	switch (wave_type(*hs))
	{
	case FIRST_PHYSICS_WAVE_TYPE:
	    return NO_PDE_BOUNDARY;
	case DIRICHLET_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	case NEUMANN_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
            return DIRICHLET_PDE_BOUNDARY;
	default:
	    (void) printf("In parab_find_state_at_crossing()\n");
	    (void) printf("Unknown wave type %s\n",
				f_wave_type_as_string(wave_type(*hs)));
	    clean_up(ERROR);
        }
}	/* end parab_find_state_at_crossing */

void Incompress_Solver_Smooth_3D_Cartesian::setParallelVelocity(void)
{
        FILE *infile;
        int i,j,id,k,l,index,G_index;
        char fname[100];
        COMPONENT comp;
        double coords[MAXD];
        int size = (int)cell_center.size();
        int myid = pp_mynode();
        int numprocs = pp_numnodes();

        int G_icoords[MAXD],pp_icoords[MAXD],icoords[MAXD];
        int local_gmax[MAXD], global_gmax[MAXD];
        int G_size, L_size;
        PP_GRID *pp_grid = front->pp_grid;
        double *local_L = pp_grid->Zoom_grid.L;
        double *local_U = pp_grid->Zoom_grid.U;
        double *GU_buff,*GV_buff, *GW_buff, *U_buff, *V_buff, *W_buff;

        for (i = 0; i < dim; i++)
        {
            global_gmax[i] = pp_grid->Global_grid.gmax[i]-1;
            local_gmax[i] = pp_grid->Zoom_grid.gmax[i]-1;
        }
        FT_MakeGridIntfc(front);
        setDomain();
        G_size = 1;
        L_size = 1;
        for (i = 0; i < dim; i++)
        {
            G_size = G_size * (global_gmax[i]+1);
            L_size = L_size * (top_gmax[i]+1);
        }
        uni_array(&U_buff,L_size,sizeof(double));
        uni_array(&V_buff,L_size,sizeof(double));
        uni_array(&W_buff,L_size,sizeof(double));
        if (myid == 0)
        {
            uni_array(&GU_buff,G_size,sizeof(double));
            uni_array(&GV_buff,G_size,sizeof(double));
            uni_array(&GW_buff,G_size,sizeof(double));

	    if (setInitialVelocity != NULL)
                (*setInitialVelocity)(comp,pp_grid->Global_grid.gmax,
				   GU_buff,GV_buff,GW_buff,
				   &(pp_grid->Global_grid),iFparams);
            for (id = 0; id < numprocs; id++)
            {
                find_Cartesian_coordinates(id,pp_grid,pp_icoords);
                for (k = kmin; k <= kmax; ++k)
                for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                {
                    icoords[0] = i;
                    icoords[1] = j;
		    icoords[2] = k;
                    G_icoords[0] = pp_icoords[0]*(local_gmax[0]+1)+icoords[0]-imin;
                    G_icoords[1] = pp_icoords[1]*(local_gmax[1]+1)+icoords[1]-jmin;
                    G_icoords[2] = pp_icoords[2]*(local_gmax[2]+1)+icoords[2]-kmin;
                    G_index = d_index(G_icoords,global_gmax,dim);
                    index = d_index(icoords,top_gmax,dim);
                    U_buff[index] = GU_buff[G_index];
                    V_buff[index] = GV_buff[G_index];
                    W_buff[index] = GW_buff[G_index];
                }
                if (id == 0)
                {
                    for (i = 0; i < L_size; i++)
                    {
                        field->vel[0][i] = U_buff[i];
                        field->vel[1][i] = V_buff[i];
                        field->vel[2][i] = W_buff[i];
                    }
                }
                else
                {
                    pp_send(1,(POINTER)(U_buff),sizeof(double)*L_size,id);
                    pp_send(2,(POINTER)(V_buff),sizeof(double)*L_size,id);
                    pp_send(3,(POINTER)(W_buff),sizeof(double)*L_size,id);
                }
            }
            FT_FreeThese(3,GU_buff,GV_buff,GW_buff);
        }
        else
        {
            pp_recv(1,0,(POINTER)(U_buff),sizeof(double)*L_size);
            pp_recv(2,0,(POINTER)(V_buff),sizeof(double)*L_size);
            pp_recv(3,0,(POINTER)(W_buff),sizeof(double)*L_size);
            for (i = 0; i < L_size; i++)
            {
                field->vel[0][i] = U_buff[i];
                field->vel[1][i] = V_buff[i];
                field->vel[2][i] = W_buff[i];
            }
        }


        m_rho[0] = iFparams->rho1;
        m_rho[1] = iFparams->rho2;
        m_mu[0] = iFparams->mu1;
        m_mu[1] = iFparams->mu2;
        m_comp[0] = iFparams->m_comp1;
        m_comp[1] = iFparams->m_comp2;
        m_smoothing_radius = iFparams->smoothing_radius;
        m_sigma = iFparams->surf_tension;
        mu_min = rho_min = HUGE;
        for (i = 0; i < 2; ++i)
        {
            if (ifluid_comp(m_comp[i]))
            {
                mu_min = std::min(mu_min,m_mu[i]);
                rho_min = std::min(rho_min,m_rho[i]);
            }
        }
        FT_FreeThese(3,U_buff,V_buff,W_buff);

        computeGradientQ();
        copyMeshStates();
        setAdvectionDt();
}

void Incompress_Solver_Smooth_3D_Cartesian::extractFlowThroughVelocity()
{
	int index,index_nb,index_op;
	int i,j,k,l,nb,icoords[MAXD],icoords_nb[MAXD],icoords_op[MAXD];
	GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	POINTER intfc_state;
        HYPER_SURF *hs;
	double crx_coords[MAXD];
	double vel_save,**vel = field->vel;
	double div;
	int status;

	/* Extract velocity for zero interior divergence */
	for (k = 0; k <= top_gmax[2]; ++k)
        for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
        {
            if (i >= imin && i <= imax && j >= jmin &&
                j <= jmax && k >= kmin && k <= kmax)
                continue;
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index])) continue;
            for (l = 0; l < dim; ++l)
                icoords_nb[l] = icoords_op[l] = icoords[l];
            for (l = 0; l < dim; ++l)
            for (nb = 0; nb < 2; ++nb)
            {
                status = (*findStateAtCrossing)(front,icoords,dir[l][nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY) continue;
                icoords_nb[l] = (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;
                icoords_op[l] = (nb == 0) ? icoords[l] - 2 : icoords[l] + 2;
                index_nb = d_index(icoords_nb,top_gmax,dim);
                index_op = d_index(icoords_op,top_gmax,dim);
                if (!ifluid_comp(top_comp[index_nb]))
                    continue;
                vel[l][index] = 0.0;
                div = computeFieldPointDiv(icoords_nb,vel);
                vel[l][index] = (nb == 0) ? vel[l][index] - 2.0*top_h[l]*div
                                : vel[l][index] + 2.0*top_h[l]*div;
            }
        }

}	/* end extractFlowThroughVelocity */

void Incompress_Solver_Smooth_3D_Cartesian::computeVelDivergence()
{
	double *div_U = field->div_U;
	double **vel = field->vel;
	int i,j,k,index,icoords[MAXD];
	double Lnorm[3];

	/* Compute velocity divergence */
	for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    icoords[0] = i;
	    icoords[1] = j;
	    icoords[2] = k;
	    index  = d_index(icoords,top_gmax,dim);
	    if (!ifluid_comp(top_comp[index]))
		div_U[index] = 0.0;
	    div_U[index] = computeFieldPointDiv(icoords,vel);
	}
}	/* end computeVelDivergence */

void Incompress_Solver_Smooth_3D_Cartesian::computeVarIncrement(
	double *var_old,
	double *var_new,
	boolean use_dual_grid)
{
	int i,j,k,index,size;
	double mag;
	double Lnorm[3];

	if (use_dual_grid)
	{
	    ;	// To add dual grid computation.
	}
	else
	{
	    Lnorm[0] = Lnorm[1] = Lnorm[2] = 0.0;
	    size = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                index = d_index3d(i,j,k,top_gmax);
		mag = fabs(var_new[index] - var_old[index]);
		Lnorm[0] += mag;
		Lnorm[1] += sqr(mag);
		if (Lnorm[2] < mag)
		{
		    Lnorm[2] = mag;
		}
            }
	    Lnorm[0] /= size;
	    Lnorm[1] = sqrt(Lnorm[1]/size);
	    (void) printf("L-1 norm = %20.14f  L-1/dt = %20.14f\n",
					Lnorm[0],Lnorm[0]/m_dt);
	    (void) printf("L-2 norm = %20.14f  L-2/dt = %20.14f\n",
					Lnorm[1],Lnorm[1]/m_dt);
	    (void) printf("L-I norm = %20.14f  L-I/dt = %20.14f\n",
					Lnorm[2],Lnorm[2]/m_dt);
	    Lnorm[0] = 0.0;
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		mag = fabs(var_new[index]);
		Lnorm[0] += mag;
	    }
	    Lnorm[0] /= size;
	    (void) printf("L-1 norm of old variable = %20.14f\n",Lnorm[0]);
	}
}	/* end computeVarIncrement */	

void Incompress_Solver_Smooth_3D_Cartesian::appendOpenEndStates()
{
        INTERFACE *intfc = front->interf;
        int dim = front->rect_grid->dim;
        int i,j,k,idir,side,ii;
        int ic[MAXD], comp, index, bdry_type;
        STATE state;

        if (debugging("trace"))
            printf("Entering appendOpenEndStates() \n");
        if (dim != 3) return;

        for (idir = 0; idir < dim; ++idir)
        for (side = 0; side < 2; ++side)
        {
            if (rect_boundary_type(intfc,idir,side) == OPEN_BOUNDARY &&
                front->open_end_func != NULL)
            {
                //TODO: Debug statement, to see if this is actually called.
                //      If not, need to make changes so that the boundary
                //      states at the outlet get updated.
                //      This is potentially being taken care of by the
                //      IF_flowThroughBoundary function, and is just dead code.
                printf("appendOpenEndStates(): OPEN_BOUNDARY\n");
                int count = 0;
                for (i = 0; i < top_gmax[0]; ++i)
                for (j = 0; j < top_gmax[1]; ++j)
                for (k = 0; k < top_gmax[2]; ++k)
                {
                    ic[0] = i; ic[1] = j; ic[2] = k;
                    if ((side == 0 && ic[idir] < lbuf[idir]) ||
                        (side == 1 && ic[idir] > top_gmax[idir]-ubuf[idir]))
                    {
                        index = d_index3d(i,j,k,top_gmax);
                        comp = top_comp[index];
                        front->open_end_func(front,front->open_end_params,
                                ic,comp,idir,side,&bdry_type,&state);
                        field->rho[index] = state.dens;
                        for (ii = 0; ii < 3; ++ii)
                        {
                            field->vel[ii][index] = state.vel[ii];
                        }
                        field->pres[index] = state.pres;
                        field->phi[index] = state.phi; 
                    }
                }
            }
        }
        if (debugging("trace"))
            printf("Leaving appendOpenEndStates() \n");
        return;
}       /* end appendOpenEndStates */

void Incompress_Solver_Smooth_Basis::setIsolatedSoln(
        int color, 
        double *soln)
{
	int i,j,k,l,ic,count;
        int ib_min,ib_max,jb_min,jb_max,kb_min,kb_max;
        double ave_soln;
        std::vector<int> iso_cell;
        COMPONENT c;

        ib_min = imax;  ib_max = imin;
        jb_min = jmax;  jb_max = jmin;
        kb_min = kmax;  kb_max = kmin;
        switch(dim)
        {
        case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                ic = d_index2d(i,j,top_gmax);
		if (color_map[ic] == color)
                    iso_cell.push_back(ic);
            }
            break;
        case 3:
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                ic = d_index3d(i,j,k,top_gmax);
		if (color_map[ic] == color)
                {
                    if (ib_min > i) ib_min = i;
                    if (jb_min > j) jb_min = j;
                    if (kb_min > k) kb_min = k;
                    if (ib_max < i) ib_max = i;
                    if (jb_max < j) jb_max = j;
                    if (kb_max < k) kb_max = k;
                    iso_cell.push_back(ic);
                }
            }
            printf("Entering setIsolatedSoln()\n");
            printf("size = %d\n",iso_cell.size());
            printf("ib_min = %d  ib_max = %d\n",ib_min,ib_max);
            printf("jb_min = %d  jb_max = %d\n",jb_min,jb_max);
            printf("kb_min = %d  kb_max = %d\n",kb_min,kb_max);
            ib_min = (ib_min == imin) ? ib_min : ib_min - 1;
            jb_min = (jb_min == jmin) ? jb_min : jb_min - 1;
            kb_min = (kb_min == kmin) ? kb_min : kb_min - 1;
            ib_max = (ib_max == imax) ? ib_max : ib_max + 1;
            jb_max = (jb_max == jmax) ? jb_max : jb_max + 1;
            kb_max = (kb_max == kmax) ? kb_max : kb_max + 1;
            printf("Revised:\n");
            printf("ib_min = %d  ib_max = %d\n",ib_min,ib_max);
            printf("jb_min = %d  jb_max = %d\n",jb_min,jb_max);
            printf("kb_min = %d  kb_max = %d\n",kb_min,kb_max);
            ave_soln = 0.0;
            count = 0;
            for (k = kb_min; k <= kb_max; k++)
            for (j = jb_min; j <= jb_max; j++)
            for (i = ib_min; i <= ib_max; i++)
            {
                ic = d_index3d(i,j,k,top_gmax);
                for (l = 0; l < iso_cell.size(); ++l)
                    if (ic == iso_cell[l]) break;
                if (l < iso_cell.size()) continue; 
                if (!ifluid_comp(top_comp[ic])) continue;
                ave_soln += soln[ic];
                count++;
            }
            ave_soln /= count;
            printf("count = %d  ave_soln = %f\n",count,ave_soln);
            //clean_up(0);
            break;
        }
}       /* end setIsolatedSoln */

void Incompress_Solver_Smooth_3D_Basis::addImmersedForce()
{
	INTERFACE *grid_intfc = front->grid_intfc;
	COMPONENT comp = grid_intfc->default_comp;
        
    CURVE **c,*curve;
    BOND *b;
    POINT *p;
    double local_vel[MAXD];
	double **f_surf = field->f_surf;
	double **vel = field->vel;
    double coords[MAXD];
    int icoords[MAXD],ic;

	double center[MAXD],point[MAXD],H,D;
	double *rho = field->rho;
	double dist,alpha;
	int range = (int)(m_smoothing_radius+1);

    if (debugging("trace"))
        printf("Entering addImmersedForce()\n");

    IF_PARAMS *iFparams = (IF_PARAMS*)front->extra1;
    double rhoF = iFparams->rho2;

    intfc_curve_loop(grid_intfc,c)
    {
        if (hsbdry_type(*c) == STRING_HSBDRY)
        {
            curve = *c;
            FINITE_STRING *params = (FINITE_STRING*)curve->extra;
            if (params == NULL) continue;
            
            double c_drag = params->c_drag;
            double radius = params->radius;
            double rhoS = params->dens;
            double ampFluidFactor = params->ampFluidFactor;

            for (b = curve->first; b != curve->last; b = b->next)
            {
                p = b->end;
                STATE* state_intfc = (STATE*)left_state(p);

                //TODO: top_grid the correct grid?
                rect_in_which(Coords(p),icoords,top_grid);
                //int index = d_index(icoords,top_gmax,3);

                //tangential direction along string BOND
                double ldir[3];
                for (int i = 0; i < 3; ++i)	
                    ldir[i] = Coords(b->end)[i] - Coords(b->start)[i];
                double length = Mag3d(ldir);
                if (length < MACH_EPS)
                {
                    printf("BOND length < MACH_EPS\n");
                    clean_up(EXIT_FAILURE);
                }
                
                for (int i = 0; i < 3; ++i)
                    ldir[i] /= length;

                double* vel_intfc = state_intfc->vel;
                double vt = 0.0;
                double vfluid[3], vrel[3];

                for (int i = 0; i < 3; ++i)
                {
                    FT_IntrpStateVarAtCoords(front,NO_COMP,Coords(p),
                            vel[i],getStateVel[i],&vfluid[i],&state_intfc->vel[i]);
                    vrel[i] = vfluid[i] - vel_intfc[i];
                    vt += vrel[i]*ldir[i];
                }

                double speed = 0.0;
                double vtan[3], vnor[3];
                for (int i = 0; i < 3; ++i)
                {
                    vnor[i] = vrel[i] - vt*ldir[i];
                    speed += sqr(vnor[i]);
                }
                speed = sqrt(speed);

                //double A_ref = 2.0*PI*radius*length;
                //double Vol = PI*radius*radius*length;
                double A_ref = 2.0*PI*radius*(0.25*length);
                double Vol = PI*radius*radius*(0.25*length);
                double mass = rhoS*Vol;

                double VolFluid = top_h[0]*top_h[1]*top_h[2];
                double massFluid = rhoF*VolFluid;
                
                double dragForce[MAXD] = {0.0};
                if (front->step > iFparams->fsi_startstep)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        dragForce[i] = 0.5*rhoF*c_drag*A_ref*speed*vnor[i];
                        dragForce[i] *= ampFluidFactor;
                        //dragForce[i] *= massFluid/mass;
                    }
                }

                /*
                //TODO: Put inside debugging string block
                printf("pt = %f %f %f \n",Coords(p)[0],Coords(p)[1],Coords(p)[2]);
                printf("\tdragForce = %g %g %g \n",dragForce[0],dragForce[1],dragForce[2]);
                printf("\tc_drag = %f  |  A_ref = %g  |  rhoF = %g \n",c_drag,A_ref,rhoF);
                printf("\tspeed = %f\n",speed);

                printf("\t\tstate->linedrag_force = %g %g %g \n",
                        state_intfc->linedrag_force[0],
                        state_intfc->linedrag_force[1],
                        state_intfc->linedrag_force[2]);
                */

                /*
                double dragForce[MAXD] = {0.0};
                double* ldragForce = state_intfc->linedrag_force;
                for (int l = 0; l < 3; ++l)
                {
                    dragForce[l] = ampFluidFactor*ldragForce[l];
                    //ldragForce[l] = 0.0;
                }
                */

                for (int k = icoords[2]-2; k <= icoords[2]+2; k++)
                for (int j = icoords[1]-2; j <= icoords[1]+2; j++)
                for (int i = icoords[0]-2; i <= icoords[0]+2; i++)
                {
                    if (i < 0 || j < 0 || k < 0) continue;
                    if (i > top_gmax[0] ||
                        j > top_gmax[1] ||
                        k > top_gmax[2]) continue;

                    coords[0] = top_L[0] + i*top_h[0];
                    coords[1] = top_L[1] + j*top_h[1];
                    coords[2] = top_L[2] + k*top_h[2];
                    ic = d_index3d(i,j,k,top_gmax);

                    //TODO: use smoothing function, see usage example
                    //      in setSmoothedProperties() in surface tension
                    //      section.
                    //
                    //TODO: radius should be factored into the smoothing operation

                    dist = distance_between_positions(Coords(p),coords,3);
                    if (dist >= top_h[0]*4.0) continue;

                    alpha = top_h[0]*4.0 - dist;
                    alpha /= top_h[0]*4.0;

                    for (int l = 0; l < dim; ++l)
                    {
                        f_surf[l][ic] -= alpha*dragForce[l]/rhoF;
                            //ldragForce[l] = 0.0;//state_intfc->linedrag_force[l] = 0.0;
                    }

                    //TODO: Put inside debugging string block
                    /*
                    printf("crds = %f %f %f d = %f alpha = %f f_surf -= %f %f %f\n",
                            coords[0],coords[1],coords[2],dist,alpha,
                            alpha*dragForce[0]/rhoF,alpha*dragForce[1]/rhoF,alpha*dragForce[2]/rhoF);
                    */
                }
            }
        }
    }

	FT_ParallelExchGridVectorArrayBuffer(f_surf,front);
}	/* end addImmersedForce in 3D */

void Incompress_Solver_Smooth_Basis::addVortexDisturbance(
        const VPARAMS& vparams)
{
	    double **vel = field->vel;
        double center[3] = {vparams.center[0],
                            vparams.center[1],
                            vparams.center[2]};

        double D = vparams.D;
        double A = vparams.A;
        double L0[MAXD],coords[MAXD];
        int i,j,k,ic;

        printf("Entering addVortexDisturbance()\n");
        for (i = 0; i < dim; ++i)
            L0[i] = center[i] - 0.5*D;
        for (k = kmin; k <= kmax; k++)
        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            coords[0] = top_L[0] + i*top_h[0];
            coords[1] = top_L[1] + j*top_h[1];
            coords[2] = top_L[2] + k*top_h[2];
            if (coords[0] < L0[0] || coords[0] > L0[0]+D ||
                coords[2] < L0[2] || coords[2] > L0[2]+D)
                continue;
            double u,w;
            u =  A*sin(PI*(coords[0] - L0[0])/D)*cos(PI*(coords[2] - L0[2])/D);
            w = -A*cos(PI*(coords[0] - L0[0])/D)*sin(PI*(coords[2] - L0[2])/D);
            u *= sin(PI*(coords[1] - L0[1])/D);
            w *= sin(PI*(coords[1] - L0[1])/D);
            ic = d_index3d(i,j,k,top_gmax);
            vel[0][ic] += u;
            vel[2][ic] += w;
        }
	FT_ParallelExchGridVectorArrayBuffer(vel,front);
}       /* end addVortexDisturbance */

