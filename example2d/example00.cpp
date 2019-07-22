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

/*
*				example00.cpp:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a slotted disk in rotation. It demonstrates
*	the geometry preservation of the front tracking method.
*
*/

#include <vector>
#include <FronTier.h>

static void test_propagate(Front*);
static int rotation_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
typedef struct {
        double omega_0;  // angular velocity
        double cen[2];  // rotation center's coordinates
} ROTATION_VEL_PARAMS;

int main(int argc, char **argv)
{
//  STEP 1) Process the command line arguments and
//  store parameters in an instance of F_BASIC_DATA
	static F_BASIC_DATA f_basic;
	FT_Init(argc,argv,&f_basic);

//  2) Store information on the computational grid in f_basic.
	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; // lower bounds on dim x and dim y
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; // upper bounds on dim x and dim y
	f_basic.gmax[0] = 128;	f_basic.gmax[1] = 128; // number of grids for each dimension
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY; // x-dim bottom, x-dim top
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY; // y-dim bottom, y-dim top
	f_basic.size_of_intfc_state = 0; // pass value, required but unused

//  3) Construct the Front object from these initial settings.
	static Front front;
	FT_StartUp(&front,&f_basic); // initialize

//	4) Initialize the disk pre-made interface
	TDISK_PARAMS disk_params;
	disk_params.x0 = 0.5; // center of disk
	disk_params.y0 = 0.5; // y coordinate of disk center
	disk_params.r = 0.3;  // disk radius
	disk_params.w = 0.1; // slot width
	disk_params.h = 0.2;  // slot height

//   Note: an interface is described by a level function and its
//   parameters.  disk_params is the paramter object for the
//	 slotted_disk_func level function (below).

//   double slotted_disk_func(
//	     POINTER func_params,
//       double *coords)
//   {

//   returns a positive number on the outside of the interface,
//   and a negative number on the inside.  Function output
//   decreases to 0.0 at the interface boundary.

//   For a 3D rendering see:
//   http://www.ams.stonybrook.edu/~jpetrill/experiments/10/index.html

//  5) attach the level function and parameters to a LEVEL_FUNC_PACK
	static LEVEL_FUNC_PACK level_func_pack;
	level_func_pack.func_params = (POINTER)&disk_params; // attach parameters
	level_func_pack.func = slotted_disk_func; // attach level function

//  label inside and outside.  Must be nonnegative, nonequal integers.
	level_func_pack.neg_component = 1; // inside of the interface
	level_func_pack.pos_component = 2; // outside of the interface

	level_func_pack.wave_type = FIRST_PHYSICS_WAVE_TYPE;  // required
	FT_InitIntfc(&front,&level_func_pack);  // initialize

//  6) Set parameters for velocity field
	ROTATION_VEL_PARAMS rv_params;
	rv_params.cen[0] = 0.5;
	rv_params.cen[1] = 0.5;
	rv_params.omega_0 = -2.0*PI;

//  7) Combine velocity field parameters with velocity function in a VELO_FUNC_PACK
	static VELO_FUNC_PACK velo_func_pack;
	velo_func_pack.func_params = (POINTER)&rv_params; // cast as POINTER (void*)
	velo_func_pack.func = rotation_vel;
	velo_func_pack.point_propagate = fourth_order_point_propagate;
	FT_InitFrontVeloFunc(&front,&velo_func_pack); // initialize

//  8) Propagate the front and save images
	test_propagate(&front);

	clean_up(0);
	return 0;
}


static  void test_propagate(
        Front *front)
{
    double CFL;

	front->max_time = 3;
	front->max_step = 10000;
	front->print_time_interval = 2.0;
	front->movie_frame_interval = 0.1; // 31 frames over time

    CFL = Time_step_factor(front);
	Frequency_of_redistribution(front,GENERAL_WAVE) = 1000;

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

    FT_RedistMesh(front);
	FT_ResetTime(front);

	// Always output the initial interface.
	FT_Save(front);
    FT_Draw(front);

	// This is a virtual propagation to get maximum front
	// speed to determine the first time step.

    FT_Propagate(front);
    FT_SetTimeStep(front);
    FT_SetOutputCounter(front);

	FT_TimeControlFilter(front);

   for (;;)
   {
		/* Propagating interface for time step dt */

		FT_Propagate(front);
		FT_AddTimeStepToCounter(front);

		//Next time step determined by maximum speed of previous
		//step, assuming the propagation is hyperbolic and
		//is not dependent on second order derivatives of
		//the interface such as curvature, and etc.

		FT_SetTimeStep(front);

		if (FT_IsSaveTime(front))
			FT_Save(front);
		if (FT_IsDrawTime(front))
			FT_Draw(front);

		if (FT_TimeLimitReached(front))
		{
		   FT_PrintTimeStamp(front);
		   break;
		}

	    /* Output section, next dt may be modified */

	    FT_TimeControlFilter(front);
        FT_PrintTimeStamp(front);
     }
     (void) delete_interface(front->interf);
} /* end test_propagate */


static int rotation_vel(
	POINTER params, // void* pointer to the input data structure
	Front *front,
	POINT *p, // function is applied pointwise.  The next two inputs are unused
	HYPER_SURF_ELEMENT *hse, // but required for the rotation function to
	HYPER_SURF *hs,          // have the proper footprint.
	double *vel)	// vel to be updated with x,y velocity of point p
{
	ROTATION_VEL_PARAMS *rv_params = (ROTATION_VEL_PARAMS*)params; // recast to data structure
	double *coords = Coords(p); // vel
	double V,xcomp,ycomp;
	double rad;
	double *cen = rv_params->cen;
	double omega_0 = rv_params->omega_0;
	double dx,dy;

	dx = coords[0] - cen[0]; 
	dy = coords[1] - cen[1];

	rad = sqrt(sqr(dx) + sqr(dy));
	if (rad == 0.0)
    {
        vel[0] = vel[1] = 0.0;
        return 1;
    }

	xcomp = fabs(coords[1]-cen[0])/rad;
    ycomp = fabs(coords[0]-cen[1])/rad;
    V = rad*(omega_0);
    if (coords[0]-cen[0] >= 0.0 &&
	coords[1]-cen[1] >= 0.0) /*1st quadrant*/
    {
        vel[0] = -V*xcomp;
        vel[1] =  V*ycomp;
    }
    else if (coords[0]-cen[0] <= 0.0 &&
	coords[1]-cen[1] >= 0.0) /*2nd quadrant*/
    {
        vel[0] = -V*xcomp;
        vel[1] = -V*ycomp;
    }
    else if (coords[0]-cen[0] <= 0.0 &&
	coords[1]-cen[1] <= 0.0) /*3rd quadrant*/
    {
        vel[0] =  V*xcomp;
        vel[1] = -V*ycomp;
    }
    else if (coords[0]-cen[0] >= 0.0 &&
	coords[1]-cen[1] <= 0.0) /*4th quadrant*/
    {
        vel[0] =  V*xcomp;
        vel[1] =  V*ycomp;
    }
}	/* end rotation_vel */
