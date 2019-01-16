#include <front/fdecs.h>		

LOCAL void cell3d_case01(struct _CELL_INFO_3D*);
LOCAL void cell3d_case02(struct _CELL_INFO_3D*);
LOCAL void cell3d_case03(struct _CELL_INFO_3D*);
LOCAL void cell3d_case04(struct _CELL_INFO_3D*);
LOCAL void cell3d_case05(struct _CELL_INFO_3D*);
LOCAL void cell3d_case06(struct _CELL_INFO_3D*);
LOCAL void cell3d_case07(struct _CELL_INFO_3D*);
LOCAL void cell3d_case08(struct _CELL_INFO_3D*);
LOCAL void cell3d_case09(struct _CELL_INFO_3D*);
LOCAL void cell3d_case10(struct _CELL_INFO_3D*);
LOCAL void cell3d_case11(struct _CELL_INFO_3D*);
LOCAL void cell3d_case12(struct _CELL_INFO_3D*);
LOCAL void cell3d_case13(struct _CELL_INFO_3D*);
LOCAL void cell3d_case14(struct _CELL_INFO_3D*);

LOCAL void set_prime_components(COMPONENT****);
LOCAL void copy_cell3d(const CELL_INFO_3D*, CELL_INFO_3D*);
LOCAL void copy_back_cell3d(const struct _CELL_INFO_3D*,
			struct _CELL_INFO_3D*);
LOCAL int compare_comp(COMPONENT (*)[2][2], COMPONENT****,int);

LOCAL void rot24(CELL_INFO_3D *cell3d, int);
LOCAL void x_rotation(struct _CELL_INFO_3D *cell3d);
LOCAL void y_rotation(struct _CELL_INFO_3D *cell3d);
LOCAL void z_rotation(struct _CELL_INFO_3D *cell3d);
LOCAL void reverse24(struct _CELL_INFO_3D*,int);
LOCAL void x_reverse(struct _CELL_INFO_3D*);
LOCAL void y_reverse(struct _CELL_INFO_3D*);
LOCAL void z_reverse(struct _CELL_INFO_3D*);


/*14 cases' volume functions*/

/*case 1, plan case volume functions*/
LOCAL double volume_plane_case01(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);
LOCAL double case01_complement(struct _CELL_INFO_3D *cell3d,double*,
		double*,double*,double*);

/*case 2, edge case volume functions*/
LOCAL double volume_edge_case02(struct _CELL_INFO_3D*,double*,double*);
LOCAL double case02_complement(struct _CELL_INFO_3D*, double*,double*,
		double*,double*,double*,double*);

/*case 3, 1 corner case volume function*/
LOCAL double volume_corner_case03(struct _CELL_INFO_3D*,double*);

/*case 4, glider case volume functions*/
LOCAL double volume_glider_case04(struct _CELL_INFO_3D*,double*,
		double*,double*);
LOCAL double case04_complement(struct _CELL_INFO_3D*,double*,
		double*,double*,double*,double*); 

/*case 5, hexagon case volume functions*/
LOCAL double volume_hexagon_case05(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);
LOCAL double case05_complement(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);

/*case 6*/
LOCAL double volume_corner1_case06(struct _CELL_INFO_3D*,double*);
LOCAL double volume_corner2_case06(struct _CELL_INFO_3D*,double*);

/*case 7, ij-twister case volume functions*/
LOCAL double volume_twister_case07(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);
LOCAL double case07_complement(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);

/*case 8*/
LOCAL double volume_corner_case08(struct _CELL_INFO_3D*,double*);
LOCAL double volume_edge_case08(struct _CELL_INFO_3D*,double*,double*);

/*case 9*/
LOCAL double volume_corner2_case09(struct _CELL_INFO_3D*,double*);

/*case 10, ji-twister case volume functions*/
LOCAL	double volume_twister_case10(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);
LOCAL	double case10_complement(struct _CELL_INFO_3D*,double*,
		double*,double*,double*);

/*case 11*/
LOCAL	double volume_edge1_case11(struct _CELL_INFO_3D*,double*,double*);
LOCAL	double volume_edge2_case11(struct _CELL_INFO_3D*,double*,double*);

/*case 12*/
/*case 12 uses the volume functions from case 6 and case9*/

/*case 13*/
/*case 13 uses the volume function from case 9*/
LOCAL 	double volume_glider_case13(struct _CELL_INFO_3D*,
		double*,double*,double*);
/*case 14*/
/*case 14 use the volume functions from case 3,8*/
LOCAL	double volume_corner1_case14(struct _CELL_INFO_3D*,double*);
LOCAL 	double volume_corner4_case14(struct _CELL_INFO_3D*,double*);


/*functions that assign volume fraction to neighbour cell, 
  which is legal for transfer and has largest contact region area */

/*case 1*/
LOCAL void cell3d_area_case01_comp0(struct _CELL_INFO_3D*,double*,
	double*,double*,double*,double); 
LOCAL void cell3d_area_case01_comp1(struct _CELL_INFO_3D*,double*,
	double*,double*,double*,double);

/*case 2*/
LOCAL void cell3d_area_case02_comp1(struct _CELL_INFO_3D*, double*,
	double*, double);
LOCAL void cell3d_area_case02_comp0(struct _CELL_INFO_3D*, double*,double*,
	double*,double*,double*,double*,double*,double*,double);

/*case 3*/
LOCAL void cell3d_area_case03_comp1(struct _CELL_INFO_3D*,double*,double);
LOCAL void cell3d_area_case03_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double*,double);

/*case 4*/
LOCAL void cell3d_area_case04_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double);
LOCAL void cell3d_area_case04_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double);

/*case 5*/
LOCAL void cell3d_area_case05_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double);
LOCAL void cell3d_area_case05_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double);

/*case 6*/
LOCAL void cell3d_area_case06_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);
LOCAL void cell3d_area_case06_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);

/*case 7*/
LOCAL void cell3d_area_case07_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double);
LOCAL void cell3d_area_case07_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double);

/*case 8*/
LOCAL void cell3d_area_case08_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double);
LOCAL void cell3d_area_case08_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double*,double);

/*case 9*/
LOCAL void cell3d_area_case09_comp1(struct _CELL_INFO_3D*,double*,double*,
	double);
LOCAL void cell3d_area_case09_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);

/*case 10*/
LOCAL void cell3d_area_case10_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double*,double);
LOCAL void cell3d_area_case10_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double);

/*case 11*/
LOCAL	void cell3d_area_case11_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double);
LOCAL	void cell3d_area_case11_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double*,double);

/*case 12*/
LOCAL	void cell3d_area_case12_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double);
LOCAL	void cell3d_area_case12_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);

/*case 13*/
LOCAL	void cell3d_area_case13_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);
LOCAL	void cell3d_area_case13_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double*,double*,double*,double);

/*case 14*/
LOCAL	void cell3d_area_case14_comp1(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double);
LOCAL	void cell3d_area_case14_comp0(struct _CELL_INFO_3D*,double*,double*,
	double*,double*,double*,double);

LOCAL const double* crx_in_idir(const struct _CELL_INFO_3D*,int,int);
LOCAL const double* crx_in_jdir(const struct _CELL_INFO_3D*,int,int);
LOCAL const double* crx_in_kdir(const struct _CELL_INFO_3D*,int,int);
LOCAL const double* crx_in_between(const int*,const int*, const struct _CELL_INFO_3D*);

EXPORT void cell_volume(
	CELL_INFO_3D* cell3d,
	int ic,
	int *top_comp)
{
	int i,j,k;
	int num_rots;
	boolean case_found;
	static COMPONENT ****prime_comp;
	CELL_INFO_3D cell3d_rot;
	void (*cell3d_comp2[14])(CELL_INFO_3D*) =
	{
	    cell3d_case01,
	    cell3d_case02,
	    cell3d_case03,
	    cell3d_case04,
	    cell3d_case05,
	    cell3d_case06,
	    cell3d_case07,
	    cell3d_case08,
	    cell3d_case09,
	    cell3d_case10,
	    cell3d_case11,
	    cell3d_case12,
	    cell3d_case13,
	    cell3d_case14,
	};
	if (prime_comp == NULL)
	{
	    quad_array(&prime_comp,14,2,2,2,sizeof(COMPONENT));
	    set_prime_components(prime_comp);
	}
	copy_cell3d(cell3d,&cell3d_rot);
	case_found = NO;
	for (j = 0; j <= 24; ++j)
	{
	    for (i = 0; i < 14; ++i)
	    {
	    	if (compare_comp(cell3d_rot.comp,prime_comp,i))
		{
		    case_found = YES;
		    cell3d_comp2[i](&cell3d_rot);
		    break;
		}
	    }
	    if (case_found == YES) break;
	    rot24(&cell3d_rot,j);
	}
	num_rots = j;

	if (case_found == NO)
	{
	    printf("ERROR: in cell_volume()  NO case is found\n");
	    clean_up(ERROR);
	}
	
	for (j = num_rots; j > 0; --j)
	{
	    reverse24(&cell3d_rot,j);	
	}
	copy_back_cell3d(&cell3d_rot,cell3d);
}


LOCAL 	void cell3d_case01(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double	(*crds)[2][2][3]; 
	double c1[3],c2[3],c3[3],c4[4];
	double c5[3],c6[3],c7[3],c8[3];
	
	/*plane case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	/*the volume you want has the same component as c5,c6,c7,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][0][0])
	{
	    vol = volume_plane_case01(cell3d,c7,c8,c5,c6);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case01_comp1(cell3d,c5,c6,c7,c8,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c4*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case01_complement(cell3d,c1,c2,c3,c4);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else	
	    	cell3d_area_case01_comp0(cell3d,c1,c2,c3,c4,vol);
	}
	else
	{
	    printf("Logic Error in case 1\n");
	    clean_up(ERROR);
	}

	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/

	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    	cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES; 
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][0][0])
	        cell3d_area_case01_comp1(cell3d,c5,c6,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case01_comp0(cell3d,c1,c2,c3,c4,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	*/
}

LOCAL 	void cell3d_case02(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double	(*crds)[2][2][3]; 
	double c1[3],c2[3],c3[3],c4[4];
	double c5[3],c6[3],c7[3],c8[3];
	
	/*edge case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	/*the volume you want has the same component as c7,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][1][0])
	{
	    vol = volume_edge_case02(cell3d,c8,c7);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else 
	    	cell3d_area_case02_comp1(cell3d,c7,c8,vol);
	}
	/*the volume has the same component as c1,c2,c3,c4,c5,c6*/
	else if(cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case02_complement(cell3d,c1,c2,c3,c4,c5,c6);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case02_comp0(cell3d,c1,c2,c3,
				c4,c5,c6,c7,c8,vol);
	}
	else 
	{
	    printf("Logic error in case 2\n");
	    clean_up(ERROR);
	}
	
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for 
	the stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][0])
	        cell3d_area_case02_comp1(cell3d,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case02_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,c7,c8,vol);
	    //cell3d->is_cell_unstalbe = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL 	void cell3d_case03(CELL_INFO_3D *cell3d)
{
	boolean met_unstable_cell = NO;
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double	(*crds)[2][2][3]; 
	double c1[3],c2[3],c3[3],c4[4];
	double c5[3],c6[3],c7[3],c8[3];
	
	/*one corner case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = volume_corner_case03(cell3d,c8);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case03_comp1(cell3d,c8,vol);
	}
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume-volume_corner_case03(cell3d,c8);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case03_comp0(cell3d,c1,c2,c3,
			c4,c5,c6,c7,c8,vol);
	}
	else
	{
	    printf("Logic error in case 3\n");
	    clean_up(ERROR);
	}

	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    		,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case03_comp1(cell3d,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case03_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,c7,c8,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL 	void cell3d_case04(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double	(*crds)[2][2][3]; 
	double c1[3],c2[3],c3[3],c4[4];
	double c5[3],c6[3],c7[3],c8[3];
	
	/*glider case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = volume_glider_case04(cell3d,c6,c8,c7);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else 
	    	cell3d_area_case04_comp1(cell3d,c5,c6,c7,c8,vol);
	}
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case04_complement(cell3d,c1,c2,c3,c4,c5);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else 
	    	cell3d_area_case04_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,c7,vol);
	}
	else
	{
	   printf("Logic error in case 4\n");
	   clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/

	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    			,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case04_comp1(cell3d,c5,c6,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case04_comp0(cell3d,c1,c2,c3,c4,c5,c6,c7,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP;
	    
	}
	*/
}

LOCAL 	void cell3d_case05(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double	(*crds)[2][2][3]; 
	double c1[3],c2[3],c3[3],c4[4];
	double c5[3],c6[3],c7[3],c8[3];
	
	/*hexagon case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	/*the volume has the same component as c4,c6,c7,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = volume_hexagon_case05(cell3d,c4,c6,c7,c8); 
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	     else
	     	cell3d_area_case05_comp1(cell3d,c1,c2,c3,c4,
				c5,c6,c7,vol);
	}
	/*the volume has the same component as c1,c2,c3,c5*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case05_complement(cell3d,c1,c2,c3,c5);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case05_comp0(cell3d,c1,c2,c3,c4,
				c5,c6,c7,vol);
	}
	else 
	{
	    printf("Logic error in case 5\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case05_comp1(cell3d,c1,c2,c3,c4,c5,c6,c7,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case05_comp0(cell3d,c1,c2,c3,c4,c5,c6,c7,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}
LOCAL 	void cell3d_case06(CELL_INFO_3D *cell3d)
{
	double vol_1,vol_2,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	/*two corners sharing the same face case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	
	vol_1 = volume_corner1_case06(cell3d,c7);
	vol_2 = volume_corner2_case06(cell3d,c6);
	/*the volume you want has the same component as c6,c7*/
	if (cell3d->soln_comp == cell3d->comp[1][0][1] &&
	    cell3d->soln_comp == cell3d->comp[1][1][0])
	{
	    vol = vol_1+vol_2; 
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case06_comp1(cell3d,c1,c2,c4,c6,c7,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c4,c5,c8*/
	else if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = cell_volume-(vol_1+vol_2);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case06_comp0(cell3d,c1,c2,c4,c6,c7,vol);
	}
	else 
	{
	    printf("Logic error in case 6\n");
	    clean_up(ERROR);
	}

	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    		,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][0][1] &&
	    	cell3d->soln_comp == cell3d->comp[1][1][0])
	        cell3d_area_case06_comp1(cell3d,c1,c2,c4,c6,c7,vol);
	    else if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case06_comp0(cell3d,c1,c2,c4,c6,c7,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP;
	}
	    */
}

LOCAL 	void cell3d_case07(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	/*ij-twister case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	/*the volume you want has the same component as c4,c5,c6,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = volume_twister_case07(cell3d,c5,c6,c8,c4);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else	
	        cell3d_area_case07_comp1(cell3d,c1,c2,c4,
			c5,c6,c7,c8,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c7*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case07_complement(cell3d,c1,c2,c3,c7);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else 
	        cell3d_area_case07_comp0(cell3d,c1,c2,c3,
			c4,c5,c7,vol);
	}
	else
	{
	    printf("Logic error in case 7\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    		,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case07_comp1(cell3d,c1,c2,c4,c5,c6,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case07_comp0(cell3d,c1,c2,c3,c4,c5,c7,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	*/
}

LOCAL 	void cell3d_case08(CELL_INFO_3D *cell3d)
{
	double vol_edge,vol_corner,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	/*edge with corner case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	vol_corner = volume_corner_case08(cell3d,c5);
	vol_edge = volume_edge_case08(cell3d,c4,c8);

	/*the volume you want has the same component as C4,C8,C5*/
	if (cell3d->soln_comp == cell3d->comp[1][0][0])
	{
	    vol = vol_corner+vol_edge;
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else 
	    	cell3d_area_case08_comp1(cell3d,c4,c5,c8,vol);
	}
	/*the volume you want has the same component as C1,C2,C3,C6,C7*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume-(vol_corner+vol_edge);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case08_comp0(cell3d,c1,c2,c3,c4,c5,c6,c7,c8,vol);
	}
	else
	{
	    printf("Logic error in case 8\n");
	    clean_up(ERROR);
	}
	
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][0][0])
	        cell3d_area_case08_comp1(cell3d,c4,c5,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case08_comp0(cell3d,c1,c2,c3,c4,c5,c6,c7,c8,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL 	void cell3d_case09(CELL_INFO_3D *cell3d)
{
	double vol_1,vol_2,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	/*two corners on c4,c5*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	vol_1 = volume_corner_case08(cell3d,c5);
	vol_2 = volume_corner2_case09(cell3d,c4);
	/*the volume you want has the same component as c4 and c5*/
	if (cell3d->soln_comp == cell3d->comp[1][0][0])
	{
	    vol = vol_1 + vol_2;
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case09_comp1(cell3d,c4,c5,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c6,c7,c8*/
	else if(cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume - (vol_1 + vol_2);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case09_comp0(cell3d,c1,c2,c3,c4,c5,vol);
	}
	else
	{
	    printf("Logic error in case 9\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = NO;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    		,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][0][0])
	        cell3d_area_case09_comp1(cell3d,c4,c5,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case09_comp0(cell3d,c1,c2,c3,c4,c5,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	
	    */
}

LOCAL 	void cell3d_case10(CELL_INFO_3D *cell3d)
{
	double vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	/*ji-twister case*/
	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	/*the volume you want has the same component as c4,c5,c7,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = volume_twister_case10(cell3d,c4,c8,c7,c5);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else 
	        cell3d_area_case10_comp1(cell3d,c1,c2,c3,
			c4,c5,c6,c7,c8,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c6*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = case10_complement(cell3d,c1,c2,c3,c6);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case10_comp0(cell3d,c1,c2,c3,
			c4,c5,c6,vol);
	}
	else
	{
	    printf("Logic error in case 10\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n"
	    		,cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case10_comp1(cell3d,c1,c2,c3,c4,
			c5,c6,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case10_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL 	void cell3d_case11(CELL_INFO_3D *cell3d)
{
	double vol1,vol2,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	
	vol1 = volume_edge1_case11(cell3d,c3,c4);
	vol2 = volume_edge2_case11(cell3d,c5,c6);
	/*the volume you want has the same component as c3,c4,c5,c6*/
	if (cell3d->soln_comp == cell3d->comp[0][1][0])
	{
	    vol = vol1+vol2;
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case11_comp1(cell3d,c3,c4,c5,c6,vol);
	}
	/*the volume you want has the same component as c1,c2,c7,c8*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume -(vol1 + vol2);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case11_comp0(cell3d,c1,c2,c3,
				c4,c5,c6,c7,c8,vol);
	}
	else
	{
	    printf("Logic error in case 11\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[0][1][0])
	        cell3d_area_case11_comp1(cell3d,c3,c4,
			c5,c6,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case11_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,c7,c8,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL 	void cell3d_case12(CELL_INFO_3D *cell3d)
{
	double vol_1,vol_2,vol_3,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	vol_1 = volume_corner1_case06(cell3d,c7);
	vol_2 = volume_corner2_case09(cell3d,c4);
	vol_3 = volume_corner2_case06(cell3d,c6);
	/*the volume you want has the same componet as c4,c6,c7*/
	if (cell3d->soln_comp == cell3d->comp[0][1][1])
	{
	    vol = vol_1 + vol_2 + vol_3;
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case12_comp1(cell3d,c4,c6,c7,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c5,c8*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume-(vol_1+vol_2+vol_3);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case12_comp0(cell3d,c1,c2,c4,c6,c7,vol);
	}
	else
	{
	    printf("Logic error in case 12\n");
	    clean_up(ERROR);
	}

	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[0][1][1])
	        cell3d_area_case12_comp1(cell3d,c4,c6,c7,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case12_comp0(cell3d,c1,c2,c4,c6,c7,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL void cell3d_case13(CELL_INFO_3D *cell3d)
{
	double vol_glider,vol_corner,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}
	vol_glider = volume_glider_case13(cell3d,c5,c6,c7);
	vol_corner = volume_corner2_case09(cell3d,c4);

	/*the volume you want has the same component as c4,c5,c6,c7*/
	if (cell3d->soln_comp == cell3d->comp[0][1][1])
	{
	    vol = vol_glider + vol_corner;
	    if (cell3d->cell_comp == cell3d->soln_comp)
	    	cell3d->nb_frac[1][1][1] = vol;
	    else
	       cell3d_area_case13_comp1(cell3d,c4,c5,c6,c7,c8,vol);
	}
	/*the volume you want has the same component as c1,c2,c3,c8*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume-(vol_glider+vol_corner);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else
	    	cell3d_area_case13_comp0(cell3d,c1,c2,
			c3,c4,c5,c6,c7,c8,vol);
	}
	else
	{    
	    printf("Logic error in case 13\n");
	    clean_up(ERROR);
	}
	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[0][1][1])
	        cell3d_area_case13_comp1(cell3d,c4,c5,c6,c7,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case13_comp0(cell3d,c1,c2,c3,c4,
			c5,c6,c7,c8,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}

LOCAL void cell3d_case14(CELL_INFO_3D *cell3d)
{
	double vol_1,vol_2,vol_3,vol_4,vol,cell_volume;
	int i;
	int (*ix)[2][2],(*iy)[2][2],(*iz)[2][2];
	double (*crds)[2][2][3];
	double c1[3],c2[3],c3[3],c4[3];
	double c5[3],c6[3],c7[3],c8[3];

	ix = cell3d->ix;
	iy = cell3d->iy;
	iz = cell3d->iz;
	crds = cell3d->crds;
	cell_volume = cell3d->full_cell_vol;
	for (i = 0; i < 3; ++i)
	{
	    c1[i] = crds[ix[0][0][0]][iy[0][0][0]][iz[0][0][0]][i];
	    c2[i] = crds[ix[0][0][1]][iy[0][0][1]][iz[0][0][1]][i];
	    c3[i] = crds[ix[0][1][0]][iy[0][1][0]][iz[0][1][0]][i];
	    c4[i] = crds[ix[0][1][1]][iy[0][1][1]][iz[0][1][1]][i];
	    c5[i] = crds[ix[1][0][0]][iy[1][0][0]][iz[1][0][0]][i];
	    c6[i] = crds[ix[1][0][1]][iy[1][0][1]][iz[1][0][1]][i];
	    c7[i] = crds[ix[1][1][0]][iy[1][1][0]][iz[1][1][0]][i];
	    c8[i] = crds[ix[1][1][1]][iy[1][1][1]][iz[1][1][1]][i];
	}

	vol_1 = volume_corner1_case14(cell3d,c3);
	vol_2 = volume_corner_case03(cell3d,c8);
	vol_3 = volume_corner_case08(cell3d,c5);
	vol_4 = volume_corner4_case14(cell3d,c2);
	
	/*the volume you want has the same component as c2,c3,c5,c8*/
	if (cell3d->soln_comp == cell3d->comp[1][1][1])
	{
	    vol = vol_1 + vol_2 + vol_3 + vol_4;
	    if (cell3d->cell_comp == cell3d->soln_comp)
		cell3d->nb_frac[1][1][1] = vol;
	    else
	        cell3d_area_case14_comp1(cell3d,c2,c3,c5,c8,vol);
	
	}
	/*the volume you want has the same component as c1,c4,c6,c6*/
	else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	{
	    vol = cell_volume - (vol_1+vol_2+vol_3+vol_4);
	    if (cell3d->cell_comp == cell3d->soln_comp)
	        cell3d->nb_frac[1][1][1] = vol;
	    else 
	    	cell3d_area_case14_comp0(cell3d,c1,c2,c3,c5,c8,vol);
	}
	else
	{
	    printf("Logic error in case 14\n");
	    clean_up(ERROR);
	}

	/* 
	The lower bound for eta_new in Csolve3d has to be >= 0.49 for the 
	stability consideration. If this cell center is SOLUTE_COMP but 
	the volume fraction for the solute component is less than 0.49 or 
	volume of solute < 0.49*full_cell_volume, then our algorithm will send 
	the solute component to one of its legal neighbour or generate 
	an orphan, also this cell center will be transfer to SWITCHED_COMP and 
	this cell will not enter the solver.
	*/
	/*
	if (cell3d->is_soln_solute == YES && 
	    cell3d->is_center_solute == YES &&
	    cell3d->nb_frac[1][1][1] > 0.0 &&
	    cell3d->nb_frac[1][1][1] < UNSTABLE_CELL_FRAC*cell_volume)
	{
	    met_unstable_cell = YES;
	    ++(*unstable_cells);
	    printf("Cell unstable, with vf = %.5f\n",
	    		cell3d->nb_frac[1][1][1]/cell_volume);
	    cell3d->unstable_cell = YES;
	    cell3d->solute_volume = cell3d->nb_frac[1][1][1];
	    cell3d->nb_frac[1][1][1] = UNSTABLE_CELL_FRAC*cell_volume;
	    */
	    /*
	    if (cell3d->soln_comp == cell3d->comp[1][1][1])
	        cell3d_area_case14_comp1(cell3d,c2,c3,c5,c8,vol);
	    else if (cell3d->soln_comp == cell3d->comp[0][0][0])
	        cell3d_area_case14_comp0(cell3d,c1,c2,c3,c5,c8,vol);
	    //cell3d->is_cell_unstable = YES;
	    top_comp[ic] = SWITCHED_COMP; 
	}
	    */
}


LOCAL double volume_plane_case01(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4;
	double D1[3], D2[3], D3[3], D4[3], D5[3], D6[3], D7[3], D8[3],
	      D9[3], D10[3], D11[3], D12[3], D13[3], D14[3], D15[3],
	      D16[3], D17[3], D18[3], D19[3]; 
	double A1[3], A2[3], A3[3], A4[3];
	double V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12;
	double volume;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_idir(cell3d,0,1);
		
	for (i = 0; i < 3; i++)
        {       
	    D1[i]=crx1[i]-crx2[i];
	    D2[i]=crx1[i]-crx3[i];
	    D3[i]=crx1[i]-crx4[i];
	    D4[i]=crx2[i]-crx3[i];
	    D5[i]=crx2[i]-crx4[i];
	    D6[i]=crx3[i]-crx4[i];
	    D7[i]=crx1[i]-crn3[i];	/* crx1 - c5   */
	    D8[i]=crx1[i]-crn4[i];	/* crx1 - c6  */
	    D9[i]=crx2[i]-crn3[i];	/* crx2 - c5  */
	    D10[i]=crx3[i]-crn3[i];	/* crx3 - c5  */
	    D11[i]=crx4[i]-crn3[i];	/* crx4 - c5  */
	    D12[i]=crx4[i]-crn4[i];	/* crx4 - c6  */
	    D13[i]=crx4[i]-crn2[i];	/* crx4 - c8   */
	    D14[i]=crn3[i]-crn4[i];	/* c5 - c6   */
	    D15[i]=crn3[i]-crn1[i];	/* c5 - c7   */
	    D16[i]=crn3[i]-crn2[i];	/* c5 - c8   */
	    D17[i]=crn4[i]-crn1[i];	/* c6 - c7   */
	    D18[i]=crn4[i]-crn2[i];	/* c6 - c8   */
	    D19[i]=crn1[i]-crn2[i];	/* c7 - c8   */
	}

	/* triangulation given by diagonal crx2-crx3  */
	V7 = fabs(Det3d(D1,D4,D10));	/* tetrahedron crx1-crx2-crx3-c5  */
	V8 = fabs(Det3d(D4,D6,D11));	/* tetrahedron crx2-crx3-crx4-c5  */
	V9 = fabs(Det3d(D5,D11,D14));	/* tetrahedron crx2-crx4-c5-c6    */
	V10 = fabs(Det3d(D9,D14,D18));	/* tetrahedron crx2-c5-c6-c8      */
	V11 = fabs(Det3d(D1,D9,D16));	/* tetrahedron crx1-crx2-c5-c8    */
	V12 = fabs(Det3d(D7,D15,D19));	/* tetrahedron crx1-c5-c7-c8      */
	
	volume = (V7+V8+V9+V10+V11+V12)/6.0; 
	
	return volume;
}	/* end volume_plane_case01 */

LOCAL   double case01_complement(
		CELL_INFO_3D *cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn4[])
{
	int i; 
	const double *crx1,*crx2,*crx3,*crx4;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],
		E7[3],E8[3],E9[3],E10[3],E11[3];
	double V1,V2,V3,V4,V5,V6;
	double volume;

	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_idir(cell3d,1,1);
	crx3 = crx_in_idir(cell3d,0,0);
	crx4 = crx_in_idir(cell3d,0,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C3->P1*/
	    E1[i] = crx1[i]-crn3[i];
	    /*C3->P2*/
	    E2[i] = crx2[i]-crn3[i];
	    /*C3->P3*/
	    E3[i] = crx3[i]-crn3[i];
	    /*C3->C4*/
	    E4[i] = crn4[i]-crn3[i];
	    /*C1->C4*/
	    E5[i] = crn4[i]-crn1[i];
	    /*P3->C4*/
	    E6[i] = crn4[i]-crx3[i];
	    /*C2->P2*/
	    E7[i] = crx2[i]-crn2[i];
	    /*C2->P3*/
	    E8[i] = crx3[i]-crn2[i];
	    /*C2->P4*/
	    E9[i] = crx4[i]-crn2[i];
	    /*C2->C4*/
	    E10[i] = crn4[i]-crn2[i];
	    /*C2->C1*/
	    E11[i] = crn1[i]-crn2[i];
	}
	/* volume of C3P1P2P3 */
	V1 = fabs(Det3d(E1,E2,E3));
	/* volume of C4C3P2P3 */
	V2 = fabs(Det3d(E2,E3,E4));
	/* volume of C1C4C3P3 */
	V3 = fabs(Det3d(E4,E5,E6));
	/* volume of C2P2P3P4 */
	V4 = fabs(Det3d(E7,E8,E9));
	/* volume of C4C2P2P3 */
	V5 = fabs(Det3d(E7,E8,E10));
	/* volume of C1C2C4P3 */
	V6 = fabs(Det3d(E8,E10,E11));

	volume = 1.0/6*(V1+V2+V3+V4+V5+V6);
	return volume;
}

LOCAL double volume_edge_case02(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4;
	double D1[3], D2[3], D3[3], D4[3], D5[3], D6[3], 
		D7[3], D8[3],D9[3], D10[3];
	double V4, V5, V6;
	double volume;
	
	crx1 = crx_in_idir(cell3d,1,1);  
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,1,0);
	crx4 = crx_in_jdir(cell3d,0,1);
	
	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
            D2[i] = crx1[i] - crx3[i];
 	    D3[i] = crx2[i] - crx3[i];
	    D4[i] = crx2[i] - crx4[i];
	    D5[i] = crx3[i] - crx4[i];
	    D6[i] = crx2[i] - crn2[i];	/* crx2 - c7  */
	    D7[i] = crx3[i] - crn1[i];  /* crx3 - c8  */
	    D8[i] = crx4[i] - crn1[i];  /* crx4 - c8  */
	    D9[i] = crx4[i] - crn2[i];	/* crx4 - c7  */
	    D10[i] = crn2[i] - crn1[i]; /* c7 -c8    */
	}
	/* triangulation given by diagonal crx2-crx3  */
	/* tetrahedron crx1-crx2-crx3-c8 */
	V4 = fabs(Det3d(D1,D3,D7)); 		
	/* tetrahedron crx2-crx3-crx4-c8 */
	V5 = fabs(Det3d(D3,D5,D8));	
	/* tetrahedron crx3-crx4-c7-c8 */
	V6 = fabs(Det3d(D5,D9,D10));	
	
	volume = (V4+V5+V6)/6.0;
	return volume;
}	/* end volume_edge_case02 */

LOCAL double case02_complement(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[])
{
	int i;
	const double *crx1,*crx2,*crx3,*crx4;
	double E1[3],E2[3],E3[3],E4[3],E5[3],
		E6[3],E7[3],E8[3],E9[3],E10[3];
	double V1,V2,V3,V4,V5,V6,V7,V8,V9;
	double volume;

	crx1 = crx_in_idir(cell3d,1,1);
	crx2 = crx_in_jdir(cell3d,1,1);
	crx3 = crx_in_idir(cell3d,1,0);
	crx4 = crx_in_jdir(cell3d,0,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C2->P2*/
	    E1[i] = crx2[i]-crn2[i];
	    /*C2->P3*/
	    E2[i] = crx3[i]-crn2[i];
	    /*C2->P4*/
	    E3[i] = crx4[i]-crn2[i];
	    /*C2->P1*/
	    E4[i] = crx1[i]-crn2[i];
	    /*C2->C1*/
	    E5[i] = crn1[i]-crn2[i];
	    /*C2->C6*/
	    E6[i] = crn6[i]-crn2[i];
	    /*C2->C5*/
	    E7[i] = crn5[i]-crn2[i];
	    /*C2->C4*/
	    E8[i] = crn4[i]-crn2[i];
	    /*C2->C3*/
	    E9[i] = crn3[i]-crn2[i];
	    /*C2->C1*/
	    E10[i] = crn1[i]-crn2[i];
	}
	
	/*volume of C2P2P3P4*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C2P1P2P3*/
	V2 = fabs(Det3d(E1,E2,E4));
	/*volume of C2P3P4C1*/
	V3 = fabs(Det3d(E2,E3,E5));
	/*volume of C6P2P4C2*/
	V4 = fabs(Det3d(E1,E3,E6));
	/*volume of C5C6C2P4*/
	V5 = fabs(Det3d(E3,E6,E7));
	/*volume of C1C2C5P4*/
	V6 = fabs(Det3d(E5,E7,E3));
	/*volume of C4C2P1P3*/
	V7 = fabs(Det3d(E2,E4,E8));
	/*volume of C3C4P3C2*/
	V8 = fabs(Det3d(E2,E8,E9));
	/*volume of C1C2C3P3*/
	V9 = fabs(Det3d(E2,E9,E10));
	
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7+V8+V9);
	return volume;
}

LOCAL double volume_corner_case03(
	CELL_INFO_3D *cell3d,
	double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;

        crx1 = crx_in_idir(cell3d,1,1); 
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_kdir(cell3d,1,1); 

	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c8 */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}	/* end volume_corner_case03 */

LOCAL double volume_glider_case04(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4, *crx5;
	double E1[3], E2[3], E3[3], E4[3], E5[3], E6[3], E7[3];
	double V1, V2, V3, V4, V5, V6, V7; 
	double volume;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
       	crx3 = crx_in_jdir(cell3d,0,1); 
	crx4 = crx_in_kdir(cell3d,1,0); 
	crx5 = crx_in_idir(cell3d,0,1); 
	
	for (i = 0; i < 3; ++i)
	{
	    E1[i] = crx4[i]-crx2[i]; 
	    E2[i] = crx3[i]-crx2[i];
	    /* c8->p2 */
	    E3[i] = crn2[i]-crx2[i];
	    /* c6->p2 */
	    E4[i] = crn1[i]-crx2[i];
	    E5[i] = crx5[i]-crx2[i];
	    /* c7->p2 */
	    E6[i] = crn3[i]-crx2[i];
	    E7[i] = crx1[i]-crx2[i];
	}

	/* volume for P2P3P4C8 */
	V1 = fabs(Det3d(E1,E2,E3));
	/* volume for P2P4C6C8 */
	V2 = fabs(Det3d(E3,E1,E4));
	/* volume for P2P4P5C6 */
	V3 = fabs(Det3d(E1,E4,E5));
	/* volume for P2P3C7C8 */
	V4 = fabs(Det3d(E2,E3,E6));
	/* volume for P2P3P1C7 */
	V5 = fabs(Det3d(E2,E6,E7));
	volume = 1.0/6*(V1+V2+V3+V4+V5);
	return volume; 
}

LOCAL	double case04_complement(
		CELL_INFO_3D *cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn4[],
		double crn5[]) 
{
	int i; 
	const double *crx1,*crx2,*crx3,*crx4,*crx5;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],
		E7[3],E8[3],E9[3],E10[3],E11[3],E12[3],
		E13[3],E14[3],E15[3];
	double V1,V2,V3,V4,V5,V6,V7,V8,V9;
	double volume;

	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_idir(cell3d,1,1);
	crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_kdir(cell3d,1,0);
	crx5 = crx_in_idir(cell3d,0,1);
	
	for (i = 0; i < 3; ++i)
	{
	    /* C4->P2 */
	    E1[i] = crx2[i]-crn4[i];
	    /* C4->P3 */
	    E2[i] = crx3[i]-crn4[i];
	    /* C4->P4 */
	    E3[i] = crx4[i]-crn4[i];
	    /* C4->C1 */
	    E4[i] = crn1[i]-crn4[i];
	    /* C1->C5 */
	    E5[i] = crn5[i]-crn1[i];
	    /* C1->P3 */
	    E6[i] = crx3[i]-crn1[i];
	    /* C1->P4 */ 
	    E7[i] = crx4[i]-crn1[i];
	    /* C4->C3 */
	    E8[i] = crn3[i]-crn4[i];
	    /* C3->P1 */
	    E9[i] = crx1[i]-crn3[i];
	    /* C3->P2 */
	    E10[i] = crx2[i]-crn3[i];
	    /* C3->P3 */
	    E11[i] = crx3[i]-crn3[i];
	    /* C4->C2 */
	    E12[i] = crn2[i]-crn4[i];
	    /* C2->P2 */
	    E13[i] = crx2[i]-crn2[i];
	    /* C2->P4 */
	    E14[i] = crx4[i]-crn2[i];
	    /* C2->P5 */
	    E15[i] = crx5[i]-crn2[i];
	}

	/*volume of C4P2P3P4*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C1C4P3P4*/
	V2 = fabs(Det3d(E2,E3,E4));
	/*volume of C5C1P3P4*/
	V3 = fabs(Det3d(E5,E6,E7));
	/*volume of C3C4C1P3*/
	V4 = fabs(Det3d(E2,E4,E8));
	/*volume of C3C4P2P3*/
	V5 = fabs(Det3d(E1,E2,E8));
	/*volume of C3P1P2P3*/
	V6 = fabs(Det3d(E9,E10,E11));
	/*volume of C2C4C1P4*/
	V7 = fabs(Det3d(E3,E4,E12));
	/*volume of C2C4P2P4*/
	V8 = fabs(Det3d(E1,E3,E12));
	/*volume of C2P2P4P5*/
	V9 = fabs(Det3d(E13,E14,E15));

	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7+V8+V9);
	return volume;
}

LOCAL double volume_hexagon_case05(
	CELL_INFO_3D *cell3d,
	double crn4[],
	double crn6[],
	double crn7[],
	double crn8[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4, *crx5, *crx6;
        double E1[3], E2[3], E3[3], E4[3], E5[3], E6[3], E7[3], E8[3], 
	      E9[3], E10[3], E11[3];
	double V1,V2,V3,V4,V5,V6,V7;
	double volume;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0);
        crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P1*/
	    E1[i] = crx1[i]-crn8[i];
	    /*C8->P2*/
	    E2[i] = crx2[i]-crn8[i];
	    /*C8->P4*/
	    E3[i] = crx4[i]-crn8[i];
	    /*C8->C7*/
	    E4[i] = crn7[i]-crn8[i];
	    /*C8->P3*/
	    E5[i] = crx3[i]-crn8[i];
	    /*C8->P5*/
	    E6[i] = crx5[i]-crn8[i];
	    /*C8->C6*/
	    E7[i] = crn6[i]-crn8[i];
	    /*C4->P4*/
	    E8[i] = crx4[i]-crn4[i];
	    /*C4->P5*/
	    E9[i] = crx5[i]-crn4[i];
	    /*C4->C8*/
	    E10[i] = crn8[i]-crn4[i];
	    /*C4->P6*/
	    E11[i] = crx6[i]-crn4[i];
	}

	/*volume of C8P1P2P4*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C8C7P1P3*/
	V2 = fabs(Det3d(E1,E4,E5));
	/*volume of P3C8P1P2*/
	V3 = fabs(Det3d(E1,E2,E5));
	/*volume of C8P2P4P5*/
	V4 = fabs(Det3d(E2,E3,E6));
	/*volume of C8C6P2P5*/
	V5 = fabs(Det3d(E2,E6,E7));
	/*volume of C4P4P5C8*/
	V6 = fabs(Det3d(E8,E9,E10));
	/*volume of C4P4P5P6*/
	V7 = fabs(Det3d(E8,E9,E11));
	
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7);
	return volume;
}	/* end volume_hexagon_case05 */

LOCAL 	double case05_complement(
		CELL_INFO_3D *cell3d,
		double crn1[],
		double crn2[],
		double crn4[],
		double crn5[])
{
	int i; 
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],E7[3],
		E8[3],E9[3],E10[3],E11[3],E12[3],E13[3];
	double V1,V2,V3,V4,V5,V6,V7;
	double volume;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0);
        crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C1->P2*/
	    E1[i] = crx2[i]-crn1[i];
	    /*C1->P4*/
	    E2[i] = crx4[i]-crn1[i];
	    /*C1->P5*/
	    E3[i] = crx5[i]-crn1[i];
	    /*C1->C2*/
	    E4[i] = crn2[i]-crn1[i];
	    /*C2->P4*/
	    E5[i] = crx4[i]-crn2[i];
	    /*C2->P5*/
	    E6[i] = crx5[i]-crn2[i];
	    /*C2->P6*/
	    E7[i] = crx6[i]-crn2[i];
	    /*C1->P1*/
	    E8[i] = crx1[i]-crn1[i];
	    /*C1->C4*/
	    E9[i] = crn4[i]-crn1[i];
	    /*C1->C5*/
	    E10[i] = crn5[i]-crn1[i];
	    /*C5->P1*/
	    E11[i] = crx1[i]-crn5[i];
	    /*C5->P2*/
	    E12[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    E13[i] = crx3[i]-crn5[i];
	}
	
	/*volume of C1P2P4P5*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C2C1P4P5*/
	V2 = fabs(Det3d(E2,E3,E4));
	/*volume of C2P4P5P6*/
	V3 = fabs(Det3d(E5,E6,E7));
	/*volume of C1P1P2P4*/
	V4 = fabs(Det3d(E1,E2,E8));
	/*volume of C4P4P1C1*/
	V5 = fabs(Det3d(E2,E8,E9));
	/*volume of C5C1P1P2*/
	V6 = fabs(Det3d(E1,E8,E10));
	/*volume of C5P1P2P3*/
	V7 = fabs(Det3d(E11,E12,E13));
	
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7);
	return volume;
}	/*end of case05_complement*/

LOCAL double volume_corner1_case06(
		CELL_INFO_3D *cell3d,
		double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;
	
	crx1 = crx_in_idir(cell3d,1,0);
        crx2 = crx_in_kdir(cell3d,1,1);
        crx3 = crx_in_jdir(cell3d,0,1);
			 
	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c7  */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}	/* end volume_corner1_case06 */


LOCAL double volume_corner2_case06(
		CELL_INFO_3D *cell3d,
		double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;

	crx1 = crx_in_jdir(cell3d,1,1);
        crx2 = crx_in_idir(cell3d,0,1);
	crx3 = crx_in_kdir(cell3d,1,0);

	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c6  */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}	/* end volume_corner2_case06 */

LOCAL	double volume_twister_case07(
		CELL_INFO_3D *cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn4[])
{
	int i;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],E7[3],
		E8[3],E9[3],E10[3],E11[3],E12[3],E13[3];
	double V1,V2,V3,V4,V5,V6,V7;
	double volume;

	crx1 = crx_in_kdir(cell3d,0,1); 
        crx2 = crx_in_jdir(cell3d,1,0);
        crx3 = crx_in_kdir(cell3d,1,1);
	crx4 = crx_in_idir(cell3d,0,1); 
	crx5 = crx_in_jdir(cell3d,0,1); 
	crx6 = crx_in_idir(cell3d,0,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P2*/
	    E1[i] = crx2[i]-crn3[i];
	    /*C8->P4*/
	    E2[i] = crx4[i]-crn3[i];
	    /*C8->P3*/
	    E3[i] = crx3[i]-crn3[i];
	    /*C8->P1*/
	    E4[i] = crx1[i]-crn3[i];
	    /*C8->C4*/
	    E5[i] = crn4[i]-crn3[i];
	    /*C8->C6*/
	    E6[i] = crn2[i]-crn3[i];
	    /*C6->P3*/
	    E7[i] = crx3[i]-crn2[i];
	    /*C6->P4*/
	    E8[i] = crx4[i]-crn2[i];
	    /*C6->P5*/
	    E9[i] = crx5[i]-crn2[i];
	    /*C6->C5*/
	    E10[i] = crn1[i]-crn2[i];
	    /*P6->P4*/
	    E11[i] = crx4[i]-crx6[i];
	    /*P6->P5*/
	    E12[i] = crx5[i]-crx6[i];
	    /*P6->C5*/
	    E13[i] = crn1[i]-crx6[i];
	}
	/* volume of C8P2P3P4 */
	V1 = fabs(Det3d(E1,E2,E3));
	/* volume of C8P1P2P3 */
	V2 = fabs(Det3d(E1,E3,E4));
	/* volume of C8P1P2C4 */
	V3 = fabs(Det3d(E1,E4,E5));
	/* volume of C8P3P4C6 */
	V4 = fabs(Det3d(E2,E3,E6));
	/* volume of C6P3P4P5 */
	V5 = fabs(Det3d(E7,E8,E9));
	/* volume of C6P4P5C5 */
	V6 = fabs(Det3d(E8,E9,E10));
	/* volume of P6P4P5C5 */
	V7 = fabs(Det3d(E11,E12,E13));

	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7);
	return volume;
}	/* end volume_twister_case07 */



LOCAL 	double case07_complement(
		CELL_INFO_3D *cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn7[])
{
	int i;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],E7[3],
		E8[3],E9[3],E10[3],E11[3],E12[3],E13[3];
	double V1,V2,V3,V4,V5,V6,V7,V8;
	double volume;

        crx1 = crx_in_kdir(cell3d,0,1); 
        crx2 = crx_in_jdir(cell3d,1,0);
        crx3 = crx_in_kdir(cell3d,1,1);
	crx4 = crx_in_idir(cell3d,0,1); 
	crx5 = crx_in_jdir(cell3d,0,1); 
	crx6 = crx_in_idir(cell3d,0,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C1->C2*/
	    E1[i] = crn2[i]-crn1[i];
	    /*C1->P2*/
	    E2[i] = crx2[i]-crn1[i];
	    /*C1->P3*/
	    E3[i] = crx3[i]-crn1[i];
	    /*C1->P4*/
	    E4[i] = crx4[i]-crn1[i];
	    /*C1->P1*/
	    E5[i] = crx1[i]-crn1[i];
	    /*C1->C3*/
	    E6[i] = crn3[i]-crn1[i];
	    /*C1->C7*/
	    E7[i] = crn7[i]-crn1[i];
	    /*C1->P5*/
	    E8[i] = crx5[i]-crn1[i];
	    /*C1->P6*/
	    E9[i] = crx6[i]-crn1[i];
	}
	/*volume of C1C2P2P4*/
	V1 = fabs(Det3d(E1,E2,E4));
	/*volume of C1P2P3P4*/
	V2 = fabs(Det3d(E2,E3,E4));
	/*volume of C1P1P2P3*/
	V3 = fabs(Det3d(E2,E3,E5));
	/*volume of C1P1P3C3*/
	V4 = fabs(Det3d(E3,E5,E6));
	/*volume of C3C7C1P3*/
	V5 = fabs(Det3d(E6,E7,E3));
	/*volume of C1P3P4P5*/
	V6 = fabs(Det3d(E3,E4,E8));
	/*volume of C1P4P5P6*/
	V7 = fabs(Det3d(E4,E8,E9));
	/*volume of C1C7P3P5*/
	V8 = fabs(Det3d(E3,E7,E8));
	
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7+V8);
	return volume;
}	/*end of case07_complement*/

LOCAL	double volume_corner_case08(
		CELL_INFO_3D* cell3d,
		double crn[])
{
	int i;
 	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;

	crx1 = crx_in_jdir(cell3d,0,1); 
        crx2 = crx_in_kdir(cell3d,1,0); 
        crx3 = crx_in_idir(cell3d,0,0);
			 
	
	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c5 */
	}
	
	volume = fabs((Det3d(D1,D2,D3)))/6.0;
	return volume;
}

LOCAL 	double volume_edge_case08(
		CELL_INFO_3D *cell3d,
		double crn4[],
		double crn8[])
{
	int i;
	const double *crx4, *crx5, *crx6, *crx7;
	double D1[3], D2[3], D3[3], D4[3], D5[3];
	double V1, V2, V3;
	double volume;

	crx4 = crx_in_kdir(cell3d,0,1);
        crx5 = crx_in_jdir(cell3d,1,0);
        crx6 = crx_in_kdir(cell3d,1,1);
	crx7 = crx_in_jdir(cell3d,1,1);
		
	for (i = 0; i < 3; ++i)
	{
	    /*C4->P4*/
	    D1[i] = crx4[i]-crn4[i];
	    /*C4->P5*/
	    D2[i] = crx5[i]-crn4[i];
	    /*C4->P6*/
	    D3[i] = crx6[i]-crn4[i];
	    /*C4->P7*/
	    D4[i] = crx7[i]-crn4[i];
	    /*C4->C8*/
	    D5[i] = crn8[i]-crn4[i];
	}

	/*volume of C4P4P5P6*/
	V1 = fabs(Det3d(D1,D2,D3));
	/*volume of C4P5P6P7*/
	V2 = fabs(Det3d(D2,D3,D4));
	/*volume of C4C8P6P7*/
	V3 = fabs(Det3d(D3,D4,D5));
	
	volume = (V1+V2+V3)/6.0;
	
	return volume;
}	/* end volume_edge_case08 */

LOCAL	double volume_corner2_case09(
		CELL_INFO_3D *cell3d,
		double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;
	
	crx1 = crx_in_kdir(cell3d,0,1);
        crx2 = crx_in_jdir(cell3d,1,0);
        crx3 = crx_in_idir(cell3d,1,1);
			 
	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c4  */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}

LOCAL	double volume_twister_case10(
		CELL_INFO_3D* cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn4[])
{
	
	int i;
	const double *crx1, *crx2, *crx3, *crx4, *crx5, *crx6;
        double E1[3], E2[3], E3[3], E4[3], E5[3], E6[3], E7[3], E8[3], 
	      E9[3], E10[3], E11[3], E12[3], E13[3]; 
	double V1,V2,V3,V4,V5,V6,V7;
	double volume;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 
	crx5 = crx_in_kdir(cell3d,0,1); 
	crx6 = crx_in_jdir(cell3d,1,0); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P4*/
	    E1[i] = crx4[i]-crn2[i];
	    /*C8->P1*/
	    E2[i] = crx1[i]-crn2[i];
	    /*C8->P5*/
	    E3[i] = crx5[i]-crn2[i];
	    /*C8->P6*/
	    E4[i] = crx6[i]-crn2[i];
	    /*C8->C4*/
	    E5[i] = crn1[i]-crn2[i];
	    /*C7->C8*/
	    E6[i] = crn2[i]-crn3[i];
	    /*C7->P4*/
	    E7[i] = crx4[i]-crn3[i];
	    /*C7->P1*/
	    E8[i] = crx1[i]-crn3[i];
	    /*C7->P2*/
	    E9[i] = crx2[i]-crn3[i];
	    /*C7->C5*/
	    E10[i] = crn4[i]-crn3[i];
	    /*C5->P1*/
	    E11[i] = crx1[i]-crn4[i];
	    /*C5->P2*/
	    E12[i] = crx2[i]-crn4[i];
	    /*C5->P3*/
	    E13[i] = crx3[i]-crn4[i];
	}
	/*volume of C8P1P4P5*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C8P4P5P6*/
	V2 = fabs(Det3d(E1,E3,E4));
	/*volume of C8P5P6C4*/
	V3 = fabs(Det3d(E3,E4,E5));
	/*volume of C7C8P1P4*/
	V4 = fabs(Det3d(E6,E7,E8));
	/*volume of C7P1P4P2*/
	V5 = fabs(Det3d(E7,E8,E9));
	/*volume of C7P1P2C5*/
	V6 = fabs(Det3d(E8,E9,E10));
	/*volume of C5P1P2P3*/
	V7 = fabs(Det3d(E11,E12,E13));
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7);
	return volume;
}

LOCAL	double case10_complement(
		CELL_INFO_3D* cell3d,
		double crn1[],
		double crn2[],
		double crn3[],
		double crn6[])
{
	int i; 
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double E1[3],E2[3],E3[3],E4[3],E5[3],E6[3],E7[3],E8[3],E9[3],E10[3],
		E11[3],E12[3],E13[3];
	double V1,V2,V3,V4,V5,V6,V7,V8;
	double volume;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 
	crx5 = crx_in_kdir(cell3d,0,1); 
	crx6 = crx_in_jdir(cell3d,1,0); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->P1*/
	    E1[i] = crx1[i]-crn2[i];
	    /*C2->P4*/
	    E2[i] = crx4[i]-crn2[i];
	    /*C2->P5*/
	    E3[i] = crx5[i]-crn2[i];
	    /*C2->P6*/
	    E4[i] = crx6[i]-crn2[i];
	    /*C2->P2*/
	    E5[i] = crx2[i]-crn2[i];
	    /*C2->C6*/
	    E6[i] = crn6[i]-crn2[i];
	    /*C2->P3*/
	    E7[i] = crx3[i]-crn2[i];
	    /*C2->C3*/
	    E8[i] = crn3[i]-crn2[i];
	    /*C1->C2*/
	    E9[i] = crn2[i]-crn1[i];
	    /*C1->C3*/
	    E10[i] = crn3[i]-crn1[i];
	    /*C1->P1*/
	    E11[i] = crx1[i]-crn1[i];
	    /*C1->P3*/
	    E12[i] = crx3[i]-crn1[i];
	}
	/*volume of C2P1P4P5*/
	V1 = fabs(Det3d(E1,E2,E3));
	/*volume of C2P4P5P6*/
	V2 = fabs(Det3d(E2,E3,E4));
	/*volume of C2P1P2P4*/
	V3 = fabs(Det3d(E1,E2,E5));
	/*volume of C2P2P4C6*/
	V4 = fabs(Det3d(E2,E5,E6));
	/*volume of C2P1P2P3*/
	V5 = fabs(Det3d(E1,E5,E7));
	/*volume of C2P1P5C3*/
	V6 = fabs(Det3d(E1,E3,E8));
	/*volume of C1C2C3P1*/
	V7 = fabs(Det3d(E9,E10,E11));
	/*volume of C1C2P1P3*/
	V8 = fabs(Det3d(E9,E11,E12));
	volume = 1.0/6*(V1+V2+V3+V4+V5+V6+V7+V8);
	return volume;
}

LOCAL	double volume_edge1_case11(
	CELL_INFO_3D *cell3d,
	double crn3[],
	double crn4[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4;
	double D1[3], D2[3], D3[3], D4[3], D5[3];
	double V1, V2, V3;
	double volume;

        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_jdir(cell3d,0,0); 
        crx3 = crx_in_jdir(cell3d,1,0); 
	crx4 = crx_in_idir(cell3d,1,1);
	
	for (i = 0; i < 3; i++)
	{
	    /*C4->P1*/
	    D1[i] = crx1[i]-crn4[i];
	    /*C4->P3*/
	    D2[i] = crx3[i]-crn4[i];
	    /*C4->P4*/
	    D3[i] = crx4[i]-crn4[i];
	    /*C4->P2*/
	    D4[i] = crx2[i]-crn4[i];
	    /*C4->C3*/
	    D5[i] = crn3[i]-crn4[i];

	}
	/*volume of C4P1P3P4*/
	V1 = fabs(Det3d(D1,D2,D3));
	/*volume of C4P1P2P3*/
	V2 = fabs(Det3d(D1,D2,D4));
	/*volume of C4C3P1P2*/
	V3 = fabs(Det3d(D1,D4,D5));

	volume = ( V1+V2+V3)/6.0;
	return volume;
}	/* end volume_edge1_case11 */

LOCAL 	double volume_edge2_case11(
	CELL_INFO_3D *cell3d,
	double crn5[],
	double crn6[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4;
	double D1[3], D2[3], D3[3], D4[3], D5[3];
	double V1, V2, V3;
	double volume;

        crx1 = crx_in_jdir(cell3d,0,1); 
        crx2 = crx_in_idir(cell3d,0,1); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 

	for (i = 0; i < 3; i++)
	{
	    /*C6->P1*/
	    D1[i] = crx1[i]-crn6[i];
	    /*C6->P2*/
	    D2[i] = crx2[i]-crn6[i];
	    /*C6->P4*/
	    D3[i] = crx4[i]-crn6[i];
	    /*C6->P3*/
	    D4[i] = crx3[i]-crn6[i];
	    /*C6->C5*/
	    D5[i] = crn5[i]-crn6[i];
	}
	/* volume of C6P1P2P4*/
	V1 = fabs(Det3d(D1,D2,D3));
	/* volume of C6P1P2P3*/
	V2 = fabs(Det3d(D1,D2,D4));
	/* volume of C6C5P1P3*/
	V3 = fabs(Det3d(D1,D4,D5));
	volume = (V1+V2+V3)/6.0;
	return volume;
}	/* end volume_edge2_case11 */

LOCAL double volume_glider_case13(
	CELL_INFO_3D *cell3d,
	double crn5[],
	double crn6[],
	double crn7[])
{
	int i;
	const double *crx1, *crx2, *crx3, *crx4, *crx5;
	double E1[3], E2[3], E3[3], E4[3], E5[3], E6[3], E7[3], E8[3], 
	      E9[3], E10[3], E11[3];
	double V1, V2, V3, V4, V5;
	double volume;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,1);
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 
	crx5 = crx_in_idir(cell3d,0,1); 

	for (i = 0; i < 3; ++i)
	{
	    /*C5->P2*/
	    E1[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    E2[i] = crx3[i]-crn5[i];
	    /*C5->P4*/
	    E3[i] = crx4[i]-crn5[i];
	    /*C5->C6*/
	    E4[i] = crn6[i]-crn5[i];
	    /*P5->C6*/
	    E5[i] = crn6[i]-crx5[i];
	    /*P5->P4*/
	    E6[i] = crx4[i]-crx5[i];
	    /*P5->P3*/
	    E7[i] = crx3[i]-crx5[i];
	    /*C5->C7*/
	    E8[i] = crn7[i]-crn5[i];
	    /*P1->P2*/
	    E9[i] = crx2[i]-crx1[i];
	    /*P1->P3*/
	    E10[i] = crx3[i]-crx1[i];
	    /*P1->C7*/
	    E11[i] = crn7[i]-crx1[i];
	}
	/* volume of C5P2P3P4 */
	V1 = fabs(Det3d(E1,E2,E3));
	/* volume of C5C6P3P4 */
	V2 = fabs(Det3d(E2,E3,E4));
	/* volume of P5P3P4C6 */
	V3 = fabs(Det3d(E5,E6,E7));
	/* volume of C5P2P3C7 */
	V4 = fabs(Det3d(E8,E1,E2));
	/* volume of P1P2P3C7 */
	V5 = fabs(Det3d(E9,E10,E11));
	
	volume = 1.0/6*(V1+V2+V3+V4+V5);
	return volume;
	
}	/* end volume_glider_case13 */

LOCAL 	double volume_corner1_case14(
	CELL_INFO_3D *cell3d,
	double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;
	crx1 = crx_in_idir(cell3d,1,0);
        crx2 = crx_in_jdir(cell3d,0,0);
        crx3 = crx_in_kdir(cell3d,0,1);

	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c3  */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}	/* end volume_corner1_case14 */



LOCAL 	double volume_corner4_case14(
	CELL_INFO_3D *cell3d,
	double crn[])
{
	int i;
	const double *crx1, *crx2, *crx3; 
	double D1[3], D2[3], D3[3];
	double volume;

	crx1 = crx_in_jdir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,0,0);
        crx3 = crx_in_idir(cell3d,0,1);
	
	for (i = 0; i < 3; i++)
	{
	    D1[i] = crx1[i] - crx2[i];
	    D2[i] = crx2[i] - crx3[i];
	    D3[i] = crx3[i] - crn[i];	/* crx3 - c2  */
	}
	
	volume = fabs(Det3d(D1,D2,D3))/6.0;
	return volume;
}	/* end volume_corner4_case14 */

LOCAL 	const double* crx_in_idir(
	const CELL_INFO_3D *cell3d,
	int j,
	int k)
{
	int ip[3],ipn[3];
	const int (*ix)[2][2] = cell3d->ix;
	const int (*iy)[2][2] = cell3d->iy;
	const int (*iz)[2][2] = cell3d->iz;
	const double *crx;

	ip[0] = ix[0][j][k];
	ip[1] = iy[0][j][k];
	ip[2] = iz[0][j][k];
	ipn[0] = ix[1][j][k];
	ipn[1] = iy[1][j][k];
	ipn[2] = iz[1][j][k];
	crx = crx_in_between(ip,ipn,cell3d);
	return crx;
}

LOCAL	const double* crx_in_jdir(
	const CELL_INFO_3D* cell3d,
	int k,
	int i)
{
	int ip[3],ipn[3];
	const int (*ix)[2][2] = cell3d->ix;
	const int (*iy)[2][2] = cell3d->iy;
	const int (*iz)[2][2] = cell3d->iz;
	const double *crx;

	ip[0] = ix[i][0][k];
	ip[1] = iy[i][0][k];
	ip[2] = iz[i][0][k];
	ipn[0] = ix[i][1][k];
	ipn[1] = iy[i][1][k];
	ipn[2] = iz[i][1][k];
	crx = crx_in_between(ip,ipn,cell3d);
	return crx;
}

LOCAL	const double* crx_in_kdir(
	const CELL_INFO_3D* cell3d,
	int i,
	int j)
{
	int ip[3],ipn[3];
	const int (*ix)[2][2] = cell3d->ix;
	const int (*iy)[2][2] = cell3d->iy;
	const int (*iz)[2][2] = cell3d->iz;
	const double *crx;

	ip[0] = ix[i][j][0];
	ip[1] = iy[i][j][0];
	ip[2] = iz[i][j][0];
	ipn[0] = ix[i][j][1];
	ipn[1] = iy[i][j][1];
	ipn[2] = iz[i][j][1];
	crx = crx_in_between(ip,ipn,cell3d);
	return crx;
}

LOCAL 	const double* crx_in_between(
	const int *ip,
	const int *ipn,
	const CELL_INFO_3D* cell3d)
{
	if (ip[0] != ipn[0])
	{
	    if (cell3d->crx_coords[0][ip[1]][ip[2]][0] == HUGE)
	    {
	        screen("ERROR in crxing_in_between(), "
		       "no crossing in i-dir between:\n"
		       "(%d %d %d) and (%d %d %d)\n",
		       ip[0],ip[1],ip[2],ipn[0],ipn[1],ipn[2]);
		clean_up(ERROR);
	    }
	    return cell3d->crx_coords[0][ip[1]][ip[2]];
	}
	else if(ip[1] != ipn[1])
	{
	    if (cell3d->crx_coords[1][ip[2]][ip[0]][0] == HUGE)
	    {
	        screen("ERROR in crxing_in_between(), "
		       "no crossing in j-dir between:\n"
		       "(%d %d %d) and (%d %d %d)\n",
		       ip[0],ip[1],ip[2],ipn[0],ipn[1],ipn[2]);
		clean_up(ERROR);
	    }
	    return cell3d->crx_coords[1][ip[2]][ip[0]];
	}
	else if(ip[2] != ipn[2])
	{
	    if (cell3d->crx_coords[2][ip[0]][ip[1]][0] == HUGE)
	    {
	        screen("ERROR in crxing_in_between(), "
		       "no crossing in k-dir between:\n"
		       "(%d %d %d) and (%d %d %d)\n",
		       ip[0],ip[1],ip[2],ipn[0],ipn[1],ipn[2]);
		clean_up(ERROR);
	    }
	    return cell3d->crx_coords[2][ip[0]][ip[1]];
	}
	screen("ERRO in crxing_in_between(), "
	       "Inconsistent values of ip and ipn arrarys\n");
	clean_up(ERROR);
	return NULL;
}

LOCAL	void set_prime_components(
	COMPONENT ****pcomp)
{
	/* Case 1: */

	pcomp[0][0][0][0] = 0;
	pcomp[0][0][0][1] = 0;
	pcomp[0][0][1][0] = 0;
	pcomp[0][0][1][1] = 0;
	pcomp[0][1][0][0] = 1;
	pcomp[0][1][0][1] = 1;
	pcomp[0][1][1][0] = 1;
	pcomp[0][1][1][1] = 1;

	/* Case 2: */

	pcomp[1][0][0][0] = 0;
	pcomp[1][0][0][1] = 0;
	pcomp[1][0][1][0] = 0;
	pcomp[1][0][1][1] = 0;
	pcomp[1][1][0][0] = 0;
	pcomp[1][1][0][1] = 0;
	pcomp[1][1][1][0] = 1;
	pcomp[1][1][1][1] = 1;

	/* Case 3: */

	pcomp[2][0][0][0] = 0;
	pcomp[2][0][0][1] = 0;
	pcomp[2][0][1][0] = 0;
	pcomp[2][0][1][1] = 0;
	pcomp[2][1][0][0] = 0;
	pcomp[2][1][0][1] = 0;
	pcomp[2][1][1][0] = 0;
	pcomp[2][1][1][1] = 1;

	/* Case 4: */

	pcomp[3][0][0][0] = 0;
	pcomp[3][0][0][1] = 0;
	pcomp[3][0][1][0] = 0;
	pcomp[3][0][1][1] = 0;
	pcomp[3][1][0][0] = 0;
	pcomp[3][1][0][1] = 1;
	pcomp[3][1][1][0] = 1;
	pcomp[3][1][1][1] = 1;

	/* Case 5: */

	pcomp[4][0][0][0] = 0;
	pcomp[4][0][0][1] = 0;
	pcomp[4][0][1][0] = 0;
	pcomp[4][0][1][1] = 1;
	pcomp[4][1][0][0] = 0;
	pcomp[4][1][0][1] = 1;
	pcomp[4][1][1][0] = 1;
	pcomp[4][1][1][1] = 1;

	/* Case 6: */

	pcomp[5][0][0][0] = 0;
	pcomp[5][0][0][1] = 0;
	pcomp[5][0][1][0] = 0;
	pcomp[5][0][1][1] = 0;
	pcomp[5][1][0][0] = 0;
	pcomp[5][1][0][1] = 1;
	pcomp[5][1][1][0] = 1;
	pcomp[5][1][1][1] = 0;

	/* Case 7: */

	pcomp[6][0][0][0] = 0;
	pcomp[6][0][0][1] = 0;
	pcomp[6][0][1][0] = 0;
	pcomp[6][0][1][1] = 1;
	pcomp[6][1][0][0] = 1;
	pcomp[6][1][0][1] = 1;
	pcomp[6][1][1][0] = 0;
	pcomp[6][1][1][1] = 1;

	/* Case 8: */

	pcomp[7][0][0][0] = 0;
	pcomp[7][0][0][1] = 0;
	pcomp[7][0][1][0] = 0;
	pcomp[7][0][1][1] = 1;
	pcomp[7][1][0][0] = 1;
	pcomp[7][1][0][1] = 0;
	pcomp[7][1][1][0] = 0;
	pcomp[7][1][1][1] = 1;

	/* Case 9: */

	pcomp[8][0][0][0] = 0;
	pcomp[8][0][0][1] = 0;
	pcomp[8][0][1][0] = 0;
	pcomp[8][0][1][1] = 1;
	pcomp[8][1][0][0] = 1;
	pcomp[8][1][0][1] = 0;
	pcomp[8][1][1][0] = 0;
	pcomp[8][1][1][1] = 0;

	/* Case 10: */

	pcomp[9][0][0][0] = 0;
	pcomp[9][0][0][1] = 0;
	pcomp[9][0][1][0] = 0;
	pcomp[9][0][1][1] = 1;
	pcomp[9][1][0][0] = 1;
	pcomp[9][1][0][1] = 0;
	pcomp[9][1][1][0] = 1;
	pcomp[9][1][1][1] = 1;

	/* Case 11: */

	pcomp[10][0][0][0] = 0;
	pcomp[10][0][0][1] = 0;
	pcomp[10][0][1][0] = 1;
	pcomp[10][0][1][1] = 1;
	pcomp[10][1][0][0] = 1;
	pcomp[10][1][0][1] = 1;
	pcomp[10][1][1][0] = 0;
	pcomp[10][1][1][1] = 0;

	/* Case 12: */

	pcomp[11][0][0][0] = 0;
	pcomp[11][0][0][1] = 0;
	pcomp[11][0][1][0] = 0;
	pcomp[11][0][1][1] = 1;
	pcomp[11][1][0][0] = 0;
	pcomp[11][1][0][1] = 1;
	pcomp[11][1][1][0] = 1;
	pcomp[11][1][1][1] = 0;

	/* Case 13: */

	pcomp[12][0][0][0] = 0;
	pcomp[12][0][0][1] = 0;
	pcomp[12][0][1][0] = 0;
	pcomp[12][0][1][1] = 1;
	pcomp[12][1][0][0] = 1;
	pcomp[12][1][0][1] = 1;
	pcomp[12][1][1][0] = 1;
	pcomp[12][1][1][1] = 0;

	/* Case 14: */

	pcomp[13][0][0][0] = 0;
	pcomp[13][0][0][1] = 1;
	pcomp[13][0][1][0] = 1;
	pcomp[13][0][1][1] = 0;
	pcomp[13][1][0][0] = 1;
	pcomp[13][1][0][1] = 0;
	pcomp[13][1][1][0] = 0;
	pcomp[13][1][1][1] = 1;
}	/* end set_prime_components */

LOCAL void copy_cell3d(
	const CELL_INFO_3D* cell1,
	CELL_INFO_3D* cell2)
{
	int i,j,k,l;
	cell2->is_corner = cell1->is_corner;
	cell2->is_side = cell1->is_side;
	cell2->is_face = cell1->is_face;
	cell2->num_comps = cell1->num_comps;
	cell2->orphan = cell1->orphan;
	cell2->full_cell_vol = cell1->full_cell_vol;

	for (i = 0; i < 2; ++i)
	{
	    cell2->comps[i] = cell1->comps[i];
	    cell2->nv[i] = cell1->nv[i];
	}
	
	for (i = 0; i < 2; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	{
	    cell2->ix[i][j][k] = cell1->ix[i][j][k];
	    cell2->iy[i][j][k] = cell1->iy[i][j][k];
	    cell2->iz[i][j][k] = cell1->iz[i][j][k];
	    for (l = 0; l < 3; ++l)
	    {
	    	cell2->icrds[i][j][k][l] = cell1->icrds[i][j][k][l];
	    	cell2->crds[i][j][k][l] = cell1->crds[i][j][k][l];
	    }
	    if (cell1->comp[i][j][k] == cell1->comps[0])
	    	cell2->comp[i][j][k] = 1;
	    else
	    	cell2->comp[i][j][k] = 0;
	}

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	for (l = 0; l < 3; ++l)
	{
	    cell2->crx_coords[i][j][k][l] = cell1->crx_coords[i][j][k][l];
	}

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    cell2->nb_flag[i][j][k] = cell1->nb_flag[i][j][k];
	    cell2->nb_frac[i][j][k] = cell1->nb_frac[i][j][k];
	    cell2->area[i][j][k] = cell1->area[i][j][k];
	}

	if (cell1->soln_comp == cell1->comps[0])
	    cell2->soln_comp = 1;
	else
	    cell2->soln_comp = 0;
	
	if (cell1->cell_comp == cell1->comps[0])
	    cell2->cell_comp = 1;
	else
	    cell2->cell_comp = 0;
	
	cell2->is_soln_solute = cell1->is_soln_solute;
	cell2->is_center_solute = cell1->is_center_solute;	
	cell2->unstable_cell = cell1->unstable_cell;
	cell2->solute_volume = cell1->solute_volume;
	return; 
}

LOCAL   int compare_comp(
	COMPONENT comp[][2][2],
	COMPONENT ****p_comp,
	int ii)
{
	int i,j,k;
        
	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
        for (k = 0; k < 2; k++) 
	if (comp[i][j][k] != p_comp[ii][i][j][k])
	{
            return NO; 
	}
	return YES; 
}       /* end compare_comp */

LOCAL	void rot24(
        CELL_INFO_3D *cell3d,
	int n)
{       
        if (n < 16)
        {
            x_rotation(cell3d);
            if ((n+1)%4 == 0)
                y_rotation(cell3d);
        }
        else if (n == 16)
        {   
            z_rotation(cell3d); 
            x_rotation(cell3d); 
        }   
        else if (n == 20)
        {   
	    z_rotation(cell3d);
	    z_rotation(cell3d);
	    x_rotation(cell3d);
	}
	else if (n < 24)
	{
	    x_rotation(cell3d);
	}
}       /* end rot24 */

LOCAL   void x_rotation(CELL_INFO_3D *cell3d)
{
        int i, temp;
	COMPONENT comp_tmp;
	boolean temp_flag;
	double temp_frac;
	int (*ix)[2][2] = cell3d->ix;
	int (*iy)[2][2] = cell3d->iy;
	int (*iz)[2][2] = cell3d->iz;
	COMPONENT (*comp)[2][2] = cell3d->comp;
	double (*nb_frac)[3][3] = cell3d->nb_frac;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;

	for (i = 0; i < 2; i++)
	{ 
	    comp_tmp = comp[i][0][0];
	    comp[i][0][0] = comp[i][1][0];
	    comp[i][1][0] = comp[i][1][1];
	    comp[i][1][1] = comp[i][0][1];
	    comp[i][0][1] = comp_tmp;
		
	    temp = ix[i][0][0];
	    ix[i][0][0] = ix[i][1][0];
	    ix[i][1][0] = ix[i][1][1];
	    ix[i][1][1] = ix[i][0][1]; 
	    ix[i][0][1] = temp;
		
		
	    temp = iy[i][0][0];
	    iy[i][0][0] = iy[i][1][0];
	    iy[i][1][0] = iy[i][1][1];
	    iy[i][1][1] = iy[i][0][1];    
	    iy[i][0][1] = temp;
		
	    temp = iz[i][0][0];
	    iz[i][0][0] = iz[i][1][0];
	    iz[i][1][0] = iz[i][1][1];
	    iz[i][1][1] = iz[i][0][1];    
	    iz[i][0][1] = temp;
	}

	for (i = 0; i < 3; ++i)
	{
	    temp_flag = nb_flag[i][0][0];
	    nb_flag[i][0][0] = nb_flag[i][2][0];
	    nb_flag[i][2][0] = nb_flag[i][2][2];
	    nb_flag[i][2][2] = nb_flag[i][0][2];
	    nb_flag[i][0][2] = temp_flag;

	    temp_flag = nb_flag[i][1][0];
	    nb_flag[i][1][0] = nb_flag[i][2][1];
	    nb_flag[i][2][1] = nb_flag[i][1][2];
	    nb_flag[i][1][2] = nb_flag[i][0][1];
	    nb_flag[i][0][1] = temp_flag;
	    
	    temp_frac = nb_frac[i][0][0];
	    nb_frac[i][0][0] = nb_frac[i][2][0];
	    nb_frac[i][2][0] = nb_frac[i][2][2];
	    nb_frac[i][2][2] = nb_frac[i][0][2];
	    nb_frac[i][0][2] = temp_frac;

	    temp_frac = nb_frac[i][1][0];
	    nb_frac[i][1][0] = nb_frac[i][2][1];
	    nb_frac[i][2][1] = nb_frac[i][1][2];
	    nb_frac[i][1][2] = nb_frac[i][0][1];
	    nb_frac[i][0][1] = temp_frac;
	    
	}
}       /* end x_rotation */

LOCAL   void y_rotation(CELL_INFO_3D *cell3d)
{
        int j, temp;
	COMPONENT comp_tmp;
	boolean temp_flag;
	double temp_frac;
	int (*ix)[2][2] = cell3d->ix;
	int (*iy)[2][2] = cell3d->iy;
	int (*iz)[2][2] = cell3d->iz;
	COMPONENT (*comp)[2][2] = cell3d->comp;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;
	double (*nb_frac)[3][3] = cell3d->nb_frac;

	for (j = 0; j < 2; j++)
	{
            comp_tmp = comp[0][j][0];
	    comp[0][j][0] = comp[0][j][1];
	    comp[0][j][1] = comp[1][j][1];
	    comp[1][j][1] = comp[1][j][0];
	    comp[1][j][0] = comp_tmp;

	    temp = ix[0][j][0];
	    ix[0][j][0] = ix[0][j][1];
	    ix[0][j][1] = ix[1][j][1];
	    ix[1][j][1] = ix[1][j][0];
	    ix[1][j][0] = temp;

	    temp = iy[0][j][0];
	    iy[0][j][0] = iy[0][j][1];
	    iy[0][j][1] = iy[1][j][1];
	    iy[1][j][1] = iy[1][j][0];
	    iy[1][j][0] = temp;

	    temp = iz[0][j][0];
	    iz[0][j][0] = iz[0][j][1];
	    iz[0][j][1] = iz[1][j][1];
	    iz[1][j][1] = iz[1][j][0];
	    iz[1][j][0] = temp;
	}

	for (j = 0; j < 3; ++j)
	{
	    temp_flag = nb_flag[0][j][0];
	    nb_flag[0][j][0] = nb_flag[0][j][2];
	    nb_flag[0][j][2] = nb_flag[2][j][2];
	    nb_flag[2][j][2] = nb_flag[2][j][0];
	    nb_flag[2][j][0] = temp_flag;

	    temp_flag = nb_flag[0][j][1];
	    nb_flag[0][j][1] = nb_flag[1][j][2];
	    nb_flag[1][j][2] = nb_flag[2][j][1];
	    nb_flag[2][j][1] = nb_flag[1][j][0];
	    nb_flag[1][j][0] = temp_flag;
	    
	    temp_frac = nb_frac[0][j][0];
	    nb_frac[0][j][0] = nb_frac[0][j][2];
	    nb_frac[0][j][2] = nb_frac[2][j][2];
	    nb_frac[2][j][2] = nb_frac[2][j][0];
	    nb_frac[2][j][0] = temp_frac;

	    temp_frac = nb_frac[0][j][1];
	    nb_frac[0][j][1] = nb_frac[1][j][2];
	    nb_frac[1][j][2] = nb_frac[2][j][1];
	    nb_frac[2][j][1] = nb_frac[1][j][0];
	    nb_frac[1][j][0] = temp_frac;

	}
}       /* y_rotation */

LOCAL   void z_rotation(CELL_INFO_3D *cell3d)
{
        int k, temp;
	boolean temp_flag;
	double temp_frac;
	COMPONENT comp_tmp;
	int (*ix)[2][2] = cell3d->ix;
	int (*iy)[2][2] = cell3d->iy;
	int (*iz)[2][2] = cell3d->iz;
	COMPONENT (*comp)[2][2] = cell3d->comp;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;
	double (*nb_frac)[3][3] = cell3d->nb_frac;

	for (k = 0; k < 2; k++)
	{
	    comp_tmp = comp[0][0][k];
	    comp[0][0][k] = comp[1][0][k];
	    comp[1][0][k] = comp[1][1][k];
	    comp[1][1][k] = comp[0][1][k];
	    comp[0][1][k] = comp_tmp;

	    temp = ix[0][0][k];
	    ix[0][0][k] = ix[1][0][k];
	    ix[1][0][k] = ix[1][1][k];
	    ix[1][1][k] = ix[0][1][k];
	    ix[0][1][k] = temp;

	    temp = iy[0][0][k];
	    iy[0][0][k] = iy[1][0][k];
	    iy[1][0][k] = iy[1][1][k];
	    iy[1][1][k] = iy[0][1][k];
	    iy[0][1][k] = temp;
		
	    temp = iz[0][0][k];
	    iz[0][0][k] = iz[1][0][k];
	    iz[1][0][k] = iz[1][1][k];
	    iz[1][1][k] = iz[0][1][k];
	    iz[0][1][k] = temp;
	}

	for (k = 0; k < 3; ++k)
	{
	    temp_flag = nb_flag[0][0][k];
	    nb_flag[0][0][k] = nb_flag[2][0][k];
	    nb_flag[2][0][k] = nb_flag[2][2][k];
	    nb_flag[2][2][k] = nb_flag[0][2][k];
	    nb_flag[0][2][k] = temp_flag;

	    temp_flag = nb_flag[1][0][k];
	    nb_flag[1][0][k] = nb_flag[2][1][k];
	    nb_flag[2][1][k] = nb_flag[1][2][k];
	    nb_flag[1][2][k] = nb_flag[0][1][k];
	    nb_flag[0][1][k] = temp_flag;

	    temp_frac= nb_frac[0][0][k];
	    nb_frac[0][0][k] = nb_frac[2][0][k];
	    nb_frac[2][0][k] = nb_frac[2][2][k];
	    nb_frac[2][2][k] = nb_frac[0][2][k];
	    nb_frac[0][2][k] = temp_frac;

	    temp_frac = nb_frac[1][0][k];
	    nb_frac[1][0][k] = nb_frac[2][1][k];
	    nb_frac[2][1][k] = nb_frac[1][2][k];
	    nb_frac[1][2][k] = nb_frac[0][1][k];
	    nb_frac[0][1][k] = temp_frac;
	}
}       /* end z_rotation */

LOCAL	void reverse24(
	CELL_INFO_3D *cell3d,
	int n)
{
	if (n-1 < 16)
	{
	    if (n%4 == 0)
	    	y_reverse(cell3d);
	    x_reverse(cell3d);
	}
	else if (n-1 == 16)
	{
	    x_reverse(cell3d);
	    z_reverse(cell3d);
	}
	else if (n-1 == 20)
	{
	    x_reverse(cell3d);
	    z_reverse(cell3d);
	    z_reverse(cell3d);
	}
	else if (n-1 < 24)
	{
	    x_reverse(cell3d);
	}
}

LOCAL	void x_reverse(CELL_INFO_3D* cell3d)
{
	int i;
	boolean temp_flag;
	double  temp_frac;
	double (*nb_frac)[3][3] = cell3d->nb_frac;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;
	
	for (i = 0; i < 3; ++i)
	{
	    temp_frac = nb_frac[i][0][0] ;
	    nb_frac[i][0][0] = nb_frac[i][0][2];
	    nb_frac[i][0][2] = nb_frac[i][2][2];
	    nb_frac[i][2][2] = nb_frac[i][2][0];
	    nb_frac[i][2][0] = temp_frac;

	    temp_frac = nb_frac[i][1][0];
	    nb_frac[i][1][0] = nb_frac[i][0][1];
	    nb_frac[i][0][1] = nb_frac[i][1][2];
	    nb_frac[i][1][2] = nb_frac[i][2][1];
	    nb_frac[i][2][1] = temp_frac;

	    temp_flag = nb_flag[i][0][0] ;
	    nb_flag[i][0][0] = nb_flag[i][0][2];
	    nb_flag[i][0][2] = nb_flag[i][2][2];
	    nb_flag[i][2][2] = nb_flag[i][2][0];
	    nb_flag[i][2][0] = temp_flag;

	    temp_flag = nb_flag[i][1][0];
	    nb_flag[i][1][0] = nb_flag[i][0][1];
	    nb_flag[i][0][1] = nb_flag[i][1][2];
	    nb_flag[i][1][2] = nb_flag[i][2][1];
	    nb_flag[i][2][1] = temp_flag;
	}
}

LOCAL	void y_reverse(CELL_INFO_3D* cell3d)
{
	int j;
	boolean temp_flag;
	double temp_frac;
	double (*nb_frac)[3][3] = cell3d->nb_frac;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;

	for (j = 0; j < 3; ++j)
	{
	     temp_frac = nb_frac[0][j][0];
	     nb_frac[0][j][0] = nb_frac[2][j][0];
	     nb_frac[2][j][0] = nb_frac[2][j][2];
	     nb_frac[2][j][2] = nb_frac[0][j][2];
	     nb_frac[0][j][2] = temp_frac;

	     temp_frac = nb_frac[0][j][1];
	     nb_frac[0][j][1] = nb_frac[1][j][0];
	     nb_frac[1][j][0] = nb_frac[2][j][1];
	     nb_frac[2][j][1] = nb_frac[1][j][2];
	     nb_frac[1][j][2] = temp_frac;
	
	     temp_flag = nb_flag[0][j][0];
	     nb_flag[0][j][0] = nb_flag[2][j][0];
	     nb_flag[2][j][0] = nb_flag[2][j][2];
	     nb_flag[2][j][2] = nb_flag[0][j][2];
	     nb_flag[0][j][2] = temp_flag;

	     temp_flag = nb_flag[0][j][1];
	     nb_flag[0][j][1] = nb_flag[1][j][0];
	     nb_flag[1][j][0] = nb_flag[2][j][1];
	     nb_flag[2][j][1] = nb_flag[1][j][2];
	     nb_flag[1][j][2] = temp_flag;
	}
}

LOCAL	void z_reverse(CELL_INFO_3D* cell3d)
{
	int k;
	boolean temp_flag;
	double temp_frac;
	double (*nb_frac)[3][3] = cell3d->nb_frac;
	boolean (*nb_flag)[3][3] = cell3d->nb_flag;

	for (k = 0; k < 3; ++k)
	{
	    temp_frac = nb_frac[0][0][k];
	    nb_frac[0][0][k] = nb_frac[0][2][k];
	    nb_frac[0][2][k] = nb_frac[2][2][k];
	    nb_frac[2][2][k] = nb_frac[2][0][k];
	    nb_frac[2][0][k] = temp_frac;

	    temp_frac = nb_frac[1][0][k];
	    nb_frac[1][0][k] = nb_frac[0][1][k];
	    nb_frac[0][1][k] = nb_frac[1][2][k];
	    nb_frac[1][2][k] = nb_frac[2][1][k];
	    nb_frac[2][1][k] = temp_frac;
	
	    temp_flag = nb_flag[0][0][k];
	    nb_flag[0][0][k] = nb_flag[0][2][k];
	    nb_flag[0][2][k] = nb_flag[2][2][k];
	    nb_flag[2][2][k] = nb_flag[2][0][k];
	    nb_flag[2][0][k] = temp_flag;

	    temp_flag = nb_flag[1][0][k];
	    nb_flag[1][0][k] = nb_flag[0][1][k];
	    nb_flag[0][1][k] = nb_flag[1][2][k];
	    nb_flag[1][2][k] = nb_flag[2][1][k];
	    nb_flag[2][1][k] = temp_flag;
	}
}

LOCAL	void copy_back_cell3d(
	const CELL_INFO_3D* cell1,
	CELL_INFO_3D* cell2)
{
	int i,j,k;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (cell1->nb_flag[i][j][k] != cell2->nb_flag[i][j][k])
	    {
	    	printf("Two cells with different nb_flag\n");
		clean_up(ERROR);
	    }
	}

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    cell2->nb_frac[i][j][k] = cell1->nb_frac[i][j][k];
	}
	cell2->unstable_cell = cell1->unstable_cell;
	cell2->solute_volume = cell1->solute_volume;
	cell2->orphan = cell1->orphan;
	return;
}

LOCAL	void cell3d_area_case01_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double d1[3],d2[3],d3[3],d4[3],d5[3],d6[3],d7[3],d8[3];
	double L1,L2,L3,L4,L5,L6,L7,L8;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = 0.0;
 	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_idir(cell3d,0,1);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C4->P2*/
	    d1[i] = crx2[i]-crn4[i]; 
	    /*C4->C2*/
	    d2[i] = crn2[i]-crn4[i];
	    /*C2->P4*/
	    d3[i] = crx4[i]-crn2[i];
	    /*C3->P1*/
	    d4[i] = crx1[i]-crn3[i];
	    /*C3->C1*/
	    d5[i] = crn1[i]-crn3[i];
	    /*C1->P3*/
	    d6[i] = crx3[i]-crn1[i];
	    /*C4->C3*/
	    d7[i] = crn3[i]-crn4[i];
	    /*C1->C2*/
	    d8[i] = crn2[i]-crn1[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(d1[i]);
	    L2 += fabs(d2[i]);
	    L3 += fabs(d3[i]);
	    L4 += fabs(d4[i]);
	    L5 += fabs(d5[i]);
	    L6 += fabs(d6[i]);
	    L7 += fabs(d7[i]);
	    L8 += fabs(d8[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L1+L3)*L2/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L7/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L3+L6)*L8/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L4+L6)*L5/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L7 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j; 
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else 
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case01_comp1(
	CELL_INFO_3D *cell3d,
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3];
	double L1,L2,L3,L4,L5,L6,L7,L8;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = 0.0;
 	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_idir(cell3d,0,1);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P2*/
	    D1[i] = crx2[i]-crn8[i]; 
	    /*C6->C8*/
	    D2[i] = crn8[i]-crn6[i]; 
	    /*C6->P4*/
	    D3[i] = crx4[i]-crn6[i];
	    /*C7->P1*/
	    D4[i] = crx1[i]-crn7[i];
	    /*C7->C5*/
	    D5[i] = crn5[i]-crn7[i];
	    /*C5->P3*/
	    D6[i] = crx3[i]-crn5[i];
	    /*C8->C7*/
	    D7[i] = crn7[i]-crn8[i];
	    /*C6->C5*/
	    D8[i] = crn5[i]-crn6[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	}
	area[1][1][2] = nb_flag[1][1][2] == YES? (L1+L3)*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L3+L6)*L8/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L4+L6)*L5/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L7/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[2][1][1] = nb_flag[2][1][1] == YES? L2*L8 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? 0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? 0 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;
	

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
	    	max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case02_comp1(
	CELL_INFO_3D *cell3d, 
	double crn7[],
	double crn8[], 
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double D1[3],D2[3],D3[3],D4[3],D5[3];
	double L1,L2,L3,L4,L5;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = 0.0;

	crx1 = crx_in_idir(cell3d,1,1);  
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,1,0);
	crx4 = crx_in_jdir(cell3d,0,1);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P1*/
	    D1[i] = crx1[i]-crn8[i];
	    /*C8->P2*/
	    D2[i] = crx2[i]-crn8[i];
	    /*C7->P3*/
	    D3[i] = crx3[i]-crn7[i];
	    /*C7->P4*/
	    D4[i] = crx4[i]-crn7[i];
	    /*C7->C8*/
	    D5[i] = crn8[i]-crn7[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L2/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L3)*L5/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L3*L4/2.0: -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	
	area[2][1][1] = nb_flag[2][1][1] == YES? (L2+L4)*L5/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
	    	max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case02_comp0(
	CELL_INFO_3D *cell3d, 
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,
		L9,L10,L11,L12,L13,L14;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8
	 = L9 = L10 = L11 = L12 = L13 = L14 = 0.0;

	crx1 = crx_in_idir(cell3d,1,1);  
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_idir(cell3d,1,0);
	crx4 = crx_in_jdir(cell3d,0,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C4->P1*/
	    D1[i] = crx1[i]-crn4[i];
	    /*C4->C2*/
	    D2[i] = crn2[i]-crn4[i];
	    /*C2->C6*/
	    D3[i] = crn6[i]-crn2[i];
	    /*C3->P3*/
	    D4[i] = crx3[i]-crn3[i];
	    /*C3->C1*/
	    D5[i] = crn1[i]-crn3[i];
	    /*C1->C5*/
	    D6[i] = crn5[i]-crn1[i];
	    /*C3->C4*/
	    D7[i] = crn4[i]-crn3[i];
	    /*C1->C2*/
	    D8[i] = crn2[i]-crn1[i];
	    /*C5->C6*/
	    D9[i] = crn6[i]-crn5[i];
	    /*C8->P1*/
	    D10[i] = crx1[i]-crn8[i];
	    /*C8->P2*/
	    D11[i] = crx2[i]-crn8[i];
	    /*C7->P3*/
	    D12[i] = crx3[i]-crn7[i];
	    /*C7->P4*/
	    D13[i] = crx4[i]-crn7[i];
	    /*C7->C8*/
	    D14[i] = crn8[i]-crn7[i];

	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L2*L3-L10*L11/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L3*L8 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L5*L6-L12*L13/2.0: -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L7/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? 0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L5*L8-(L11+L13)*L14/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L7*L5: -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? 0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
	    	max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;

}

LOCAL void cell3d_area_case03_comp1(
	CELL_INFO_3D *cell3d,
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double D1[3],D2[3],D3[3];
	double L1,L2,L3;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = 0.0;

	crx1 = crx_in_idir(cell3d,1,1); 
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_kdir(cell3d,1,1); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C8->P1*/
	    D1[i] = crx1[i]-crn8[i];
	    /*C8->P2*/
	    D2[i] = crx2[i]-crn8[i];
	    /*C8->P3*/
	    D3[i] = crx3[i]-crn8[i];
	}
	
	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	}

	area[1][2][1] = nb_flag[1][2][1] == YES? L1*L3/2.0 : -1.0;
	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L2/2.0 : -1.0;
	area[2][1][1] = nb_flag[2][1][1] == YES? L2*L3/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case03_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{	
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],
		D7[3],D8[3],D9[3],D10[3],D11[3],D12[3];
	double d1[3],d2[3],d3[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,
		L9,L10,L11,L12;
	double l1,l2,l3;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 
	 = L8 = L9 = L10 = L11 = L12 = 0.0;
	l1 = l2 = l3 = 0.0;
	
	crx1 = crx_in_idir(cell3d,1,1); 
        crx2 = crx_in_jdir(cell3d,1,1); 
        crx3 = crx_in_kdir(cell3d,1,1); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C4->C8*/
	    D1[i] = crn8[i]-crn4[i]; 
	    /*C4->C2*/
	    D2[i] = crn2[i]-crn4[i];
	    /*C2->C6*/
	    D3[i] = crn6[i]-crn2[i];
	    /*C6->C8*/
	    D4[i] = crn8[i]-crn6[i];
	    /*C3->C7*/
	    D5[i] = crn7[i]-crn3[i];
	    /*C3->C1*/
	    D6[i] = crn1[i]-crn3[i];
	    /*C1->C5*/
	    D7[i] = crn5[i]-crn1[i];
	    /*C7->C5*/
	    D8[i] = crn5[i]-crn7[i];
	    /*C3->C4*/
	    D9[i] = crn4[i]-crn3[i];
	    /*C1->C2*/
	    D10[i] = crn2[i]-crn1[i];
	    /*C5->C6*/
	    D11[i] = crn6[i]-crn5[i];
	    /*C7->C8*/
	    D12[i] = crn8[i]-crn7[i];
	    /*P1->C8*/
	    d1[i] = crn8[i]-crx1[i];
	    /*P2->C8*/
	    d2[i] = crn8[i]-crx2[i];
	    /*P3->C8*/
	    d3[i] = crn8[i]-crx3[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L4-(l1*l2/2.0) : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L1*L12-(l1*l3/2.0) : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L3*L11 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L6*L7 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? 0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? 0 : -1.0;
	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L9 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	area[2][1][1] = nb_flag[2][1][1] == YES? L4*L12-l2*l3/2.0 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? 0 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;
	area[2][2][0] = (nb_flag[2][2][0] == YES) ? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case04_comp1(
	CELL_INFO_3D *cell3d,
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],
		D7[3],D8[3],D9[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 
	 = L8 = L9  =  0.0;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
       	crx3 = crx_in_jdir(cell3d,0,1); 
	crx4 = crx_in_kdir(cell3d,1,0); 
	crx5 = crx_in_idir(cell3d,0,1); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*P2->C8*/
	    D1[i] = crn8[i]-crx2[i];
	    /*C8->C6*/
	    D2[i] = crn6[i]-crn8[i];
	    /*C6->P5*/
	    D3[i] = crx5[i]-crn6[i];
	    /*P1->C7*/
	    D4[i] = crn7[i]-crx1[i];
	    /*C7->P3*/
	    D5[i] = crx3[i]-crn7[i];
	    /*C7->C8*/
	    D6[i] = crn8[i]-crn7[i];
	    /*P4->C6*/
	    D7[i] = crn6[i]-crx4[i];
	    /*C5->P4*/
	    D8[i] = crx4[i]-crn5[i];
	    /*P3->C5*/
	    D9[i] = crn5[i]-crx3[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? (L1+L3)*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L3*L7/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L4*L5/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L6/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	
	area[2][1][1] = nb_flag[2][1][1] == YES? L2*L6-L8*L9/2.0: -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? 0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;
	

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case04_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],
		D7[3],D8[3],D9[3],D10[3],D11[3],D12[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 
	 = L8 = L9 = L10 =  L11 = L12 = 0.0;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_idir(cell3d,1,1); 
       	crx3 = crx_in_jdir(cell3d,0,1); 
	crx4 = crx_in_kdir(cell3d,1,0); 
	crx5 = crx_in_idir(cell3d,0,1); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C4->P2*/
	    D1[i] = crx2[i]-crn4[i];
	    /*C4->C2*/
	    D2[i] = crn2[i]-crn4[i];
	    /*C2->P5*/
	    D3[i] = crx5[i]-crn2[i];
	    /*C3->P1*/
	    D4[i] = crx1[i]-crn3[i];
	    /*C3->C1*/
	    D5[i] = crn1[i]-crn3[i];
	    /*C1->C5*/
	    D6[i] = crn5[i]-crn1[i];
	    /*C3->C4*/
	    D7[i] = crn4[i]-crn3[i];
	    /*C1->C2*/
	    D8[i] = crn2[i]-crn1[i];
	    /*C6->P5*/
	    D9[i] = crx5[i]-crn6[i];
	    /*C6->P4*/
	    D10[i] = crx4[i]-crn6[i];
	    /*C7->P1*/
	    D11[i] = crx1[i]-crn7[i];
	    /*C7->P3*/
	    D12[i] = crx3[i]-crn7[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L1+L3)*L2/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L7/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L8*L6-L9*L10/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L5*L6-L11*L12/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	
	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L7 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? 0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case05_comp1(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3],D15[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 
	= L10 =  L11 = L12 = L13 = L14 = L15 = 0.0;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0);
        crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);

	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    D1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    D2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C2->P6*/
	    D4[i] = crx6[i]-crn2[i];
	    /*C2->P5*/
	    D5[i] = crx5[i]-crn2[i];
	    /*C6->P5*/
	    D6[i] = crx5[i]-crn6[i];
	    /*C6->P2*/
	    D7[i] = crx2[i]-crn6[i];
	    /*C5->P2*/
	    D8[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D9[i] = crx3[i]-crn5[i];
	    /*C3->P4*/
	    D10[i] = crx4[i]-crn3[i];
	    /*C3->P1*/
	    D11[i] = crx1[i]-crn3[i];
	    /*C4->P4*/
	    D12[i] = crx4[i]-crn4[i];
	    /*C4->P6*/
	    D13[i] = crx6[i]-crn4[i];
	    /*C7->P1*/
	    D14[i] = crx1[i]-crn7[i];
	    /*C7->P3*/
	    D15[i] = crx3[i]-crn7[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	    L15 += fabs(D15[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L2-L4*L5/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L2*L3-L10*L11/2.0: -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L6*L7/2.0: -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L14*L15/2.0: -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? 0 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	
	area[2][1][1] = nb_flag[2][1][1] == YES? L1*L3-L8*L9/2.0: -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? 0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05: -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;
	
	area[0][1][1] = nb_flag[0][1][1] == YES? L12*L13/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case05_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3],D15[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 
	= L10 =  L11 = L12 = L13 = L14 = L15 = 0.0;
	
        crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,0);
        crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);

	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    D1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    D2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C2->P6*/
	    D4[i] = crx6[i]-crn2[i];
	    /*C2->P5*/
	    D5[i] = crx5[i]-crn2[i];
	    /*C3->P4*/
	    D6[i] = crx4[i]-crn3[i];
	    /*C3->P1*/
	    D7[i] = crx1[i]-crn3[i];
	    /*C5->P2*/
	    D8[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D9[i] = crx3[i]-crn5[i];
	    /*C7->P1*/
	    D10[i] = crx1[i]-crn7[i];
	    /*C7->P3*/
	    D11[i] = crx3[i]-crn7[i];
	    /*C6->P5*/
	    D12[i] = crx5[i]-crn6[i];
	    /*C6->P2*/
	    D13[i] = crx2[i]-crn6[i];
	    /*C4->P4*/
	    D14[i] = crx4[i]-crn4[i];
	    /*C4->P6*/
	    D15[i] = crx6[i]-crn4[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	    L15 += fabs(D15[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L4*L5/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L6*L7/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L1*L2-L10*L11/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L2*L3-L12*L13/2.0 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	
	area[0][1][1] = nb_flag[0][1][1] == YES? L1*L3-L14*L15/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L8*L9/2.0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case06_comp1(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn4[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double d1[3],d2[3],d3[3],D1[3],D2[3],D3[3],D4[3],D5[3],D6[3];
	double l1,l2,l3,L1,L2,L3,L4,L5,L6;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	l1 = l2 = l3 = 0.0;
	L1 = L2 = L3 = L4 = L5 = L6 = 0.0;

	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_kdir(cell3d,1,1);
	crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_jdir(cell3d,1,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_kdir(cell3d,1,0);

	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    d1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    d2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    d3[i] = crn1[i]-crn2[i];
	    /*C7->P1*/
	    D1[i] = crx1[i]-crn7[i];
	    /*C7->P2*/
	    D2[i] = crx2[i]-crn7[i];
	    /*C7->P3*/
	    D3[i] = crx3[i]-crn7[i];
	    /*C6->P4*/
	    D4[i] = crx4[i]-crn6[i];
	    /*C6->P5*/
	    D5[i] = crx5[i]-crn6[i];
	    /*C6->P6*/
	    D6[i] = crx6[i]-crn6[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L4*L5/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L5*L6/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L1*L3/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L1*L2/2.0 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L4*L6 + L2*L3)/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case06_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn4[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double d1[3],d2[3],d3[3],D1[3],D2[3],D3[3],D4[3],D5[3],D6[3];
	double l1,l2,l3,L1,L2,L3,L4,L5,L6;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

 	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
 	l1 = l2 = l3 = 0.0;
	L1 = L2 = L3 = L4 = L5 = L6 = 0.0;

	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_kdir(cell3d,1,1);
	crx3 = crx_in_jdir(cell3d,0,1);
	crx4 = crx_in_jdir(cell3d,1,1);
	crx5 = crx_in_idir(cell3d,0,1);
	crx6 = crx_in_kdir(cell3d,1,0);

	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    d1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    d2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    d3[i] = crn1[i]-crn2[i];
	    /*C7->P1*/
	    D1[i] = crx1[i]-crn7[i];
	    /*C7->P2*/
	    D2[i] = crx2[i]-crn7[i];
	    /*C7->P3*/
	    D3[i] = crx3[i]-crn7[i];
	    /*C6->P4*/
	    D4[i] = crx4[i]-crn6[i];
	    /*C6->P5*/
	    D5[i] = crx5[i]-crn6[i];
	    /*C6->P6*/
	    D6[i] = crx6[i]-crn6[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? l1*l2-L4*L5/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? l2*l3-L5*L6/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? l1*l2-L1*L3/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? l2*l3-L1*L2/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? 0 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? l1*l3-(L2*L3+L4*L6)/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? l1*l3 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? 0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void cell3d_area_case07_comp1(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3],D15[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 
	L10 = L11 = L12 = L13 = L14 = L15 = 0.0;

	crx1 = crx_in_kdir(cell3d,0,1); 
        crx2 = crx_in_jdir(cell3d,1,0);
        crx3 = crx_in_kdir(cell3d,1,1);
	crx4 = crx_in_idir(cell3d,0,1); 
	crx5 = crx_in_jdir(cell3d,0,1); 
	crx6 = crx_in_idir(cell3d,0,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->P2*/
	    D1[i] = crx2[i]-crn2[i];
	    /*C2->P4*/
	    D2[i] = crx4[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C1->P6*/
	    D4[i] = crx6[i]-crn1[i];
	    /*C5->P6*/
	    D5[i] = crx6[i]-crn5[i];
	    /*C5->P5*/
	    D6[i] = crx5[i]-crn5[i];
	    /*C7->P5*/
	    D7[i] = crx5[i]-crn7[i];
	    /*C7->P3*/
	    D8[i] = crx3[i]-crn7[i];
	    /*C5->C6*/
	    D9[i] = crn6[i]-crn5[i];
	    /*C6->C8*/
	    D10[i] = crn8[i]-crn6[i];
	    /*C8->C4*/
	    D11[i] = crn4[i]-crn8[i];
	    /*C8->P3*/
	    D12[i] = crx3[i]-crn8[i];
	    /*C4->P1*/
	    D13[i] = crx1[i]-crn4[i];
	    /*C4->P2*/
	    D14[i] = crx2[i]-crn4[i];
	    /*C1->C5*/
	    D15[i] = crn5[i]-crn1[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	    L15 += fabs(D15[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L10*L11-L1*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L15*L9-(L2+L4)*L3/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L5*L6/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L12+L13)*L11/2.0 : -1.0;
	
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? 0 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L9*L10-L7*L8/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? 0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L13*L14/2.0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case07_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 
	L10 = L11 = L12 = L13 = L14 = 0.0;

	crx1 = crx_in_kdir(cell3d,0,1); 
        crx2 = crx_in_jdir(cell3d,1,0);
        crx3 = crx_in_kdir(cell3d,1,1);
	crx4 = crx_in_idir(cell3d,0,1); 
	crx5 = crx_in_jdir(cell3d,0,1); 
	crx6 = crx_in_idir(cell3d,0,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->P2*/
	    D1[i] = crx2[i]-crn2[i];
	    /*C2->P4*/
	    D2[i] = crx4[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C1->P6*/
	    D4[i] = crx6[i]-crn1[i];
	    /*C1->C3*/
	    D5[i] = crn3[i]-crn1[i];
	    /*C4->P2*/
	    D6[i] = crx2[i]-crn4[i];
	    /*C4->P1*/
	    D7[i] = crx1[i]-crn4[i];
	    /*C5->P6*/
	    D8[i] = crx6[i]-crn5[i];
	    /*C5->P5*/
	    D9[i] = crx5[i]-crn5[i];
	    /*C1->C5*/
	    D10[i] = crn5[i]-crn1[i];
	    /*C3->C7*/
	    D11[i] = crn7[i]-crn3[i];
	    /*C3->P1*/
	    D12[i] = crx1[i]-crn3[i];
	    /*C7->P3*/
	    D13[i] = crx3[i]-crn7[i];
	    /*C7->P5*/
	    D14[i] = crx5[i]-crn7[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	}
	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L2+L4)*L3/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L5*L10-L8*L9/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L12+L13)*L11/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? 0 : -1.0;
	
	area[2][1][1] = nb_flag[2][1][1] == YES? L13*L14/2.0 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L3*L5-L6*L7/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0; 
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;

	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case08_comp1(
	CELL_INFO_3D *cell3d,
	double crn4[],
	double crn5[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3];
	double L1,L2,L3,L4,L5,L6,L7,L8;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = 0.0; 
	
	crx1 = crx_in_jdir(cell3d,0,1);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);

	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_jdir(cell3d,1,0);
	crx6 = crx_in_kdir(cell3d,1,1);
	crx7 = crx_in_jdir(cell3d,1,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C8->P7*/
	    D1[i] = crx7[i]-crn8[i];
	    /*C8->C4*/
	    D2[i] = crn4[i]-crn8[i];
	    /*C4->P5*/
	    D3[i] = crx5[i]-crn4[i];
	    /*C8->P6*/
	    D4[i] = crx6[i]-crn8[i];
	    /*C4->P4*/
	    D5[i] = crx4[i]-crn4[i];
	    /*C5->P2*/
	    D6[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D7[i] = crx3[i]-crn5[i];
	    /*C5->P1*/
	    D8[i] = crx1[i]-crn5[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? (L1+L3)*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L6*L7/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L7*L8/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L4+L5)*L2/2.0 : -1.0;
	
	area[1][2][2] = nb_flag[1][2][2] == YES? 0.0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L1*L4/2.0+L6*L8/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L3*L5/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;

}


LOCAL 	void cell3d_area_case08_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = L10 
		= L11 = L12 = L13 = L14 = 0.0; 
	
	crx1 = crx_in_jdir(cell3d,0,1);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);

	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_jdir(cell3d,1,0);
	crx6 = crx_in_kdir(cell3d,1,1);
	crx7 = crx_in_jdir(cell3d,1,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C1->C5*/
	    D1[i] = crn5[i]-crn1[i];
	    /*C1->C3*/
	    D2[i] = crn3[i]-crn1[i];
	    /*C1->C2*/
	    D3[i] = crn2[i]-crn1[i];
	    /*C2->P5*/
	    D4[i] = crx5[i]-crn2[i];
	    /*C6->P7*/
	    D5[i] = crx7[i]-crn6[i];
	    /*C5->P2*/
	    D6[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D7[i] = crx3[i]-crn5[i];
	    /*C5->P1*/
	    D8[i] = crx1[i]-crn5[i];
	    /*C3->P4*/
	    D9[i] = crx4[i]-crn3[i];
	    /*C7->P6*/
	    D10[i] = crx6[i]-crn7[i];
	    /*C4->P4*/
	    D11[i] = crx4[i]-crn4[i];
	    /*C4->P5*/
	    D12[i] = crx5[i]-crn4[i];
	    /*C8->P7*/
	    D13[i] = crx7[i]-crn8[i];
	    /*C8->P6*/
	    D14[i] = crx6[i]-crn8[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L4+L5)*L1/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L1*L3-L6*L7/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L1*L2-L7*L8/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L9+L10)*L1/2.0 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? 0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? 0 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L2*L3-L6*L8/2.0-L13*L14/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L3-L11*L12/2.0 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case09_comp1(
	CELL_INFO_3D *cell3d,
	double crn4[],
	double crn5[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3];
	double L1,L2,L3,L4,L5,L6;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = 0.0; 
	
	crx1 = crx_in_jdir(cell3d,0,1);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);

	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_jdir(cell3d,1,0);
	crx6 = crx_in_idir(cell3d,1,1);

	for (i = 0; i < 3; ++i)
	{	
	    /*C5->P2*/
	    D1[i] = crx2[i]-crn5[i];
	    /*C5->P1*/
	    D2[i] = crx1[i]-crn5[i];
	    /*C5->P3*/
	    D3[i] = crx3[i]-crn5[i];
	    /*C4->P6*/
	    D4[i] = crx6[i]-crn4[i];
	    /*C4->P5*/
	    D5[i] = crx5[i]-crn4[i];
	    /*C4->P4*/
	    D6[i] = crx4[i]-crn4[i];
	}
	
	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L4*L5/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L1*L3/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L2*L3/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L4*L6/2.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L1*L2/2.0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L5*L6/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case09_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 0.0; 

	crx1 = crx_in_jdir(cell3d,0,1);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);

	crx4 = crx_in_kdir(cell3d,0,1);
	crx5 = crx_in_jdir(cell3d,1,0);
	crx6 = crx_in_idir(cell3d,1,1);

	for (i = 0; i < 3; ++i)
	{
	    /*C1->C5*/
	    D1[i] = crn5[i]-crn1[i];
	    /*C1->C3*/
	    D2[i] = crn3[i]-crn1[i];
	    /*C1->C2*/
	    D3[i] = crn2[i]-crn1[i];
	    /*C4->P6*/
	    D4[i] = crx6[i]-crn4[i];
	    /*C4->P5*/
	    D5[i] = crx5[i]-crn4[i];
	    /*C4->P4*/
	    D6[i] = crx4[i]-crn4[i];
	    /*C5->P2*/
	    D7[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D8[i] = crx3[i]-crn5[i];
	    /*C5->P1*/
	    D9[i] = crx1[i]-crn5[i];
	}
	
	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);	
	    L2 += fabs(D2[i]);	
	    L3 += fabs(D3[i]);	
	    L4 += fabs(D4[i]);	
	    L5 += fabs(D5[i]);	
	    L6 += fabs(D6[i]);	
	    L7 += fabs(D7[i]);	
	    L8 += fabs(D8[i]);	
	    L9 += fabs(D9[i]);	
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? L1*L2-L4*L5/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L1*L3-L7*L8/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L1*L2-L8*L9/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L1*L3-L4*L6/2.0 : -1.0;
	area[1][0][2] = nb_flag[1][0][2] == YES? 0.0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? 0.0 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L2*L3-L7*L9/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? 0 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L3-L5*L6/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}


LOCAL	void cell3d_area_case10_comp1(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3],D15[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 
	L10 = L11 = L12 = L13 = L14 = L15 = 0.0;
	
	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);
	crx4 = crx_in_jdir(cell3d,1,1);
	crx5 = crx_in_kdir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    D1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    D2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C4->P6*/
	    D4[i] = crx6[i]-crn4[i];
	    /*C4->C8*/
	    D5[i] = crn8[i]-crn4[i];
	    /*C8->P4*/
	    D6[i] = crx4[i]-crn8[i];
	    /*C5->P2*/
	    D7[i] = crx2[i]-crn5[i];
	    /*C5->P3*/
	    D8[i] = crx3[i]-crn5[i];
	    /*C5->C7*/
	    D9[i] = crn7[i]-crn5[i];
	    /*C7->P1*/
	    D10[i] = crx1[i]-crn7[i];
	    /*C3->p5*/
	    D11[i] = crx5[i]-crn3[i];
	    /*C3->P1*/
	    D12[i] = crx1[i]-crn3[i];
	    /*C6->P4*/
	    D13[i] = crx4[i]-crn6[i];
	    /*C6->P2*/
	    D14[i] = crx2[i]-crn6[i];
	    /*C4->P5*/
	    D15[i] = crx5[i]-crn4[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	    L15 += fabs(D15[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L4+L6)*L5/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L7*L8/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L8+L10)*L9/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L2*L3-L11*L12/2.0 : -1.0;
	
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? 0 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L1*L3-L13*L14/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? 0 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0 : -1.0;

	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L14*L15/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL 	void cell3d_area_case10_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3],D11[3],D12[3],D13[3],D14[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 
	L10 = L11 = L12 = L13 = L14 = 0.0;
	
	crx1 = crx_in_idir(cell3d,1,0);
	crx2 = crx_in_kdir(cell3d,1,0);
	crx3 = crx_in_idir(cell3d,0,0);
	crx4 = crx_in_jdir(cell3d,1,1);
	crx5 = crx_in_kdir(cell3d,0,1);
	crx6 = crx_in_jdir(cell3d,1,0);
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->C4*/
	    D1[i] = crn4[i]-crn2[i];
	    /*C2->C6*/
	    D2[i] = crn6[i]-crn2[i];
	    /*C2->C1*/
	    D3[i] = crn1[i]-crn2[i];
	    /*C2->P6*/
	    D4[i] = crx6[i]-crn2[i];
	    /*C6->P4*/
	    D5[i] = crx4[i]-crn6[i];
	    /*C5->P3*/
	    D6[i] = crx3[i]-crn5[i];
	    /*C5->P2*/
	    D7[i] = crx2[i]-crn5[i];
	    /*C1->P3*/
	    D8[i] = crx3[i]-crn1[i];
	    /*C3->P1*/
	    D9[i] = crx1[i]-crn3[i];
	    /*C1->C3*/
	    D10[i] = crn3[i]-crn1[i];
	    /*C3->P5*/
	    D11[i] = crx5[i]-crn3[i];
	    /*C6->P2*/
	    D12[i] = crx2[i]-crn6[i];
	    /*C4->P5*/
	    D13[i] = crx5[i]-crn4[i];
	    /*C4->P6*/
	    D14[i] = crx6[i]-crn4[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? (L4+L5)*L2/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L2*L3-L6*L7/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L8+L9)*L10/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? L9*L11/2.0 : -1.0;
	
	area[1][0][2] = nb_flag[1][0][2] == YES? 0 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? L5*L12/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L1*L3-L13*L14/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case11_comp1(
	CELL_INFO_3D *cell3d,
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],
		D9[3],D10[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = L10 = 0.0;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_jdir(cell3d,0,0); 
        crx3 = crx_in_jdir(cell3d,1,0); 
	crx4 = crx_in_idir(cell3d,1,1);
	crx5 = crx_in_jdir(cell3d,0,1); 
        crx6 = crx_in_idir(cell3d,0,1); 
        crx7 = crx_in_idir(cell3d,0,0); 
	crx8 = crx_in_jdir(cell3d,1,1); 

	for (i = 0; i < 3; ++i)
	{
	    /*C4->P4*/
	    D1[i] = crx4[i]-crn4[i];
	    /*C4->P3*/
	    D2[i] = crx3[i]-crn4[i];
	    /*C4->C3*/
	    D3[i] = crn3[i]-crn4[i];
	    /*C3->P1*/
	    D4[i] = crx1[i]-crn3[i];
	    /*C3->P2*/
	    D5[i] = crx2[i]-crn3[i];
	    /*C6->P8*/
	    D6[i] = crx8[i]-crn6[i];
	    /*C6->P6*/
	    D7[i] = crx6[i]-crn6[i];
	    /*C6->C5*/
	    D8[i] = crn5[i]-crn6[i];
	    /*C5->P5*/
	    D9[i] = crx5[i]-crn5[i];
	    /*C5->P7*/
	    D10[i] = crx7[i]-crn5[i];
	}
	
	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	}
	
	area[1][1][2] = nb_flag[1][1][2] == YES? (L1*L2+L6*L7)/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L7+L10)*L8/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L4*L5+L9*L10)/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1+L4)*L3/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L6+L9)*L8/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0.0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? (L2+L5)*L3/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case11_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8;
	double d1[3],d2[3],d3[3];
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3],
		D10[3],D11[3],D12[3],D13[3],D14[3],D15[3],D16[3];
	double l1,l2,l3;
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	l1 = l2 = l3 = 0.0;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = L10 = 	
		L11 = L12 = L13 = L14 = L15 = L16 = 0.0;

	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_jdir(cell3d,0,0); 
        crx3 = crx_in_jdir(cell3d,1,0); 
	crx4 = crx_in_idir(cell3d,1,1);
	crx5 = crx_in_jdir(cell3d,0,1); 
        crx6 = crx_in_idir(cell3d,0,1); 
        crx7 = crx_in_idir(cell3d,0,0); 
	crx8 = crx_in_jdir(cell3d,1,1); 

	for (i = 0; i < 3; ++i)
	{
	    /* C1->C5 */
	    d1[i] = crn5[i]-crn1[i];
	    /* C1->C3 */
	    d2[i] = crn3[i]-crn1[i];
	    /* C1->C2 */
	    d3[i] = crn2[i]-crn1[i];
	    /* C4->P4 */
	    D1[i] = crx4[i]-crn4[i];
	    /* C4->P3 */
	    D2[i] = crx3[i]-crn4[i];
	    /* C6->P8 */
	    D3[i] = crx8[i]-crn6[i];
	    /* C6->P6 */
	    D4[i] = crx6[i]-crn6[i];
	    /* C3->P2 */
	    D5[i] = crx2[i]-crn3[i];
	    /* C3->P1 */
	    D6[i] = crx1[i]-crn3[i];
	    /* C5->P7 */
	    D7[i] = crx7[i]-crn5[i];
	    /* C5->P5 */
	    D8[i] = crx5[i]-crn5[i];
	    /* C2->P6 */
	    D9[i] = crx6[i]-crn2[i];
	    /* C1->P7 */
	    D10[i] = crx7[i]-crn1[i];
	    /* C2->P3 */
	    D11[i] = crx3[i]-crn2[i];
	    /* C1->P2 */
	    D12[i] = crx2[i]-crn1[i];
	    /* C8->P8 */
	    D13[i] = crx8[i]-crn8[i];
	    /* C7->P5 */
	    D14[i] = crx5[i]-crn7[i];
	    /* C8->P4 */
	    D15[i] = crx4[i]-crn8[i];
	    /* C7->P1 */
	    D16[i] = crx1[i]-crn7[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	    L13 += fabs(D13[i]);
	    L14 += fabs(D14[i]);
	    L15 += fabs(D15[i]);
	    L16 += fabs(D16[i]);
	}
	
	area[1][1][2] = 
		nb_flag[1][1][2] == YES? (l1*l2-L1*L2/2.0-L3*L4/2.0) : -1.0;	
	area[1][0][1] = nb_flag[1][0][1] == YES? (L9+L10)*l3/2.0 : -1.0;
	area[1][1][0] = 
		nb_flag[1][1][0] == YES? (l1*l2-L5*L6/2.0-L7*L8/2.0) : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L15+L16)*l3/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L13+L14)*l3/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? 0.0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? (L11+L12)*l3/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? -0.05 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL 	void cell3d_area_case12_comp1(
	CELL_INFO_3D *cell3d,
	double crn4[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8,*crx9;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 0.0; 	
	
	/* the first corner_tri */
	crx1 = crx_in_idir(cell3d,1,0);
        crx2 = crx_in_kdir(cell3d,1,1); 
        crx3 = crx_in_jdir(cell3d,0,1);

	/* the second corner_tri */
        crx4 = crx_in_kdir(cell3d,0,1);
        crx5 = crx_in_jdir(cell3d,1,0); 
        crx6 = crx_in_idir(cell3d,1,1); 

	/* the third corner_tri */
        crx7 = crx_in_jdir(cell3d,1,1); 
        crx8 = crx_in_idir(cell3d,0,1); 
        crx9 = crx_in_kdir(cell3d,1,0); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C4->P6*/
	    D1[i] = crx6[i]-crn4[i];
	    /*C4->P4*/
	    D2[i] = crx4[i]-crn4[i];
	    /*C4->P5*/
	    D3[i] = crx5[i]-crn4[i];
	    /*C7->P2*/
	    D4[i] = crx2[i]-crn7[i];
	    /*C7->P1*/
	    D5[i] = crx1[i]-crn7[i];
	    /*C7->P3*/
	    D6[i] = crx3[i]-crn7[i];
	    /*C6->P7*/
	    D7[i] = crx7[i]-crn6[i];
	    /*C6->P9*/
	    D8[i] = crx9[i]-crn6[i];
	    /*C6->P8*/
	    D9[i] = crx8[i]-crn6[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L1*L3+L7*L9)/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? L8*L9/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? L5*L6/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1*L2+L4*L5)/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L7*L8+L4*L6)/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L2*L3/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case12_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn4[],
	double crn6[],
	double crn7[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8,*crx9;
	double d1[3],d2[3],d3[3];
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3];
	double l1,l2,l3,L1,L2,L3,L4,L5,L6,L7,L8,L9;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	l1 = l2 = l3 = L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = 0.0; 	
	
	/* the first corner_tri */
	crx1 = crx_in_idir(cell3d,1,0);
        crx2 = crx_in_kdir(cell3d,1,1); 
        crx3 = crx_in_jdir(cell3d,0,1);

	/* the second corner_tri */
        crx4 = crx_in_kdir(cell3d,0,1);
        crx5 = crx_in_jdir(cell3d,1,0); 
        crx6 = crx_in_idir(cell3d,1,1); 

	/* the third corner_tri */
        crx7 = crx_in_jdir(cell3d,1,1); 
        crx8 = crx_in_idir(cell3d,0,1); 
        crx9 = crx_in_kdir(cell3d,1,0); 
	
	for (i = 0; i < 3; ++i)
	{
	    /*C2->C6*/
	    d1[i] = crn6[i]-crn2[i];
	    /*C2->C4*/
	    d2[i] = crn4[i]-crn2[i];
	    /*C2->C1*/
	    d3[i] = crn1[i]-crn2[i];
	    /*C4->P6*/
	    D1[i] = crx6[i]-crn4[i];
	    /*C4->P4*/
	    D2[i] = crx4[i]-crn4[i];
	    /*C4->P5*/
	    D3[i] = crx5[i]-crn4[i];
	    /*C7->P2*/
	    D4[i] = crx2[i]-crn7[i];
	    /*C7->P1*/
	    D5[i] = crx1[i]-crn7[i];
	    /*C7->P3*/
	    D6[i] = crx3[i]-crn7[i];
	    /*C6->P7*/
	    D7[i] = crx7[i]-crn6[i];
	    /*C6->P9*/
	    D8[i] = crx9[i]-crn6[i];
	    /*C6->P8*/
	    D9[i] = crx8[i]-crn6[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	}

	area[1][1][2] = 
		nb_flag[1][1][2] == YES? l1*l2-L1*L3/2.0-L7*L9/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? l1*l3-L8*L9/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? l1*l2-L5*L6/2.0 : -1.0;
	area[1][2][1] = 
		nb_flag[1][2][1] == YES? l1*l3-L1*L2/2.0-L4*L5/2.0 : -1.0;
	
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? 0.0 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = 
		nb_flag[2][1][1] == YES? l2*l3-L7*L8/2.0-L4*L6/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? l2*l3-L2*L3/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0.0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0.0 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;


	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case13_comp1(
	CELL_INFO_3D *cell3d,
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3],
		D10[3],D11[3],D12[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = L10 = L11 = L12 = 0.0; 	
	
	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,1);
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 
	crx5 = crx_in_idir(cell3d,0,1);
        crx6 = crx_in_kdir(cell3d,0,1); 
        crx7 = crx_in_jdir(cell3d,1,0); 
        crx8 = crx_in_idir(cell3d,1,1); 

	for (i = 0; i < 3; ++i)
	{
	    /*C6->P4*/
	    D1[i] = crx4[i]-crn6[i];
	    /*C6->P5*/
	    D2[i] = crx5[i]-crn6[i];
	    /*C6->C5*/
	    D3[i] = crn5[i]-crn6[i];
	    /*C5->P3*/
	    D4[i] = crx3[i]-crn5[i];
	    /*C7->P1*/
	    D5[i] = crx1[i]-crn7[i];
	    /*C5->C7*/
	    D6[i] = crn7[i]-crn5[i];
	    /*C7->P2*/
	    D7[i] = crx2[i]-crn7[i];
	    /*C8->P4*/
	    D8[i] = crx4[i]-crn8[i];
	    /*C8->P2*/
	    D9[i] = crx2[i]-crn8[i];
	    /*C4->P8*/
	    D10[i] = crx8[i]-crn4[i];
	    /*C4->P6*/
	    D11[i] = crx6[i]-crn4[i];
	    /*C4->P7*/
	    D12[i] = crx7[i]-crn4[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L1*L2+L10*L12)/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L2+L4)*L3/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L4+L5)*L6/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L5*L7+L10*L11)/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L3*L6-L8*L9/2.0) : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? 0 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? 0 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;
	area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? L11*L12/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case13_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn4[],
	double crn5[],
	double crn6[],
	double crn7[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8;
	double d1[3],d2[3],d3[3];
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3],
		D10[3],D11[3],D12[3];
	double l1,l2,l3,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	l1 = l2 = l3 = L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 
		= L9 = L10 = L11 = L12 = 0.0; 	
	
	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_kdir(cell3d,1,1);
        crx3 = crx_in_idir(cell3d,0,0); 
	crx4 = crx_in_jdir(cell3d,1,1); 
	crx5 = crx_in_idir(cell3d,0,1);
        crx6 = crx_in_kdir(cell3d,0,1); 
        crx7 = crx_in_jdir(cell3d,1,0); 
        crx8 = crx_in_idir(cell3d,1,1); 

	for (i = 0; i < 3; ++i)
	{
	    /*C1->C5*/
	    d1[i] = crn5[i]-crn1[i];
	    /*C1->C3*/
	    d2[i] = crn3[i]-crn1[i];
	    /*C1->C2*/
	    d3[i] = crn2[i]-crn1[i];
	    /*C6->P4*/
	    D1[i] = crx4[i]-crn6[i];
	    /*C6->P5*/
	    D2[i] = crx5[i]-crn6[i];
	    /*C4->P8*/
	    D3[i] = crx8[i]-crn4[i];
	    /*C4->P7*/
	    D4[i] = crx7[i]-crn4[i];
	    /*C2->P5*/
	    D5[i] = crx5[i]-crn2[i];
	    /*C1->P3*/
	    D6[i] = crx3[i]-crn1[i];
	    /*C3->P1*/
	    D7[i] = crx1[i]-crn3[i];
	    /*C4->P6*/
	    D8[i] = crx6[i]-crn4[i];
	    /*C8->P4*/
	    D9[i] = crx4[i]-crn8[i];
	    /*C8->P2*/
	    D10[i] = crx2[i]-crn8[i];
	    /*C7->P2*/
	    D11[i] = crx2[i]-crn7[i];
	    /*C7->P1*/
	    D12[i] = crx1[i]-crn7[i];
	}
	
	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);
	    L2 += fabs(D2[i]);
	    L3 += fabs(D3[i]);
	    L4 += fabs(D4[i]);
	    L5 += fabs(D5[i]);
	    L6 += fabs(D6[i]);
	    L7 += fabs(D7[i]);
	    L8 += fabs(D8[i]);
	    L9 += fabs(D9[i]);
	    L10 += fabs(D10[i]);
	    L11 += fabs(D11[i]);
	    L12 += fabs(D12[i]);
	}

	area[1][1][2] = 
		nb_flag[1][1][2] == YES? l1*l2-L1*L2/2.0-L3*L4/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L5+L6)*l3/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L6+L7)*l2/2.0 : -1.0;
	area[1][2][1] = 
		nb_flag[1][2][1] == YES? l1*l3-L3*L8/2.0-L11*L12/2.0 : -1.0;
	
	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;
	 
	area[2][1][1] = nb_flag[2][1][1] == YES? L9*L10/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? l2*l3-L4*L8/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? 0 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? 0 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	 
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;
	area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}


LOCAL	void cell3d_area_case14_comp1(
	CELL_INFO_3D *cell3d,
	double crn2[],
	double crn3[],
	double crn5[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8,
			*crx9,*crx10,*crx11,*crx12;
	double D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],D7[3],D8[3],D9[3],
		D10[3],D11[3],D12[3];
	double L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 = L9 = L10 = L11 = L12 = 0.0; 	
	        
	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_jdir(cell3d,0,0); 
        crx3 = crx_in_kdir(cell3d,0,1); 
        crx4 = crx_in_kdir(cell3d,1,1); 
        crx5 = crx_in_idir(cell3d,1,1); 
        crx6 = crx_in_jdir(cell3d,1,1); 
        crx7 = crx_in_jdir(cell3d,0,1); 
        crx8 = crx_in_kdir(cell3d,1,0); 
        crx9 = crx_in_idir(cell3d,0,0); 
        crx10 = crx_in_jdir(cell3d,1,0); 
        crx11 = crx_in_kdir(cell3d,0,0); 
        crx12 = crx_in_idir(cell3d,0,1); 
	
	for (i = 0 ; i < 3; ++i)
	{
	    /*C3->P1*/
	    D1[i] = crx1[i]-crn3[i];
	    /*C3->P2*/
	    D2[i] = crx2[i]-crn3[i];
	    /*C3->P3*/
	    D3[i] = crx3[i]-crn3[i];
	    /*C8->P4*/
	    D4[i] = crx4[i]-crn8[i];
	    /*C8->P5*/
	    D5[i] = crx5[i]-crn8[i];
	    /*C8->P6*/
	    D6[i] = crx6[i]-crn8[i];
	    /*C5->P7*/
	    D7[i] = crx7[i]-crn5[i];
	    /*C5->P8*/
	    D8[i] = crx8[i]-crn5[i];
	    /*C5->P9*/
	    D9[i] = crx9[i]-crn5[i];
	    /*C2->P10*/
	    D10[i] = crx10[i]-crn2[i];
	    /*C2->P11*/
	    D11[i] = crx11[i]-crn2[i];
	    /*C2->P12*/
	    D12[i] = crx12[i]-crn2[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    L1 += fabs(D1[i]);		
	    L2 += fabs(D2[i]);		
	    L3 += fabs(D3[i]);		
	    L4 += fabs(D4[i]);		
	    L5 += fabs(D5[i]);		
	    L6 += fabs(D6[i]);		
	    L7 += fabs(D7[i]);		
	    L8 += fabs(D8[i]);		
	    L9 += fabs(D9[i]);		
	    L10 += fabs(D10[i]);		
	    L11 += fabs(D11[i]);		
	    L12 += fabs(D12[i]);		
	}

	area[1][1][2] = nb_flag[1][1][2] == YES? (L5*L6+L10*L12)/2.0 : -1.0;
	area[1][0][1] = nb_flag[1][0][1] == YES? (L8*L9+L11*L12)/2.0 : -1.0;
	area[1][1][0] = nb_flag[1][1][0] == YES? (L7*L9+L1*L2)/2.0 : -1.0;
	area[1][2][1] = nb_flag[1][2][1] == YES? (L1*L3+L4*L5)/2.0 : -1.0;

	area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	area[2][1][1] = nb_flag[2][1][1] == YES? (L4*L6+L7*L8)/2.0 : -1.0;
	area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	area[2][2][2] = nb_flag[2][2][2] == YES? -0.1 : -1.0;
	area[2][0][0] = nb_flag[2][0][0] == YES? -0.1 : -1.0;

	area[0][1][1] = nb_flag[0][1][1] == YES? (L10*L11+L2*L3)/2.0 : -1.0;
	area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	area[0][0][1] = nb_flag[0][0][1] == YES? -0.05 : -1.0;
	area[0][1][0] = nb_flag[0][1][0] == YES? -0.05 : -1.0;
	area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	area[0][2][0] = nb_flag[0][2][0] == YES? -0.1 : -1.0;
	area[0][0][2] = nb_flag[0][0][2] == YES? -0.1 : -1.0;


	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL	void cell3d_area_case14_comp0(
	CELL_INFO_3D *cell3d,
	double crn1[],
	double crn2[],
	double crn3[],
	double crn5[],
	double crn8[],
	double vol)
{
	int i,j,k;
	const double *crx1,*crx2,*crx3,*crx4,*crx5,*crx6,*crx7,*crx8,
			*crx9,*crx10,*crx11,*crx12;
	double d1[3],d2[3],d3[3],D1[3],D2[3],D3[3],D4[3],D5[3],D6[3],
		D7[3],D8[3],D9[3],D10[3],D11[3],D12[3];
	double l1,l2,l3,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12;
	double (*area)[3][3];
	boolean (*nb_flag)[3][3];
	double max_area = -0.5;
	int max_index[3];
	boolean find_nb = NO;

	area = cell3d->area;
	nb_flag = cell3d->nb_flag;
	l1 = l2 = l3 = L1 = L2 = L3 = L4 = L5 = L6 = L7 = L8 
		= L9 = L10 = L11 = L12 = 0.0; 	
	        
	crx1 = crx_in_idir(cell3d,1,0); 
        crx2 = crx_in_jdir(cell3d,0,0); 
        crx3 = crx_in_kdir(cell3d,0,1); 
        crx4 = crx_in_kdir(cell3d,1,1); 
        crx5 = crx_in_idir(cell3d,1,1); 
        crx6 = crx_in_jdir(cell3d,1,1); 
        crx7 = crx_in_jdir(cell3d,0,1); 
        crx8 = crx_in_kdir(cell3d,1,0); 
        crx9 = crx_in_idir(cell3d,0,0); 
        crx10 = crx_in_jdir(cell3d,1,0); 
        crx11 = crx_in_kdir(cell3d,0,0); 
        crx12 = crx_in_idir(cell3d,0,1); 
	
	for (i = 0 ; i < 3; ++i)
	{
	    /*C1->C5*/
	    d1[i] = crn5[i]-crn1[i];
	    /*C1->C3*/
	    d2[i] = crn3[i]-crn1[i];
	    /*C1->C2*/
	    d3[i] = crn2[i]-crn1[i];
	    /*C3->P1*/
	    D1[i] = crx1[i]-crn3[i];
	    /*C3->P2*/
	    D2[i] = crx2[i]-crn3[i];
	    /*C3->P3*/
	    D3[i] = crx3[i]-crn3[i];
	    /*C8->P4*/
	    D4[i] = crx4[i]-crn8[i];
	    /*C8->P5*/
	    D5[i] = crx5[i]-crn8[i];
	    /*C8->P6*/
	    D6[i] = crx6[i]-crn8[i];
	    /*C5->P7*/
	    D7[i] = crx7[i]-crn5[i];
	    /*C5->P8*/
	    D8[i] = crx8[i]-crn5[i];
	    /*C5->P9*/
	    D9[i] = crx9[i]-crn5[i];
	    /*C2->P10*/
	    D10[i] = crx10[i]-crn2[i];
	    /*C2->P11*/
	    D11[i] = crx11[i]-crn2[i];
	    /*C2->P12*/
	    D12[i] = crx12[i]-crn2[i];
	}

	for (i = 0; i < 3; ++i)
	{
	    l1 += fabs(d1[i]);
	    l2 += fabs(d2[i]);
	    l3 += fabs(d3[i]);
	    L1 += fabs(D1[i]);		
	    L2 += fabs(D2[i]);		
	    L3 += fabs(D3[i]);		
	    L4 += fabs(D4[i]);		
	    L5 += fabs(D5[i]);		
	    L6 += fabs(D6[i]);		
	    L7 += fabs(D7[i]);		
	    L8 += fabs(D8[i]);		
	    L9 += fabs(D9[i]);		
	    L10 += fabs(D10[i]);		
	    L11 += fabs(D11[i]);		
	    L12 += fabs(D12[i]);		
	}

	for (i = 0; i < 3; ++i)
	{
	    area[1][1][2] = 
	    	nb_flag[1][1][2] == YES? (l1*l2-L10*L12/2.0-L5*L6/2.0) : -1.0;
	    area[1][0][1] = 
	    	nb_flag[1][0][1] == YES? (l1*l3-L11*L12/2.0-L8*L9/2.0) : -1.0;
	    area[1][1][0] = 
	    	nb_flag[1][1][0] == YES? (l1*l2-L7*L9/2.0-L1*L2/2.0) : -1.0;
	    area[1][2][1] =
	    	nb_flag[1][2][1] == YES? (l1*l3-L1*L3/2.0-L4*L5/2.0) : -1.0;
	    area[1][0][2] = nb_flag[1][0][2] == YES? -0.05 : -1.0;
	    area[1][0][0] = nb_flag[1][0][0] == YES? -0.05 : -1.0;
	    area[1][2][0] = nb_flag[1][2][0] == YES? -0.05 : -1.0;
	    area[1][2][2] = nb_flag[1][2][2] == YES? -0.05 : -1.0;

	    area[2][1][1] = 
	    	nb_flag[2][1][1] == YES? (l2*l3-L4*L6/2.0-L7*L8/2.0) : -1.0;
	    area[2][1][2] = nb_flag[2][1][2] == YES? -0.05 : -1.0;
	    area[2][0][1] = nb_flag[2][0][1] == YES? -0.05 : -1.0;
	    area[2][1][0] = nb_flag[2][1][0] == YES? -0.05 : -1.0;
	    area[2][2][1] = nb_flag[2][2][1] == YES? -0.05 : -1.0;
	    area[2][0][2] = nb_flag[2][0][2] == YES? -0.1 : -1.0;
	    area[2][2][0] = nb_flag[2][2][0] == YES? -0.1 : -1.0;

	    area[0][1][1] = 
	    	nb_flag[0][1][1] == YES? (l2*l3-L10*L11/2.0-L2*L3/2.0) : -1.0;
	    area[0][1][2] = nb_flag[0][1][2] == YES? -0.05 : -1.0;
	    area[0][0][1] = nb_flag[0][0][1] == YES? -0.05 : -1.0;
	    area[0][1][0] = nb_flag[0][1][0] == YES? -0.05 : -1.0;
	    area[0][2][1] = nb_flag[0][2][1] == YES? -0.05 : -1.0;
	    area[0][2][2] = nb_flag[0][2][2] == YES? -0.1 : -1.0;
	    area[0][0][0] = nb_flag[0][0][0] == YES? -0.1 : -1.0;
	}

	
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 3; ++j)
	for (k = 0; k < 3; ++k)
	{
	    if (area[i][j][k] > max_area)
	    {
	    	find_nb = YES;
		max_area = area[i][j][k];
		max_index[0] = i;
		max_index[1] = j;
		max_index[2] = k;
	    }
	}

	if (find_nb)
	    cell3d->nb_frac[max_index[0]][max_index[1]][max_index[2]] = vol;
	else
	    cell3d->orphan = vol;
}

LOCAL void gview_cell3d_skeleton(FILE*,const CELL_INFO_3D*);
LOCAL void gview_cell3d_points(FILE*,double**,int,const char*);
LOCAL void FT_GetTrisInCell(INTERFACE*,const int*,const int*,TRI***,int*);
LOCAL void gview_cell3d_tris(FILE*,TRI**,int);

EXPORT	void gview_plot_intfc_within_cell3d(
	char *out_name,
	INTERFACE *grid_intfc,
	const int *icoords,
	const int *shift,
	const CELL_INFO_3D *cell3d,
	int step)
{
	int i,j,k,num1,num2,num3;
	FILE* file;
	static const char *indent = "	";
	char fname[256],dname[256],path[256];
	const double *L,*U;
	static const char *color1 = "0 0 1 1",
			  *color2 = "0 1 0 1",
			  *color3 = "1 0 0 1";
	int comp3_count,comp0_count,crxing_count;
	const int (*comp)[2][2] = cell3d->comp;
	const double (*crx_coords)[2][2][3] = cell3d->crx_coords;
	double **comp3_coords,**comp0_coords,**crxing_coords;
	TRI **tris;
	int num_tris;

	comp3_count = comp0_count = crxing_count = 0;

	for (i = 0; i < 2; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	{
	    if (comp[i][j][k] == 3)
	    	comp3_count++;
	    else if (comp[i][j][k] == 0)
	    	comp0_count++;
	}

	for (i = 0; i < 3; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	{
	    if (crx_coords[i][j][k][0] != HUGE)
	    	crxing_count++;
	}
	FT_MatrixMemoryAlloc((POINTER*)&comp3_coords,
			comp3_count,3,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&comp0_coords,
			comp0_count,3,sizeof(double));
	FT_MatrixMemoryAlloc((POINTER*)&crxing_coords,
			crxing_count,3,sizeof(double));
	num1 = num2 = 0;
	for (i = 0; i < 2; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	{
	    if (comp[i][j][k] == 3)
	    {
	    	comp3_coords[num1][0] = cell3d->crds[i][j][k][0];
		comp3_coords[num1][1] = cell3d->crds[i][j][k][1];
		comp3_coords[num1][2] = cell3d->crds[i][j][k][2];
		++num1;
	    }
	    else if (comp[i][j][k] == 0)
	    {

	    	comp0_coords[num2][0] = cell3d->crds[i][j][k][0];
		comp0_coords[num2][1] = cell3d->crds[i][j][k][1];
		comp0_coords[num2][2] = cell3d->crds[i][j][k][2];
		++num2;
	    }
	}

	num3 = 0; 
	for (i = 0; i < 3; ++i)
	for (j = 0; j < 2; ++j)
	for (k = 0; k < 2; ++k)
	{
	    if (crx_coords[i][j][k][0] != HUGE)
	    {
	    	crxing_coords[num3][0] = crx_coords[i][j][k][0];
	    	crxing_coords[num3][1] = crx_coords[i][j][k][1];
	    	crxing_coords[num3][2] = crx_coords[i][j][k][2];
		++num3;
	    }
	}
	
	L = cell3d->crds[0][0][0];
	U = cell3d->crds[1][1][1];
	sprintf(dname,"%s/gview/gv.ts%s",out_name,right_flush(step,7));
	mkdir(dname,0777);

	sprintf(fname,"cell-%d-%d-%d.off",icoords[0],icoords[1],icoords[2]);
	sprintf(path,"%s/%s",dname,fname);
	if ((file = fopen(path,"w")) == NULL)
	{
	    screen("Cannot open %s in gview_plot_intfc_within_cell3d()\n",path);
	    return;
	}
	fprintf(file,"%s{LIST%s\n",indent,indent);
	gview_bounding_box(file,L,U,0,indent);
	gview_cell3d_skeleton(file,cell3d);
	gview_cell3d_points(file,comp3_coords,num1,color1);
	gview_cell3d_points(file,comp0_coords,num2,color2);
	gview_cell3d_points(file,crxing_coords,num3,color3);
	FT_GetTrisInCell(grid_intfc,icoords,shift,&tris,&num_tris);	
	gview_cell3d_tris(file,tris,num_tris);
	fprintf(file,"%s}\n",indent);
	fclose(file);
}

LOCAL	void gview_cell3d_skeleton(
	FILE *file,
	const CELL_INFO_3D *cell3d)
{
	const double (*crds)[2][2][3] = cell3d->crds;
	const char *indent = "	";	
	int i;

	fprintf(file, "{\n%sVECT\n",indent);
	fprintf(file, "12 24 1\n");
	fprintf(file, "2 2 2 2 2 2 2 2 2 2 2 2\n");
	fprintf(file, "1 0 0 0 0 0 0 0 0 0 0 0\n\n");

	fprintf(file, "%f %f %f%s",crds[0][0][0][0],
		crds[0][0][0][1],crds[0][0][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[1][0][0][0],
		crds[1][0][0][1],crds[1][0][0][2]);
	fprintf(file, "%f %f %f%s",crds[1][0][0][0],
		crds[1][0][0][1],crds[1][0][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[1][0][1][0],
		crds[1][0][1][1],crds[1][0][1][2]);
	fprintf(file, "%f %f %f%s",crds[1][0][1][0],
		crds[1][0][1][1],crds[1][0][1][2],indent);
	fprintf(file, "%f %f %f\n",crds[0][0][1][0],
		crds[0][0][1][1],crds[0][0][1][2]);
	fprintf(file, "%f %f %f%s",crds[0][0][1][0],
		crds[0][0][1][1],crds[0][0][1][2],indent);
	fprintf(file, "%f %f %f\n\n",crds[0][0][0][0],
		crds[0][0][0][1],crds[0][0][0][2]);
	
	fprintf(file, "%f %f %f%s",crds[0][1][0][0],
		crds[0][1][0][1],crds[0][1][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[1][1][0][0],
		crds[1][1][0][1],crds[1][1][0][2]);
	fprintf(file, "%f %f %f%s",crds[1][1][0][0],
		crds[1][1][0][1],crds[1][1][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[1][1][1][0],
		crds[1][1][1][1],crds[1][1][1][2]);
	fprintf(file, "%f %f %f%s",crds[1][1][1][0],
		crds[1][1][1][1],crds[1][1][1][2],indent);
	fprintf(file, "%f %f %f\n",crds[0][1][1][0],
		crds[0][1][1][1],crds[0][1][1][2]);
	fprintf(file, "%f %f %f%s",crds[0][1][1][0],
		crds[0][1][1][1],crds[0][1][1][2],indent);
	fprintf(file, "%f %f %f\n\n",crds[0][1][0][0],
		crds[0][1][0][1],crds[0][1][0][2]);

	fprintf(file, "%f %f %f%s",crds[0][0][0][0],
		crds[0][0][0][1],crds[0][0][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[0][1][0][0],
		crds[0][1][0][1],crds[0][1][0][2]);
	fprintf(file, "%f %f %f%s",crds[1][0][0][0],
		crds[1][0][0][1],crds[1][0][0][2],indent);
	fprintf(file, "%f %f %f\n",crds[1][1][0][0],
		crds[1][1][0][1],crds[1][1][0][2]);
	fprintf(file, "%f %f %f%s",crds[0][0][1][0],
		crds[0][0][1][1],crds[0][0][1][2],indent);
	fprintf(file, "%f %f %f\n",crds[0][1][1][0],
		crds[0][1][1][1],crds[0][1][1][2]);
	fprintf(file, "%f %f %f%s",crds[1][0][1][0],
		crds[1][0][1][1],crds[1][0][1][2],indent);
	fprintf(file, "%f %f %f\n\n",crds[1][1][1][0],
		crds[1][1][1][1],crds[1][1][1][2]);

	fprintf(file,"1 0 0 1\n");	
	fprintf(file,"%s}\n",indent);
}


LOCAL	void gview_cell3d_points(
	FILE *file,
	double **coords,
	int num,
	const char *color)
{
	static const char *indent = "	";
	int i;
	
	fprintf(file,"appearance{linewidth 10}\n");
	fprintf(file,"%s{\n%s%sVECT\n",indent,indent,indent);

	/*number of polygones, number of points, number of colors*/
	fprintf(file,"%d %d %d\n\n",num,num,num);
	
	for (i = 0; i < num; ++i)
	    fprintf(file, "1 "); /*num of points for each polyline*/
	fprintf(file,"\n\n");

	for (i = 0; i < num; ++i)
	    fprintf(file,"1 ");
	fprintf(file,"\n\n"); /*num of colors for each polyline*/

	for (i = 0; i < num; ++i)
	    fprintf(file,"%f %f %f\n", coords[i][0],coords[i][1],coords[i][2]);
	fprintf(file,"\n");

	for (i = 0; i < num; ++i)
	    fprintf(file,"%s\n",color); /*color for each point*/
	 
	fprintf(file,"%s}\n",indent);
}

LOCAL void FT_GetTrisInCell(
	INTERFACE *grid_intfc,
	const int *icoords,
	const int *shift,
	TRI ***tris,
	int *num_tris)
{
	struct Table    *T = table_of_interface(grid_intfc);
	int i = icoords[0]-shift[0];
	int j = icoords[1]-shift[1];
	int k = icoords[2]-shift[2];

	*tris = T->tris[k][j][i];
	*num_tris = T->num_of_tris[k][j][i];
}	/* end FT_GetTrisInCell */

LOCAL	void gview_cell3d_tris(
	FILE *file,
	TRI **tris,
	int num_tris)
{
	int i,j;
	POINT *p;
	char *indent = "	";
	
	fprintf(file,"%s{\n%s%sOFF\n%s%s%6d %6d %6d\n",
			indent,indent,indent,
			indent,indent,3*num_tris,num_tris,0);
	for (i = 0; i < num_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tris[i])[j];
		fprintf(file, "%s%s%-9g %-9g %-9g\n",
			indent,indent,
			Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    }
	}
	for (i = 0; i < num_tris; ++i)
	{
	    fprintf(file,"%s%s%-4d %-4d %-4d %-4d\n",indent,indent,
		3,3*i,3*i+1,3*i+2);
	}
	fprintf(file,"%s}\n",indent);
}
