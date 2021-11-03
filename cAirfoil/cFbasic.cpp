#include "cFluid.h"


void CFABRIC_CARTESIAN::applicationSetComponent()
{
	int i,icrd[MAXD],ic;
        int size = (int)cell_center.size();

        // cell center components
        for (i = 0; i < size; i++)
        {
            cell_center[i].comp =
                        getComponent(cell_center[i].icoords);
        }
	if (debugging("set_shifted_states"))
	{
	    printf("Sample component in applicationSetComponent()\n");
	    if (dim == 3)
	    {
            	icrd[0] = top_gmax[0]/2;
            	for (icrd[2] = 0; icrd[2] <= top_gmax[2]; ++icrd[2])
            	{
                    for (icrd[1] = 0; icrd[1] <= top_gmax[1]; ++icrd[1])
                    {
                        ic = d_index(icrd,top_gmax,dim);
                        printf("%d",top_comp[ic]);
                    }
                    printf("\n");
            	}
	    }
        }
}	/* end applicationSetComponent */

void CFABRIC_CARTESIAN::applicationSetStates()
{
	double coords[MAXD];
	int *icoords;
	int i,j,size = (int)cell_center.size();
	int id;
	STATE state;
	int ave_comp;
	double p_intfc[MAXD],t[MAXD];
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	double dist;
	double **vel = field.vel;
	double *pres = field.pres;
	
	setDomain();
	for (i = 0; i < size; i++)
        {
            icoords = cell_center[i].icoords;
            if (cell_center[i].comp != -1 &&
                cell_center[i].comp != top_comp[i])
            {
		for (j = 0; j < dim; ++j)
		    coords[j] = top_L[j] + icoords[j]*top_h[j];
		id = d_index(icoords,top_gmax,dim);
		if (fabs(cell_center[i].comp - top_comp[i]) != 2)
		    continue;

		if (debugging("set_crossed_state"))
		{
		    double r;
		    printf("\n");
		    printf("Shifted component:\n");
		    printf("icoords = %d %d %d\n",icoords[0],icoords[1],
					icoords[2]);
		    printf("old comp = %d  new comp = %d\n",
					cell_center[i].comp,top_comp[i]);
		    r = sqrt(sqr(coords[0] - 7.0) + sqr(coords[1] - 7.0));
		    printf("Radius = %f\n",r);
		}

		ave_comp = (cell_center[i].comp + top_comp[i])/2;
		if (!FT_FindNearestIntfcPointInRange(front,ave_comp,
			coords,NO_BOUNDARIES,p_intfc,t,&hse,&hs,2))
		    continue;

		dist = 0.0;
		for (j = 0; j < dim; ++j)
		    dist += sqr(coords[j] - p_intfc[j]);
		dist = sqrt(dist);
		if (debugging("set_crossed_state"))
		{
		    printf("coords  = %f %f %f\n",coords[0],coords[1],
					coords[2]);
		    printf("p_intfc = %f %f %f\n",p_intfc[0],p_intfc[1],
					p_intfc[2]);
		}
		if (dist > top_h[0]*Time_step_factor(front))
		{
		    if (debugging("set_crossed_state"))
			printf("external point: dist = %f\n",dist);
		    continue;
		}

		FrontNearestIntfcState(front,coords,ave_comp,(POINTER)&state);

		if (debugging("set_crossed_state"))
		{
		    printf("Old velocity  : %f %f %f\n",vel[0][id],
				vel[1][id],vel[2][id]);
		    printf("Intfc velocity: %f %f %f\n",state.vel[0],
			state.vel[1],state.vel[2]);
		    printf("Old pressure   = %f  \n",
				pres[id]);
		    printf("Intfc pressure = %f  \n",
				state.pres);
		}
		for (j = 0; j < dim; ++j)
		    vel[j][id] = state.vel[j];
	    }
        }
	FT_FreeGridIntfc(front);
	FT_MakeGridIntfc(front);
}	/* end applicationSetStates */

