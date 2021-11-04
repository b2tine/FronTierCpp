#include "cFabric.h"

static double (*getStateVel[MAXD])(Locstate) =
               {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateMom[MAXD])(Locstate) =
               {getStateXmom,getStateYmom,getStateZmom};


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

void CFABRIC_CARTESIAN::addFluxAlongGridLine(
	int idir,
	int *grid_icoords,
	double dt,
	SWEEP *m_vst,
    FSWEEP *m_flux)
{
	int i,l,n,index;
	SCHEME_PARAMS scheme_params;
	EOS_PARAMS	*eos;
	
    COMPONENT comp;
	int seg_min,seg_max;
	static int icoords[MAXD];
	int icoords_next[MAXD];
	EQN_PARAMS *eqn_params = (EQN_PARAMS*)(front->extra1);
	INTERFACE *grid_intfc=front->grid_intfc;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	SURFACE **s;
	STATE *state;
	
    GRID_DIRECTION	ldir[3]={WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3]={EAST,NORTH,UPPER};

//	printf("icoords={%d,%d,%d}\n",icoords[0],icoords[1],icoords[2]);
//	printf("icoords=%d,%d,%d, wave_type=%d\n",grid_icoords[0],grid_icoords[1],grid_icoords[2],wave_type(hs));
	
    SWEEP vst;
    FSWEEP vflux;
    allocDirVstFlux(&vst,&vflux);
	/*if (first)
    {
        //TODO: more efficient to free and release?
        //      freeing and reallocating could keep from fragmenting/page faulting
        //      and improve efficiency...
        first = NO;
        allocDirVstFlux(&vst,&vflux);
    }*/

	scheme_params.lambda = dt/top_h[idir];
    scheme_params.beta = 0.0;
	scheme_params.artificial_compression = eqn_params->articomp;

    for (i = 0; i < dim; ++i)
	    icoords[i] = grid_icoords[i];


    //TODO: use these below if checking 2 crossings..
    //      if fabric not seperated well, will be able to tell from 
    //      distance between the crossings
	double ldir_crx_coords[MAXD];
	double rdir_crx_coords[MAXD];
    
    //TODO: See the TODO needBufferFromIntfc() regarding gas_comp().
    //      Need to make sure both work correctly with index coating
    //      algorithm used for ELASTIC_BOUNDARYs
    
    seg_min = imin[idir];	
	while (seg_min <= imax[idir])
	{
	    for (; seg_min <= imax[idir]; ++seg_min)
	    {
		    icoords[idir] = seg_min;
	    	index = d_index(icoords,top_gmax,dim);
	    	comp = top_comp[index];
	    	if (gas_comp(comp)) break;
	    }

	    if (seg_min > imax[idir]) break;
	    
        //TODO: is this good enough zeroing?
        //      what about the +7 in size value in allocDirVstFlux()??
        //      For now should be safe since allocating and freeing everytime.
        for (i = 0; i <= top_gmax[idir]; ++i)
	    {
	    	vst.dens[i] = 0.0; 
	    	vst.pres[i] = 0.0; 
	    	vst.engy[i] = 0.0; 
	    	vst.momn[0][i] = vst.momn[1][i] = vst.momn[2][i] = 0.0;
	    }
	    
        i = seg_min;
	    icoords[idir] = i;
	    index = d_index(icoords,top_gmax,dim);
	    comp = top_comp[index];
	    n = 0;

        vst.dens[n+nrad] = m_vst->dens[index];
        vst.engy[n+nrad] = m_vst->engy[index];
        vst.pres[n+nrad] = m_vst->pres[index];
	    
        //index cycling so 0-th component aligned along the idir-th gridline,
        //and remaining chosen to form right hand coordinate system
        for (l = 0; l < dim; ++l)
            vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
	    for (l = dim; l < 3; ++l)
            vst.momn[l][n+nrad] = 0.0;
	    
        seg_max = i;
	    n++;

//	    printf("Component=%d\n",comp);

	    for (i = seg_min + 1; i <= imax[idir]; i++)
	    {
            icoords[idir] = i;
            index = d_index(icoords,top_gmax,dim);
	
            for (int ii = 0; ii < dim; ++ii)
                icoords_next[ii] = icoords[ii];
            icoords_next[idir]++;
            
            boolean status1, status2;
            
            //TODO: use ldir_crx_coords and rdir_crx_coords to differentiate crossings
            status1 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords,rdir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);

            //TODO: status2 never gets checked...
            //      Guessing is for case that fabric is folded, but the bugs never
            //      quite got worked out so it was disabled.
            status2 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords_next,ldir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);

            //TODO: why did cxxu abandon/comment out these 2 blocks??
            //      Dead code or unfinished code?

//		For the following part, if needBufferFromIntfc is true, which means
// 		it meets a boundary, we have a break. If needBufferFromIntfc is not true,
// 		we check the two statuses. If one of them is true, we still have a break.
// 		If all three of them are false, we assume that it does not meet a boundary
   
            //if (needBufferFromIntfc(comp,top_comp[index]))
            //{
            //    printf("get boundary \n");
            //    break;
		    //}
//		if (status1){
//		    
//		    printf("Outside::: The wave_type is %d\n", wave_type(hs));
//		    if (!needBufferFromIntfc(comp,top_comp[index])){
//			seg_max=i;
//			n++;
//		    }
//		    break;
//		}
//
            if (needBufferFromIntfc(comp,top_comp[index]))
            {
                printf("get boundary \n");
                break;
            }
		    else
            {
                vst.dens[n+nrad] = m_vst->dens[index];
                vst.engy[n+nrad] = m_vst->engy[index];
                vst.pres[n+nrad] = m_vst->pres[index];

                for (l = 0; l < dim; ++l)
                    vst.momn[l][n+nrad] = m_vst->momn[(l+idir)%dim][index];
                for (l = dim; l < 3; ++l)
                    vst.momn[l][n+nrad] = 0.0;

                n++;
                if (status1)
                {

    //			printf("icoords[0]=%d icoords[1]=%d, icoords[2]=%d, wave_type=%d\n",icoords[0],icoords[1], icoords[2],wave_type(hs));

   //			if (wave_type(hs)==7) //NOTE: 7 is NUEMANN_BOUNDARY
   //  			{
   //                 for (int ii=0;ii<dim;ii++)
   //                 {
   //                    printf("icoords[%d]=%d,",ii,icoords[ii]);
   //                 }
   //                 printf("    found\n");
   //                 
   //                 printf("crx_coords[0]=%f\n",crx_coords[0]);
   //                 for (int ii=0;ii<dim;ii++)
   //                 {
   //                     printf("crx_coords[%d]=%f\n",ii,crx_coords[ii]);
   //                 }

   //             }
   // 
                    
                    seg_max = i++;
                    break;
                 }
            }
    
            seg_max = i;
        }

        //Elastic Boundary and Porosity accounted for in appendGhostBuffer()
	    icoords[idir] = seg_min;
	    appendGhostBuffer(&vst,m_vst,0,icoords,idir,0);
	    icoords[idir] = seg_max;
	    appendGhostBuffer(&vst,m_vst,n,icoords,idir,1);
	    
	    eos = &(eqn_params->eos[comp]);
	    EosSetTVDParams(&scheme_params, eos);
	    numericalFlux((POINTER)&scheme_params,&vst,&vflux,n);
		    
	    n = 0;
	    for (i = seg_min; i <= seg_max; ++i)
	    {
            icoords[idir] = i;
            index = d_index(icoords,top_gmax,dim);

            m_flux->dens_flux[index] += vflux.dens_flux[n+nrad];
            m_flux->engy_flux[index] += vflux.engy_flux[n+nrad];

            for (l = 0; l < dim; ++l)
            {
                m_flux->momn_flux[(l+idir)%dim][index] +=
                    vflux.momn_flux[l][n+nrad];
            }
            for (l = dim; l < 3; ++l)
                m_flux->momn_flux[l][index] = 0.0;
            
            n++;
	    }

	    seg_min = seg_max + 1;
	}
        
    freeDirVstFlux(&vst,&vflux);
}	/* end addFluxAlongGridLine */


void CFABRIC_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
	int		i,j,k,index,ic[MAXD];
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};
	HYPER_SURF 	*hs;
	HYPER_SURF	*hs1;
	COMPONENT 	comp;
	double 		crx_coords[MAXD];
	STATE 		*state,ghost_st;
	int		ind2[2][2] = {{0,1},{1,0}};
	int		ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int 		ic_next[MAXD];
	INTERFACE	*grid_intfc = front->grid_intfc;
	static int count = 0;
	count++;
	boolean Debug = NO;

	SURFACE **s;


	if (debugging("append_buffer"))
		printf("Entering appendGhostBuffer()\n");

	for (i = 0; i < dim; ++i)
        ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	comp = cell_center[index].comp;

	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
	    {
		ic[idir] = icoords[idir] - i;
		index = d_index(ic,top_gmax,dim);

//	 	The following is for debugging		    
		boolean status;
//		check neighbor in ldir[idir] 
		for (k = 0; k < dim; ++k)
		    ic_next[k] = ic[k];
		ic_next[idir]++;
		status = FT_StateStructAtGridCrossing(front,grid_intfc,
			ic_next,ldir[idir],comp,(POINTER*)&state,
			&hs,crx_coords);
/*
		if (status)
		if (status && wave_type(hs) != 6)
		    printf("HI,status=%d, wave_type(*hs)=%d\n",status,wave_type(hs));
*/

		if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status)
		{
		    vst->dens[nrad-i] = m_vst->dens[index];
		    vst->engy[nrad-i] = m_vst->engy[index];
		    vst->pres[nrad-i] = m_vst->pres[index];

		    for (j = 0; j < 3; j++)
			vst->momn[j][nrad-i] = 0.0;
		    if (dim == 1)
			vst->momn[0][nrad-i] = m_vst->momn[0][index];
		    else if (dim == 2)
                for(j = 0; j < 2; j++)
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3){
			for (j = 0; j < 3; j++){
			    vst->momn[j][nrad-i] = 
			    	 	m_vst->momn[ind3[idir][j]][index];
			}

		    }

		}
		else
		{
		    if (!status) /* extreme cases */
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			/* check neighbor in the opposite direction */
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]--;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,rdir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    /* check second neighbor in the same direction */
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] += 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,ldir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				/* must be something wrong */
		    	    	printf("In appendGhostBuffer() Case 0\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				/* check if crossing is close enough */
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k]+ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    /* check if crossing is close enough */
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                {
                    coords[k] = top_L[k] + ic[k] * top_h[k];
                    wtol[k] = crx_coords[k] - coords[k];
                    if (fabs(wtol[k]) > tol[k])
                        close_enough = NO;
                }
                if (!close_enough)
                {
                    (void) printf("ERROR: Not close enough!\n");
                    clean_up(ERROR);
			    }
            }
		    
            }

		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic_next,idir,
						nb,0,i,comp);
			break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,
					idir,nb,0,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	GFMGhostState(ic,comp,&ghost_st);
		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[nrad-k] = ghost_st.dens;
		    	    vst->engy[nrad-k] = ghost_st.engy;
		    	    vst->pres[nrad-k] = ghost_st.pres;
			
			    for (j=0; j < 3; j++)
			    	    vst->momn[j][nrad-k] = 0.0;
			    if (dim == 1)
				vst->momn[0][nrad-k] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for (j=0; j < 2; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for (j = 0; j < 3; j++)
				    vst->momn[j][nrad-k] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
					wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=", icoords,3,"\n");
		    	clean_up(ERROR);
		    }
		    break;
		}
	    
        }
	    break;
	case 1:
	    for (i = 1; i <= nrad; ++i)  //nrad=3
	    {
		ic[idir] = icoords[idir] + i;
		index = d_index(ic,top_gmax,dim);

//		For debugging
		boolean status;
//		check neighbor in rdir[idir] 
		for (k = 0; k < dim; ++k)
			ic_next[k] = ic[k];
		ic_next[idir]--;
		status = FT_StateStructAtGridCrossing(front,grid_intfc,
				icoords,rdir[idir],comp,(POINTER*)&state,
				&hs,crx_coords);

//For the needBufferFromIntfc function, if the two component are different,
//YES is returned. Then for the following, the if statement is satisfied when 
//the two component are the same, which means it does not meet the rectangle
//boundary. It may meet the elastic boundary or does not meet any boundary.
//Then !status exclude the possibility of meeting elastic boundary.


		if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status )
		{
//		    if (status && wave_type(hs) == 13)
//			printf("233: target boundary found.\n");
		    vst->dens[n+nrad+i-1] = m_vst->dens[index];
		    vst->engy[n+nrad+i-1] = m_vst->engy[index];
		    vst->pres[n+nrad+i-1] = m_vst->pres[index];
		    
		    for (j = 0; j < 3; j++)
			vst->momn[j][n+nrad+i-1] = 0.0;
		    if (dim == 1)
			vst->momn[0][n+nrad+i-1] = 
			         	m_vst->momn[0][index];
		    else if (dim == 2)
			for(j = 0; j < 2; j++)
			    	vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind2[idir][j]][index];
		    else if (dim == 3)
			for (j = 0; j < 3; j++)
			    vst->momn[j][n+nrad+i-1] = 
			         	m_vst->momn[ind3[idir][j]][index];
		}
		else
		{


		    if (!status) /* extreme cases */
		    {
			double coords[MAXD], wtol[MAXD], tol[MAXD];
			int ic_tmp[MAXD];
			for (k = 0; k < dim; ++k)
                            tol[k] = 2.0 * IG_TOL * top_h[k];
			/* check neighbor in the opposite direction */
			for (k = 0; k < dim; ++k)
                            ic_tmp[k] = ic[k];
                        ic_tmp[idir]++;
			status = FT_StateStructAtGridCrossing(front,grid_intfc,
					ic_tmp,ldir[idir],comp,(POINTER*)&state,
					&hs,crx_coords);
			if(!status)
			{
			    /* check second neighbor in the same direction */
			    for (k = 0; k < dim; ++k)
                                ic_tmp[k] = ic[k];
                            ic_tmp[idir] -= 2;
                            status = FT_StateStructAtGridCrossing(front,
                                        grid_intfc,ic_tmp,rdir[idir],comp,
                                        (POINTER*)&state,&hs,crx_coords);
			    if (!status)
			    {
				/* must be something wrong */
		    	    	printf("In appendGhostBuffer() Case 1\n");
		    	    	printf("ERROR: No crossing found!\n");
		    	    	print_int_vector("ic=",ic,dim,"\n");
		    	    	printf("direction: %s side %d\n",
		           		grid_direction_name(ldir[idir]), nb);
				clean_up(ERROR);
			    }
			    else
			    {
				/* check if crossing is close enough */
			        boolean close_enough = YES;
                                for (k = 0; k < dim; ++k)
                                {
                                    coords[k] = top_L[k] + ic_next[k]*top_h[k];
                                    wtol[k] = crx_coords[k] - coords[k];
                                    if (fabs(wtol[k]) > tol[k])
                                        close_enough = NO;
                                }
                                if (!close_enough)
                                {
                                    (void) printf("ERROR: Not close enough!\n");
                                    clean_up(ERROR);
                                }
			    }
			}
			else
			{
			    /* check if crossing is close enough */
			    boolean close_enough = YES;
			    for (k = 0; k < dim; ++k)
                            {
                                coords[k] = top_L[k] + ic[k] * top_h[k];
                                wtol[k] = crx_coords[k] - coords[k];
			        if (fabs(wtol[k]) > tol[k])
				    close_enough = NO;
			    }
			    if (!close_enough)
			    {
			        (void) printf("ERROR: Not close enough!\n");
			        clean_up(ERROR);
			    }
			}
		    }

		    switch (wave_type(hs))
		    {
		    case NEUMANN_BOUNDARY:
		    case MOVABLE_BODY_BOUNDARY:
		    	setNeumannStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case ELASTIC_BOUNDARY:
		    	setElasticStates(vst,m_vst,hs,state,ic_next,idir,
						nb,n,i,comp);
		    	break;
		    case DIRICHLET_BOUNDARY:
		    	setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,
						n,i);
		    	break;
		    case FIRST_PHYSICS_WAVE_TYPE:
		    	GFMGhostState(ic,comp,&ghost_st);

		    	for (k = i; k <= nrad; ++k)
		    	{
		    	    vst->dens[n+nrad+k-1] = ghost_st.dens;
		    	    vst->engy[n+nrad+k-1] = ghost_st.engy;
		    	    vst->pres[n+nrad+k-1] = ghost_st.pres;
			
			    for(j=0; j<3; j++)
			    	vst->momn[j][n+nrad+k-1] = 0.0;
			    if (dim == 1)
				vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
			    else if (dim == 2)
			    	for(j = 0; j < 2; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind2[idir][j]];
			    else if (dim == 3)
			    	for(j = 0; j < 3; j++)
				    vst->momn[j][n+nrad+k-1] = 
				     	    ghost_st.momn[ind3[idir][j]];
		    	}
		    	break;
		    default:
		    	(void) printf("In appendGhostBuffer(): ");
		    	(void) print_wave_type("Unknown wave type ",
				wave_type(hs),"\n",front->interf);
		    	(void) print_int_vector("icoords=",icoords,3,"\n");
		    	(void) printf("nb = %d\n",nb);
		    	clean_up(ERROR);
		    }
		    break;
		}
	    }
	}
}	/* end appendGhostBuffer */

void CFABRIC_CARTESIAN::setElasticStates(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
    switch (eqn_params->poro_scheme)
    {
        case PORO_SCHEME::NORMAL_REFLECTION:
            setElasticStatesRFB_normal(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
            break;
        case PORO_SCHEME::REFLECTION:
            setElasticStatesRFB(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
            break;
        case PORO_SCHEME::RIEMANN:
            setElasticStatesRiem(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
            break;
        default:
            printf("\n\nERROR setElasticStates(): unrecognized PORO_SCHEME\n\n");
            LOC(); clean_up(EXIT_FAILURE);
    }
}

//Reflection Boundary Formulation of Porosity -- No relative tangential velocity
void CFABRIC_CARTESIAN::setElasticStatesRFB_normal(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ind2[2][2] = {{0,1},{1,0}};
    int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
	double* vel_intfc = state->vel;
	double poro = eqn_params->porosity;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE st_tmp_real;
	STATE st_tmp_ghost;	

	st_tmp_real.dim = dim;
	st_tmp_real.eos = &eqn_params->eos[comp];

    st_tmp_ghost.dim = dim;
    st_tmp_ghost.eos = &eqn_params->eos[comp];
	//st_tmp_ghost.eos = state->eos;

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRFB_normal():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
	}

	    //if nb = 0, the point is above the boundary, and we
        //           select three points below the boundary
	    //if nb = 1, the point is below the boundary, and we
        //           select three points above the boundary

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

	    index_ghost = d_index(ic_ghost,top_gmax,dim);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
            coords_ref[j] = coords_ghost[j];
	    }
        
        /* Reflect ghost point through intfc-mirror at crossing */
        //first reflect across the grid line containing the intfc crossing 
	    double vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

	    for (j = 0; j < dim; ++j)
	    {
		    v[j] = coords_ref[j] - crx_coords[j];
		    vn += v[j]*nor[j];
	    }
           
        //reflect v across the line containing the normal vector
	    for (j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];
	    
        //desired reflected point
        for (j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];
			
        /* Interpolate the state at the reflected point */
	    
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->dens,getStateDens,&st_tmp_ghost.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->pres,getStatePres,&st_tmp_ghost.pres,&m_vst->pres[index]);
	    
        for (j = 0; j < dim; ++j)
        {
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->momn[j],
                    getStateMom[j],&st_tmp_ghost.momn[j],&m_vst->momn[j][index]);
        }
        
        //Compute relative normal velocity in frame of interface crossing.
        vn = 0.0;
        double v_reflect[3], v_rel[3];
        for (j = 0; j < dim; j++)
	    {
            v_reflect[j] = st_tmp_ghost.momn[j]/st_tmp_ghost.dens;
            vn += (v_reflect[j] - vel_intfc[j])*nor[j];
	    }
	    
        //Ghost vel has relative normal velocity component equal in magnitude to
        //reflected point's relative normal velocity and going in the opposite direction.
        //TODO: We leave the question of relative tangential velocity wrt to the intfc
        //      to be determined
        //
        //      NEEDS TO BE TESTED
        for (j = 0; j < dim; j++)
        {
            /*
            //allow relative tangential velocity
            v_ghost[j] = v_reflect[j] - 2.0*vn*nor[j];
            */

            //zero relative tangential velocity
            v_ghost[j] = vel_intfc[j] - vn*nor[j];
        }

	    st_tmp_real.dens = m_vst->dens[index_ghost];
	    st_tmp_real.pres = m_vst->pres[index_ghost];
	    
        st_tmp_ghost.dens = (1.0 - poro)*st_tmp_ghost.dens + poro*st_tmp_real.dens;
	    st_tmp_ghost.pres = (1.0 - poro)*st_tmp_ghost.pres + poro*st_tmp_real.pres;
	  
        for (j = 0; j < dim; ++j)
        {
            st_tmp_real.momn[j] = m_vst->momn[j][index_ghost];
            v_real[j] = st_tmp_real.momn[j]/st_tmp_real.dens;
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            st_tmp_ghost.momn[j] = v_ghost[j]*st_tmp_ghost.dens;
        }
	    
	    st_tmp_ghost.engy = EosEnergy(&st_tmp_ghost);

	    /* debugging printout */
	    if (st_tmp_ghost.engy < 0.0 || st_tmp_ghost.eos->gamma < 0.001)
	    {
            printf("negative engrgy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",st_tmp_ghost.dens,st_tmp_ghost.momn[0],
                st_tmp_ghost.momn[1],st_tmp_ghost.momn[2],st_tmp_ghost.pres,
                st_tmp_ghost.engy);
            printf("st_tmp_ghost.dim = %d, idir = %d, nb = %d \n",
                st_tmp_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp_ghost.eos->gamma,
                st_tmp_ghost.eos->einf,st_tmp_ghost.eos->pinf);
            printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
                            coords_ref[2]);
            clean_up(EXIT_FAILURE);
	    }

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = st_tmp_ghost.dens;
            vst->engy[nrad-i] = st_tmp_ghost.engy;
            vst->pres[nrad-i] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = st_tmp_ghost.momn[0];
            }
	    }
	    else
	    {
            /* Debug selectively!
            if (debugging("crx_reflection"))
            {
                    sprintf(fname,"intfc-%d-%d",count,i);
                    sprintf(fname,"intfc-xx");
                    xgraph_2d_reflection(fname,front->grid_intfc,coords,
                    crx_coords,coords_ref,nor);
            }
            */
            vst->dens[n+nrad+i-1] = st_tmp_ghost.dens;
            vst->engy[n+nrad+i-1] = st_tmp_ghost.engy;
            vst->pres[n+nrad+i-1] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = st_tmp_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRFB_normal()\n");
}	/* end setElasticStatesRFB_normal */

//Reflection Boundary Formulation of Porosity -- Allow relative tangential velocity
void CFABRIC_CARTESIAN::setElasticStatesRFB(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ind2[2][2] = {{0,1},{1,0}};
    int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
	double* vel_intfc = state->vel;
	double poro = eqn_params->porosity;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE st_tmp_real;
	STATE st_tmp_ghost;	

	st_tmp_real.dim = dim;
	st_tmp_real.eos = &eqn_params->eos[comp];

    st_tmp_ghost.dim = dim;
    st_tmp_ghost.eos = &eqn_params->eos[comp];
	//st_tmp_ghost.eos = state->eos;

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
	
    dir = (nb == 0) ? ldir[idir] : rdir[idir];
	FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);

	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRFB():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	    (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
	}

	    //if nb = 0, the point is above the boundary, and we
        //           select three points below the boundary
	    //if nb = 1, the point is below the boundary, and we
        //           select three points above the boundary

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

	    index_ghost = d_index(ic_ghost,top_gmax,dim);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
            coords_ref[j] = coords_ghost[j];
	    }
        
        /* Reflect ghost point through intfc-mirror at crossing */
        //first reflect across the grid line containing the intfc crossing 
	    double vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

	    for (j = 0; j < dim; ++j)
	    {
		    v[j] = coords_ref[j] - crx_coords[j];
		    vn += v[j]*nor[j];
	    }
           
        //reflect v across the line containing the normal vector
	    for (j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];
	    
        //desired reflected point
        for (j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];
			
        /* Interpolate the state at the reflected point */
	    
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->dens,getStateDens,&st_tmp_ghost.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->pres,getStatePres,&st_tmp_ghost.pres,&m_vst->pres[index]);
	    
        for (j = 0; j < dim; ++j)
        {
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->momn[j],
                    getStateMom[j],&st_tmp_ghost.momn[j],&m_vst->momn[j][index]);
        }
        
        vn = 0.0;
        double v_reflect[3], v_rel[3];
        for (j = 0; j < dim; j++)
	    {
            v_reflect[j] = st_tmp_ghost.momn[j]/st_tmp_ghost.dens;
            vn += (v_reflect[j] - vel_intfc[j])*nor[j];
	    }
	    
        //Ghost vel has relative normal velocity component equal in magnitude to
        //reflected point's relative normal velocity and going in the opposite direction.
        //TODO: We leave the question of relative tangential velocity wrt to the intfc
        //      to be determined
        //
        //      NEEDS TO BE TESTED
        for (j = 0; j < dim; j++)
        {
            //allow relative tangential velocity
            v_ghost[j] = v_reflect[j] - 2.0*vn*nor[j];
            
            /*
            //zero relative tangential velocity
            v_ghost[j] = vel_intfc[j] - vn*nor[j];
            */
        }

	    st_tmp_real.dens = m_vst->dens[index_ghost];
	    st_tmp_real.pres = m_vst->pres[index_ghost];
	    
        st_tmp_ghost.dens = (1.0 - poro)*st_tmp_ghost.dens + poro*st_tmp_real.dens;
	    st_tmp_ghost.pres = (1.0 - poro)*st_tmp_ghost.pres + poro*st_tmp_real.pres;
	  
        for (j = 0; j < dim; ++j)
        {
            st_tmp_real.momn[j] = m_vst->momn[j][index_ghost];
            v_real[j] = st_tmp_real.momn[j]/st_tmp_real.dens;
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            st_tmp_ghost.momn[j] = v_ghost[j]*st_tmp_ghost.dens;
        }
	    
	    st_tmp_ghost.engy = EosEnergy(&st_tmp_ghost);

	    /* debugging printout */
	    if (st_tmp_ghost.engy < 0.0 || st_tmp_ghost.eos->gamma < 0.001)
	    {
            printf("negative engrgy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",st_tmp_ghost.dens,st_tmp_ghost.momn[0],
                st_tmp_ghost.momn[1],st_tmp_ghost.momn[2],st_tmp_ghost.pres,
                st_tmp_ghost.engy);
            printf("st_tmp_ghost.dim = %d, idir = %d, nb = %d \n",
                st_tmp_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp_ghost.eos->gamma,
                st_tmp_ghost.eos->einf,st_tmp_ghost.eos->pinf);
            printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
                            coords_ref[2]);
            clean_up(EXIT_FAILURE);
	    }

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = st_tmp_ghost.dens;
            vst->engy[nrad-i] = st_tmp_ghost.engy;
            vst->pres[nrad-i] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = st_tmp_ghost.momn[0];
            }
	    }
	    else
	    {
            /* Debug selectively!
            if (debugging("crx_reflection"))
            {
                    sprintf(fname,"intfc-%d-%d",count,i);
                    sprintf(fname,"intfc-xx");
                    xgraph_2d_reflection(fname,front->grid_intfc,coords,
                    crx_coords,coords_ref,nor);
            }
            */
            vst->dens[n+nrad+i-1] = st_tmp_ghost.dens;
            vst->engy[n+nrad+i-1] = st_tmp_ghost.engy;
            vst->pres[n+nrad+i-1] = st_tmp_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = st_tmp_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRFB()\n");
}	/* end setElasticStatesRFB */

//Riemann Problem Formulation of Porosity
void CFABRIC_CARTESIAN::setElasticStatesRiem(
	SWEEP		*vst,
	SWEEP		*m_vst,
	HYPER_SURF 	*hs,
	STATE		*state,
	int		*icoords,
	int		idir,
	int		nb,
	int		n,
	int		istart,
	COMPONENT	comp)
{
	int i,j,index,index_ghost;
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD],crx_coords[MAXD];
	double	nor[MAXD],v[MAXD],v_ghost[MAXD],v_real[MAXD];
	
    double vn, vn_intfc;
	
	GRID_DIRECTION  dir;
	GRID_DIRECTION 	ldir[3] = {WEST,SOUTH,LOWER};
	GRID_DIRECTION 	rdir[3] = {EAST,NORTH,UPPER};

	STATE sl, sr, state_ghost;

	sl.dim = sr.dim = state_ghost.dim = dim;
	sl.eos = sr.eos = state_ghost.eos = &eqn_params->eos[comp];

	index = d_index(icoords,top_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    ic_ghost[i] = icoords[i];
	}
    dir = (nb == 0) ? ldir[idir] : rdir[idir];

    /*
	if (debugging("elastic_buffer"))
	{
	    (void) printf("\nEntered setElasticStatesRiem():\n");
	    (void) printf("comp = %d\n",comp);
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],
				icoords[2]);
	    (void) printf("idir = %d nb = %d\n",idir,nb);
	    (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
	    (void) print_general_vector("coords = ",coords,dim,"\n");
	    (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
	    (void) print_general_vector("nor = ",nor,dim,"\n");
	}
    */

	for (i = istart; i <= nrad; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

        //ghost point coords
	    for (j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
	    }
        
        index_ghost = d_index(ic_ghost,top_gmax,dim);
	    COMPONENT comp_ghost = cell_center[index_ghost].comp;

        //Find the closest interface point
        boolean status;
        double intrp_a[3];
        HYPER_SURF_ELEMENT* nearHse;
        HYPER_SURF* nearHs;

        status = FT_FindNearestIntfcPointInRange(front,comp_ghost,
                coords_ghost,INCLUDE_BOUNDARIES,crx_coords,intrp_a,
                &nearHse,&nearHs,5);
        
        if (!status) 
        {
            LOC();
            printf("ERROR: could not find interface point\n");
            clean_up(EXIT_FAILURE);
        }
        
        //Get 2 points straddling interface in normal direction
        double pl[MAXD], pr[MAXD], nor[MAXD];

        //TODO: Use all three points of triangle in approximation,
        //      and use the triangle normal like in compute_total_canopy_force3d()
        TRI* nearTri = Tri_of_hse(nearHse);
        FT_NormalAtPoint(Point_of_tri(nearTri)[0],front,nor,comp);
        STATE* state_intfc = (STATE*)left_state(Point_of_tri(nearTri)[0]);
            //nor = Tri_normal_vector(nearTri);//USE THIS
        double h = FT_GridSizeInDir(nor,front);

        for (j = 0; j < 3; ++j)
        {
            pl[j] = crx_coords[j] + 1.5*h*nor[j];
            pr[j] = crx_coords[j] - 1.5*h*nor[j];//on ghost point side
        }

        //Interpolate states for the 2 points
        FT_IntrpStateVarAtCoords(front,comp,pl,
                m_vst->dens,getStateDens,&sl.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,pl,
                m_vst->pres,getStatePres,&sl.pres,&m_vst->pres[index]);
	    
        FT_IntrpStateVarAtCoords(front,comp_ghost,pr,
                m_vst->dens,getStateDens,&sr.dens,&m_vst->dens[index_ghost]);
	    FT_IntrpStateVarAtCoords(front,comp_ghost,pr,
                m_vst->pres,getStatePres,&sr.pres,&m_vst->pres[index_ghost]);
        
        //TODO: retain the tangential velocities of the interpolated states
        //      and add to the normal ghost velocity obtained from the
        //      riemann problem.
        //
        //Using relative velocity wrt to interface velocity
        //TODO: can we do better than this interpolation?
        /*
        double intfc_dens;
        double intfc_momn[3], vel_intfc[3];
        FT_IntrpStateVarAtCoords(front,comp,crx_coords,
                m_vst->dens,getStateDens,&intfc_dens,&m_vst->dens[index]);
        */

        //TODO: Should the left state velocity be set to the right state velocity
        //      as in GFM for strong shock interaction part 1&2 (or is that a typo)
        //      Don't think it is a typo as results are much better.
        //      Need to understand why.
        double vl[3], vr[3];
        for (j = 0; j < dim; ++j)
        {
            /*
            FT_IntrpStateVarAtCoords(front,comp,crx_coords,m_vst->momn[j],
                    getStateMom[j],&intfc_momn[j],&m_vst->momn[j][index]);
            vel_intfc[j] = intfc_momn[j]/intfc_dens;
            */
            
            FT_IntrpStateVarAtCoords(front,comp,pl,m_vst->momn[j],
                    getStateMom[j],&sl.momn[j],&m_vst->momn[j][index]);
            vl[j] = sl.momn[j]/sl.dens - state_intfc->vel[j];//relative left state vel
            //vl[j] = sl.momn[j]/sl.dens - vel_intfc[j];//relative left state vel
                //vl[j] = sl.momn[j]/sl.dens;

            FT_IntrpStateVarAtCoords(front,comp_ghost,pr,m_vst->momn[j],
                    getStateMom[j],&sr.momn[j],&m_vst->momn[j][index_ghost]);
            vr[j] = sr.momn[j]/sr.dens - state_intfc->vel[j];//relative right state vel
            //vr[j] = sr.momn[j]/sr.dens - vel_intfc[j];//relative right state vel
                //vr[j] = sr.momn[j]/sr.dens;
        }

        double nor_vl = 0.0;
        double nor_vr = 0.0;
        for (j = 0; j < 3; ++j)
        {
            //relative normal velocities
            nor_vl += vl[j]*nor[j];
            nor_vr += vr[j]*nor[j];
        }

        /*
        //DEBUG: assuming centerline is through (0.5, 0.5, z)
        double ctrlinedist = sqrt(sqr(coords[0] - 0.5) + sqr(coords[1] - 0.5));
        double tolprint = sqrt(sqr(top_h[0]) + sqr(top_h[1]));
        if (ctrlinedist < tolprint)
        {
            printf("comp = %d\n",comp);
            printf("comp_ghost = %d\n",comp_ghost);
            print_general_vector("coords = ",coords,dim,"\n");
            print_general_vector("coords_ghost = ",coords_ghost,dim,"\n");
            print_general_vector("crx_coords = ",crx_coords,dim,"\n");
            print_general_vector("nor = ",nor,dim,"\n");
            print_general_vector("pr = ",pr,dim,"\n");
            print_general_vector("pl = ",pl,dim,"\n");

            printf("input states: sl sr\n");
            printf("\tdens: %f %f\n",sl.dens,sr.dens);
            printf("\tvn: %f %f\n",nor_vl,nor_vr);
            printf("\tpres: %f %f\n",sl.pres,sr.pres);
        }
        */

        //solve 1d riemann problem in interface normal direction
        RIEMANN_INPUT riem_input;
        RIEMANN_SOLN riem_soln;

        riem_input.left_state.d = sl.dens;
        riem_input.left_state.p = sl.pres;
        riem_input.left_state.u = nor_vr;//Typo or correct? UPDATE: not a typo
            //riem_input.left_state.u = nor_vl;
        riem_input.left_state.gamma = sl.eos->gamma;

        riem_input.right_state.d = sr.dens;
        riem_input.right_state.p = sr.pres;
        riem_input.right_state.u = nor_vr;
        riem_input.right_state.gamma = sr.eos->gamma;

        bool rp_status;
        rp_status = RiemannSolution(riem_input,&riem_soln);
        if (!rp_status)
        {
            printf("\nERROR: RiemannSolution()\n");
            printf("input states: sl sr\n");
            printf("\tpres: %f %f\n",sl.pres,sr.pres);
            printf("\tdens: %f %f\n",sl.dens,sr.dens);
            printf("\tvn: %f %f\n",nor_vl,nor_vr);
            clean_up(EXIT_FAILURE);
        }

        /*
        //DEBUG
        if (ctrlinedist < tolprint)
        {
            printf("\nsoln states, velocity\n");
            printf("right_state.u = %g\n",riem_soln.right_state.u);
            printf("right_center_state.u = %g\n",riem_soln.right_center_state.u);
            printf("left_center_state.u = %g\n",riem_soln.left_center_state.u);
            printf("left_state.u = %g\n",riem_soln.left_state.u);

            printf("\nsoln states, density\n");
            printf("right_state.d = %g\n",riem_soln.right_state.d);
            printf("right_center_state.d = %g\n",riem_soln.right_center_state.d);
            printf("left_center_state.d = %g\n",riem_soln.left_center_state.d);
            printf("left_state.d = %g\n",riem_soln.left_state.d);
           
            printf("\nsoln states, pressure\n");
            printf("right_state.p = %g\n",riem_soln.right_state.p);
            printf("right_center_state.p = %g\n",riem_soln.right_center_state.p);
            printf("left_center_state.p = %g\n",riem_soln.left_center_state.p);
            printf("left_state.p = %g\n",riem_soln.left_state.p);
        }
        */
       
        ///////////////////////////////////////////////////////////////////////
        /*
        ///////////////////////////////////////////////////////////////////////
        //TODO: This was first attempt that used the center state
        //          (centered on interface) soln.
        RIEM_STATE riem_soln_intfc;
        rp_status = RiemannSolnAtXi(&riem_soln,&riem_soln_intfc,0.0);
        if (!rp_status)
        {
            printf("ERROR: RiemannSolnAtXi()\n");
            clean_up(EXIT_FAILURE);
        }

        //printf("riem_soln_intfc:\n");
        //printf("\t(d,u,p) = %f %f %f\n",riem_soln_intfc.d,
          //      riem_soln_intfc.u,riem_soln_intfc.p);
        
        //Assign the midpoint solution state to the ghost state
        double dens_ghost = riem_soln_intfc.d;
        double pres_ghost = riem_soln_intfc.p;
        double vn_ghost = riem_soln_intfc.u;
        */

        /*
        RIEM_STATE left_center_state = riem_soln.left_center_state;
        double dens_ghost = left_center_state.d;
        double pres_ghost = left_center_state.p;
        double vn_ghost = left_center_state.u;
        */
        
        RIEM_STATE right_center_state = riem_soln.right_center_state;
        double dens_ghost = right_center_state.d;
        double pres_ghost = right_center_state.p;
        double vn_ghost = right_center_state.u;
        //TODO: Try using the interpolated left state pressure and density values,
        //      and only derive the velocity from the riemann solution.
            //double dens_ghost = sl.dens; //copy from the interpolated left state
            //double pres_ghost = sl.pres;  //copy from the interpolated left state

        //take weighted average using porosity to get the modified ghost point
	    double poro = eqn_params->porosity;
        state_ghost.dens = (1.0 - poro)*dens_ghost + poro*m_vst->dens[index_ghost];
        state_ghost.pres = (1.0 - poro)*pres_ghost + poro*m_vst->pres[index_ghost];
        
        double v_real[3];
        for (j = 0; j < dim; ++j)
        {
            v_ghost[j] = vn_ghost*nor[j] + state_intfc->vel[j];
                //v_ghost[j] = vn_ghost*nor[j] + vel_intfc[j];

            //reflect normal velocity component and convert back to world frame
                //v_ghost[j] = vel_intfc[j] - vn_ghost*nor[j];
                    //v_ghost[j] = vn_ghost*nor[j];
            v_real[j] = m_vst->momn[j][index_ghost]/m_vst->dens[index_ghost];
            v_ghost[j] = (1.0 - poro)*v_ghost[j] + poro*v_real[j];
            state_ghost.momn[j] = v_ghost[j]*state_ghost.dens;
        }
        state_ghost.engy = EosEnergy(&state_ghost);

	    // debugging printout
	    if (state_ghost.engy < 0.0 || state_ghost.eos->gamma < 0.001)
	    {
            printf("ERROR: Negative Energy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("%f %f %f %f %f %f \n",state_ghost.dens,state_ghost.momn[0],
                state_ghost.momn[1],state_ghost.momn[2],state_ghost.pres,
                state_ghost.engy);
            printf("state_ghost.dim = %d, idir = %d, nb = %d \n",
                state_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",state_ghost.eos->gamma,
                state_ghost.eos->einf,state_ghost.eos->pinf);
            printf("coords_ghost = %f %f %f \n",coords_ghost[0],coords_ghost[1],
                            coords_ghost[2]);
            clean_up(EXIT_FAILURE);
	    }

	    int ind2[2][2] = {{0,1},{1,0}};
        int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};

	    if (nb == 0)
	    {
            vst->dens[nrad-i] = state_ghost.dens;
            vst->engy[nrad-i] = state_ghost.engy;
            vst->pres[nrad-i] = state_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = state_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = state_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = state_ghost.momn[0];
            }
	    }
	    else
	    {
            vst->dens[n+nrad+i-1] = state_ghost.dens;
            vst->engy[n+nrad+i-1] = state_ghost.engy;
            vst->pres[n+nrad+i-1] = state_ghost.pres;
	    	for (j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
	    	if (dim == 3)
            {
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = state_ghost.momn[ind3[idir][j]];
            }
	    	else if (dim == 2)
            {
                for (j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = state_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = state_ghost.momn[0];
            }
	    }
	}

	if (debugging("elastic_buffer"))
        (void) printf("Leaving setElasticStatesRiem()\n");
}	/* end setElasticStatesRiem */


