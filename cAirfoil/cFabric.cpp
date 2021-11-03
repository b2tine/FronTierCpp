#include "cFabric.h"


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


//TODO: ADD
//CFABRIC_CARTESIAN::setElasticStates() etc.
//


