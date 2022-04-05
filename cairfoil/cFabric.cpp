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
        cell_center[i].comp = getComponent(cell_center[i].icoords);
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
        if (cell_center[i].comp != -1 && cell_center[i].comp != top_comp[i])
        {
            for (j = 0; j < dim; ++j)
                coords[j] = top_L[j] + icoords[j]*top_h[j];
    
            id = d_index(icoords,top_gmax,dim);
            if (fabs(cell_center[i].comp - top_comp[i]) != 2) continue;

            if (debugging("set_crossed_state"))
            {
                double r;
                printf("\n");
                printf("Shifted component:\n");
                printf("icoords = %d %d %d\n",icoords[0],icoords[1],icoords[2]);
                printf("old comp = %d  new comp = %d\n",cell_center[i].comp,top_comp[i]);
                r = sqrt(sqr(coords[0] - 7.0) + sqr(coords[1] - 7.0));
                printf("Radius = %f\n",r);
            }

            ave_comp = (cell_center[i].comp + top_comp[i])/2;
            if (!FT_FindNearestIntfcPointInRange(front,ave_comp,
                        coords,NO_BOUNDARIES,p_intfc,t,&hse,&hs,2)) continue;

            dist = 0.0;
            for (j = 0; j < dim; ++j)
                dist += sqr(coords[j] - p_intfc[j]);
            dist = sqrt(dist);
            
            if (debugging("set_crossed_state"))
            {
                printf("coords  = %f %f %f\n",coords[0],coords[1],coords[2]);
                printf("p_intfc = %f %f %f\n",p_intfc[0],p_intfc[1],p_intfc[2]);
            }

            if (dist > hmin*Time_step_factor(front))
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
	int i,l,n,index,index_next;
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
    
	if (debugging("trace"))
		printf("Entering addFluxAlongGridLine()\n");

    
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

	    for (i = seg_min + 1; i <= imax[idir]; i++)
	    {
            icoords[idir] = i;
            index = d_index(icoords,top_gmax,dim);
	
            //switch from far component to near elastic intfc component
            if (std::abs(comp - top_comp[index]) == 1)
            {
                comp = top_comp[index];
            }

            for (int ii = 0; ii < dim; ++ii)
                icoords_next[ii] = icoords[ii];
            icoords_next[idir]++;
            
            index_next = d_index(icoords_next,top_gmax,dim);

            
            boolean status1;
            
            //TODO: use ldir_crx_coords and rdir_crx_coords to differentiate crossings
            status1 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords,rdir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);

            /*
            //if (status1 && wave_type(hs) != DIRICHLET_BOUNDARY)
            if (status1)
            {
                LOC();
                printf("In addFluxAlongGridLine() -- FT_StateStructAtGridCrossing()\n");
                printf("comp = %d\n",comp);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("positive_component(hs) = %d\n",positive_component(hs));
                printf("negative_component(hs) = %d\n",negative_component(hs));
                printf("icoords = %d %d %d\n",icoords[0],icoords[1],icoords[2]);
                printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
                printf("wave_type(hs) = %d\n\n",wave_type(hs));
            }
            else if (!status1)
            {
                LOC();
                printf("In addFluxAlongGridLine() -- !FT_StateStructAtGridCrossing()\n");
                printf("comp = %d\n",comp);
                printf("cell_center[index].comp = %d\n",cell_center[index].comp);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("top_comp[index_next] = %d\n",top_comp[index_next]);
                double check_coords[MAXD];
                for (int j = 0; j < dim; ++j)
                {
                    check_coords[j] = top_L[j] + icoords[j]*top_h[j];
                }
                printf("coords = %f %f %f\n\n",
                        check_coords[0], check_coords[1], check_coords[2]);
            }
            */

            /*
            //TODO: status2 never gets checked...
            //      Guessing is for case that fabric is folded, but the bugs never
            //      quite got worked out so it was disabled.
            
            boolean status2;
            status2 = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords_next,ldir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);
            */
            
            //if (comp != top_comp[index] || !gas_comp(top_comp[index]))
            //if (status1 || !gas_comp(top_comp[index]))
            //if (!gas_comp(top_comp[index]))
            if (needBufferFromIntfc(comp,top_comp[index]))
            {
                //printf("get boundary \n");
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
	    
	    eos = &(eqn_params->eos[GAS_COMP2]);
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

	if (debugging("trace"))
		printf("Leaving addFluxAlongGridLine()\n");

}	/* end addFluxAlongGridLine */


void CFABRIC_CARTESIAN::appendGhostBuffer(
	SWEEP *vst,
	SWEEP *m_vst,
	int n,
	int *icoords,
	int idir,
	int nb)
{
    EQN_PARAMS *eqn_params = (EQN_PARAMS*)(front->extra1);
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


	if (debugging("trace"))
		printf("Entering appendGhostBuffer()\n");

	for (i = 0; i < dim; ++i)
        ic[i] = icoords[i];
	
	index = d_index(ic,top_gmax,dim);
	//comp = cell_center[index].comp;
	comp = top_comp[index];

    bool nobreak = false;

	switch(nb)
	{
	case 0:
	    for (i = 1; i <= nrad; ++i)
        {
            ic[idir] = icoords[idir] - i;
            index = d_index(ic,top_gmax,dim);
            
            //check neighbor in ldir[idir] 
            for (k = 0; k < dim; ++k)
                ic_next[k] = ic[k];
            
            ic_next[idir]++;
            int index_next = d_index(ic_next,top_gmax,dim);
            
            boolean status;

            /*
            status = FT_StateStructAtGridCrossing(front,grid_intfc,
                icoords,ldir[idir],comp,(POINTER*)&state,
                &hs,crx_coords);
            */

            status = FT_StateStructAtGridCrossing(front,grid_intfc,
                ic_next,ldir[idir],comp,(POINTER*)&state,
                &hs,crx_coords);

            /*
            if (status && wave_type(hs) != DIRICHLET_BOUNDARY)
            {
                LOC();
                printf("In appendGhostBuffer() -- FT_StateStructAtGridCrossing()\n");
                printf("nb = %d\n",nb);
                printf("comp = %d n = %d\n",comp,n);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("top_comp[index_next] = %d\n",top_comp[index_next]);
                printf("positive_component(hs) = %d\n",positive_component(hs));
                printf("negative_component(hs) = %d\n",negative_component(hs));
                printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
                printf("wave_type(hs) = %d\n",wave_type(hs));
                printf("ic = %d %d %d\n",ic[0], ic[1], ic[2]);
                printf("icoords = %d %d %d\n",icoords[0], icoords[1], icoords[2]);
                printf("ic_next = %d %d %d\n\n",ic_next[0], ic_next[1], ic_next[2]);
            }
            else if (!status)
            {
                LOC();
                printf("In appendGhostBuffer() -- !FT_StateStructAtGridCrossing()\n");
                printf("nb = %d \t comp = %d\n",nb,comp);
                printf("cell_center[index].comp = %d\n",cell_center[index].comp);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("top_comp[index_next] = %d\n",top_comp[index_next]);
                printf("ic = %d %d %d\n",ic[0], ic[1], ic[2]);
                printf("icoords = %d %d %d\n",icoords[0], icoords[1], icoords[2]);
                printf("ic_next = %d %d %d\n\n",ic_next[0], ic_next[1], ic_next[2]);
            }
            */

            //if ((comp == top_comp[index] && gas_comp(cell_center[index].comp)) || !status)
            if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status)
            {
                vst->dens[nrad-i] = m_vst->dens[index];
                vst->engy[nrad-i] = m_vst->engy[index];
                vst->pres[nrad-i] = m_vst->pres[index];

                for (j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = 0.0;
                
                if (dim == 1)
                {
                    vst->momn[0][nrad-i] = m_vst->momn[0][index];
                }
                else if (dim == 2)
                {
                    for(j = 0; j < 2; j++)
                    {
                        vst->momn[j][nrad-i] = m_vst->momn[ind2[idir][j]][index];
                    }
                }
                else if (dim == 3)
                {
                    for (j = 0; j < 3; j++)
                    {
                        vst->momn[j][nrad-i] = m_vst->momn[ind3[idir][j]][index];
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
                        
                        status = FT_StateStructAtGridCrossing(front,grid_intfc,
                                ic_tmp,ldir[idir],comp,(POINTER*)&state,
                                &hs,crx_coords);
                        
                        if (!status)
                        {
                            /* must be something wrong */
                            printf("In appendGhostBuffer() Case 0\n");
                            printf("ERROR: No crossing found!\n");
                            print_int_vector("ic=",ic,dim,"\n");
                            print_int_vector("ic_next=",ic_next,dim,"\n");
                            printf("direction: %s side %d comp %d\n",
                            grid_direction_name(ldir[idir]), nb, comp);
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
                    setNeumannStates(vst,m_vst,hs,state,ic_next,idir,nb,0,i,comp);
                    break;

                case ELASTIC_BOUNDARY:
                    setElasticStates(vst,m_vst,hs,state,ic_next,idir,nb,0,i,comp);
                    break;
    
                case DIRICHLET_BOUNDARY:
                    setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,0,i,comp);
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
                        {
                            vst->momn[0][nrad-k] = ghost_st.momn[0];
                        }
                        else if (dim == 2)
                        {
                            for (j=0; j < 2; j++)
                            vst->momn[j][nrad-k] = ghost_st.momn[ind2[idir][j]];
                        }
                        else if (dim == 3)
                        {    
                            for (j = 0; j < 3; j++)
                                vst->momn[j][nrad-k] = ghost_st.momn[ind3[idir][j]];
                        }
                    }
                    break;
                default:
                    (void) printf("In appendGhostBuffer(): ");
                    (void) print_wave_type("Unknown wave type ",
                        wave_type(hs),"\n",front->interf);
                    (void) print_int_vector("icoords=", icoords,3,"\n");
                    LOC(); clean_up(ERROR);
                }

                break; //breaks out of for i loop
            }
            
        }
	    break;

	case 1:
	    for (i = 1; i <= nrad; ++i)  //nrad=3
	    {
            ic[idir] = icoords[idir] + i;
            index = d_index(ic,top_gmax,dim);

            //For debugging
            boolean status;
            
            //check neighbor in rdir[idir] 
            for (k = 0; k < dim; ++k)
                ic_next[k] = ic[k];

            ic_next[idir]--;
            int index_next = d_index(ic_next,top_gmax,dim);

            status = FT_StateStructAtGridCrossing(front,grid_intfc,
                    ic_next,rdir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);
            /*
            status = FT_StateStructAtGridCrossing(front,grid_intfc,
                    icoords,rdir[idir],comp,(POINTER*)&state,
                    &hs,crx_coords);
            */

            /*
            if (status && wave_type(hs) != DIRICHLET_BOUNDARY)
            {
                LOC();
                printf("In appendGhostBuffer() -- FT_StateStructAtGridCrossing()\n");
                printf("nb = %d\n",nb);
                printf("comp = %d n = %d\n",comp,n);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("top_comp[index_next] = %d\n",top_comp[index_next]);
                printf("positive_component(hs) = %d\n",positive_component(hs));
                printf("negative_component(hs) = %d\n",negative_component(hs));
                printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
                printf("wave_type(hs) = %d\n",wave_type(hs));
                printf("ic = %d %d %d\n",ic[0], ic[1], ic[2]);
                printf("icoords = %d %d %d\n",icoords[0], icoords[1], icoords[2]);
                printf("ic_next = %d %d %d\n\n",ic_next[0], ic_next[1], ic_next[2]);
            }
            else if (!status)
            {
                LOC();
                printf("In appendGhostBuffer() -- !FT_StateStructAtGridCrossing()\n");
                printf("nb = %d \t comp = %d\n",nb,comp);
                printf("cell_center[index].comp = %d\n",cell_center[index].comp);
                printf("top_comp[index] = %d\n",top_comp[index]);
                printf("top_comp[index_next] = %d\n",top_comp[index_next]);
                printf("ic = %d %d %d\n",ic[0], ic[1], ic[2]);
                printf("icoords = %d %d %d\n",icoords[0], icoords[1], icoords[2]);
                printf("ic_next = %d %d %d\n\n",ic_next[0], ic_next[1], ic_next[2]);
            }
            */

            //if ((comp == top_comp[index] && gas_comp(cell_center[index].comp)) || !status)
		    if (!needBufferFromIntfc(comp,cell_center[index].comp) && !status)
            {
                vst->dens[n+nrad+i-1] = m_vst->dens[index];
                vst->engy[n+nrad+i-1] = m_vst->engy[index];
                vst->pres[n+nrad+i-1] = m_vst->pres[index];
                
                for (j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = 0.0;
                
                if (dim == 1)
                {
                    vst->momn[0][n+nrad+i-1] = m_vst->momn[0][index];
                }
                else if (dim == 2)
                {
                    for(j = 0; j < 2; j++)
                    {
                        vst->momn[j][n+nrad+i-1] = m_vst->momn[ind2[idir][j]][index];
                    }
                }
                else if (dim == 3)
                {
                    for (j = 0; j < 3; j++)
                    {
                        vst->momn[j][n+nrad+i-1] = m_vst->momn[ind3[idir][j]][index];
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
                        
                        status = FT_StateStructAtGridCrossing(front,grid_intfc,
                                ic_tmp,rdir[idir],comp,(POINTER*)&state,
                                &hs,crx_coords);
                        
                        if (!status)
                        {
                            /* must be something wrong */
                            printf("In appendGhostBuffer() Case 1\n");
                            printf("ERROR: No crossing found!\n");
                            print_int_vector("ic=",ic,dim,"\n");
                            print_int_vector("ic_next=",ic_next,dim,"\n");
                            printf("direction: %s side %d comp %d\n",
                            grid_direction_name(rdir[idir]), nb, comp);
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
                    setNeumannStates(vst,m_vst,hs,state,ic_next,idir,nb,n,i,comp);
                    break;

                case ELASTIC_BOUNDARY:
                    setElasticStates(vst,m_vst,hs,state,ic_next,idir,nb,n,i,comp);
                    break;

                case DIRICHLET_BOUNDARY:
                    setDirichletStates(state,vst,m_vst,hs,ic_next,idir,nb,n,i,comp);
                    break;

                case FIRST_PHYSICS_WAVE_TYPE:
                    GFMGhostState(ic,comp,&ghost_st);
                    for (k = i; k <= nrad; ++k)
                    {
                        vst->dens[n+nrad+k-1] = ghost_st.dens;
                        vst->engy[n+nrad+k-1] = ghost_st.engy;
                        vst->pres[n+nrad+k-1] = ghost_st.pres;
                
                        for (j = 0; j < 3; j++)
                            vst->momn[j][n+nrad+k-1] = 0.0;

                        if (dim == 1)
                        {
                            vst->momn[0][n+nrad+k-1] = ghost_st.momn[0];
                        }
                        else if (dim == 2)
                        {
                            for (j = 0; j < 2; j++)
                                vst->momn[j][n+nrad+k-1] = ghost_st.momn[ind2[idir][j]];
                        }
                        else if (dim == 3)
                        {
                            for (j = 0; j < 3; j++)
                                vst->momn[j][n+nrad+k-1] = ghost_st.momn[ind3[idir][j]];
                        }
                    }
                    break;

                default:
                    (void) printf("In appendGhostBuffer(): ");
                    (void) print_wave_type("Unknown wave type ",
                    wave_type(hs),"\n",front->interf);
                    (void) print_int_vector("icoords=",icoords,3,"\n");
                    (void) printf("nb = %d\n",nb);
                    LOC(); clean_up(ERROR);
                }
                
                break; //breaks out of for i loop
            }
        }
	}

	if (debugging("trace"))
		printf("Leaving appendGhostBuffer()\n");
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
    if (eqn_params->porosity == 0 || !eqn_params->with_porosity)
    {
        setNeumannStates(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
        return;
    }

    switch (eqn_params->poro_scheme)
    {
        case PORO_SCHEME::DARCY:
            setElasticStatesDarcy(vst,m_vst,hs,state,icoords,idir,nb,n,istart,comp);
            break;
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

void CFABRIC_CARTESIAN::setElasticStatesDarcy(
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
    double crx_coords[MAXD] = {0.0};
    INTERFACE* grid_intfc = front->grid_intfc;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
    HYPER_SURF_ELEMENT* hse;
    POINTER sl, sr;
	
    STATE st_tmp_real;
	STATE st_tmp_ghost;	
    
    double vel_fluid[MAXD] = {0.0};
    
    int index = d_index(icoords,top_gmax,dim);
    
    double nor[MAXD];
    FT_NormalAtGridCrossing(front,icoords,dir[idir][nb],NO_COMP,nor,&hs,crx_coords);
	    //FT_NormalAtGridCrossing(front,icoords,dir[idir][nb],comp,nor,&hs,crx_coords);

    double nor_save[MAXD];
    for (int i = 0; i < dim; ++i)
        nor_save[i] = nor[i];


    double* vel_intfc = state->vel;
    

	int index_ghost;
	int ind2[2][2] = {{0,1},{1,0}};
    int ind3[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
	int ic_ghost[MAXD];

	double	coords[MAXD],coords_ref[MAXD],coords_ghost[MAXD];
	double	v[MAXD],v_ghost[MAXD];
	
	
	st_tmp_real.dim = dim;
	st_tmp_real.eos = &eqn_params->eos[GAS_COMP2];
	    //st_tmp_real.eos = &eqn_params->eos[comp];

    st_tmp_ghost.dim = dim;
    st_tmp_ghost.eos = &eqn_params->eos[GAS_COMP2];
        //st_tmp_ghost.eos = &eqn_params->eos[comp];

    int icoords_ghost[MAXD] = {0.0};
	for (int i = 0; i < dim; ++i)
	{
	    coords[i] = top_L[i] + icoords[i]*top_h[i];
	    icoords_ghost[i] = icoords[i];
	    ic_ghost[i] = icoords[i];
	}


    icoords_ghost[idir] = (nb == 0) ?  icoords[idir] - 1 : icoords[idir] + 1;
	for (int i = 0; i < dim; ++i)
	{
        coords_ghost[i] = top_L[i] + icoords_ghost[i]*top_h[i];
    }

    double vec_ghost[MAXD] = {0.0};
	for (int i = 0; i < dim; ++i)
	{
        vec_ghost[i] = coords_ghost[i] - crx_coords[i];
    }

    double side = Dotd(vec_ghost,nor,dim);
    if (side < 0)
    {
        for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
    }

    index_ghost = d_index(icoords_ghost,top_gmax,dim);
    if (side > 0)
    {
        cell_center[index_ghost].comp = 4;
    }
    else
    {
        cell_center[index_ghost].comp = 2;
    }


    //TODO: Can get rid of this loop.
    //      We compute a single ghost state and fill the entire stencil with it.
	for (int i = istart; i < istart + 1; ++i)
	{
	    //ghost point icoords and index
	    ic_ghost[idir] = (nb == 0) ?
            icoords[idir] - (i - istart + 1) : icoords[idir] + (i - istart + 1);

	    index_ghost = d_index(ic_ghost,top_gmax,dim);
        COMPONENT ghost_comp = top_comp[index_ghost];
            //COMPONENT ghost_comp = cell_center[index_ghost].comp;

        //ghost point coords
	    for (int j = 0; j < dim; ++j)
        {
            coords_ghost[j] = top_L[j] + ic_ghost[j]*top_h[j];
            coords_ref[j] = coords_ghost[j];
	    }

        /* Reflect ghost point through intfc-mirror at crossing */
        //first reflect across the grid line containing the intfc crossing 
	    double vn = 0.0;
	    coords_ref[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

	    for (int j = 0; j < dim; ++j)
	    {
		    v[j] = coords_ref[j] - crx_coords[j];
		    vn += v[j]*nor[j];
	    }
           
        //reflect v across the line containing the normal vector
	    for (int j = 0; j < dim; ++j)
		    v[j] = 2.0*vn*nor[j] - v[j];
	    
        //desired reflected point
        for (int j = 0; j < dim; ++j)
		    coords_ref[j] = crx_coords[j] + v[j];

        /*
        if (debugging("trace"))
        {
            (void) printf("\nEntered setElasticStatesDarcy():\n");
            (void) printf("comp = %d\n",comp);
            (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],icoords[2]);
            (void) printf("idir = %d nb = %d\n",idir,nb);
            (void) printf("istart = %d nrad = %d n = %d\n",istart,nrad,n);
            (void) print_general_vector("coords = ",coords,dim,"\n");
            (void) print_general_vector("coords_reflect = ",coords_ref,dim,"\n");
            (void) print_general_vector("crx_coords = ",crx_coords,dim,"\n");
            (void) print_general_vector("nor = ",nor,dim,"\n");
            (void) print_general_vector("vel_intfc = ",vel_intfc,dim,"\n");
        }
        */

        /* Interpolate the state at the reflected point */ 
        
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->dens,getStateDens,&st_tmp_ghost.dens,&m_vst->dens[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->pres,getStatePres,&st_tmp_ghost.pres,&m_vst->pres[index]);
	    
        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->temp,getStateTemp,&st_tmp_ghost.temp,&m_vst->temp[index]);

        FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->engy,getStateEngy,&st_tmp_ghost.engy,&m_vst->engy[index]);
	    
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                m_vst->mu,getStateMu,&st_tmp_ghost.mu,&m_vst->mu[index]);
	    //FT_IntrpStateVarAtCoords(front,comp,coords_ref,
          //      m_vst->mu_turb,getStateMuTurb,&st_tmp_ghost.mu_turb,&m_vst->mu_turb[index]);
	    FT_IntrpStateVarAtCoords(front,comp,coords_ref,
                field.mu_turb,getStateMuTurb,&st_tmp_ghost.mu_turb,&field.mu_turb[index]);

        
        double mu_total = st_tmp_ghost.mu + st_tmp_ghost.mu_turb;
        
        double v_reflect[3];
        for (int j = 0; j < dim; ++j)
        {
            FT_IntrpStateVarAtCoords(front,comp,coords_ref,m_vst->momn[j],
                    getStateMom[j],&st_tmp_ghost.momn[j],&m_vst->momn[j][index]);
            v_reflect[j] = st_tmp_ghost.momn[j]/st_tmp_ghost.dens;
        }

        /*
        if (front->step < eqn_params->fsi_startstep)
        {
            break;
        }
        */

	    st_tmp_real.dens = m_vst->dens[index_ghost];
	    st_tmp_real.pres = m_vst->pres[index_ghost];
	    st_tmp_real.temp = m_vst->temp[index_ghost];

        double mu_total_real = m_vst->mu[index_ghost] + field.mu_turb[index_ghost];
        //double mu_total_real = m_vst->mu[index_ghost] + m_vst->mu_turb[index_ghost];
        
        double v_real[MAXD];
        for (int j = 0; j < dim; j++)
        {
            st_tmp_real.momn[j] = m_vst->momn[j][index_ghost];
            v_real[j] = st_tmp_real.momn[j]/st_tmp_real.dens;
        }
        
        double vel_rel_real[MAXD] = {0.0};
        double vn_real = 0.0;
        for (int j = 0; j < dim; j++)
	    {
            vel_rel_real[j] = v_real[j] - vel_intfc[j];
            vn_real += (v_real[j] - vel_intfc[j])*nor[j];
	    }

        double vel_rel_real_tan[MAXD] = {0.0};
        for (int j = 0; j < dim; j++)
        {
            vel_rel_real_tan[j] = vel_rel_real[j] - vn_real*nor[j];
        }

        double pl = st_tmp_ghost.pres;
        double rhol = st_tmp_ghost.dens;

        double pr = st_tmp_real.pres;
        double rhor = st_tmp_real.dens;
        
        double gamma = st_tmp_ghost.eos->gamma;

        double Msqr = gamma/(gamma + 1.0)*std::abs(rhor*pr - rhol*pl);

        double poro = eqn_params->porosity;
        double alpha = eqn_params->porous_coeff[0];
        double beta = eqn_params->porous_coeff[1];
        
        //TODO: OR SHOULD INTERPOLATE TEMPERATURES AND USE TO COMPUTE VISCOSITIES

        //double A = mu_total*alpha;

        double A = 0.5*(sqr(m_vst->mu[index_ghost]) - sqr(st_tmp_ghost.mu))*alpha;
            //double A = 0.5*(sqr(mu_total_real) - sqr(mu_total))/k_perm;
        //double B = rhol*beta;
        
        double sgn = (rhor*pr - rhol*pl >= 0) ? 1.0 : -1.0;

        double mdot = -2.0*sgn*Msqr/(A + std::sqrt(A*A + 4.0*beta*Msqr));

        double ghost_nor_vel = 0.0;
        if (std::abs(mdot) > MACH_EPS)
        {
            ghost_nor_vel = mdot/rhor;
        }

        st_tmp_ghost.dens = rhor;
        st_tmp_ghost.pres = pr;

        double velo[MAXD] = {0.0};
        for (int j = 0; j < dim; ++j)
        {
            //TODO: NEED TO HANDLE TANGENTIAL JUMP?
            velo[j] = vel_rel_real_tan[j] + ghost_nor_vel*nor[j] + vel_intfc[j];
            st_tmp_ghost.vel[j] = velo[j];
            st_tmp_ghost.momn[j] = st_tmp_ghost.dens*velo[j];
        }

        st_tmp_ghost.temp = st_tmp_real.temp;
        st_tmp_ghost.engy = EosEnergy(&st_tmp_ghost);
        
	    
        ////////////////////////////////////////////////////////////////
        if (debugging("elastic_buffer"))
        {
            double debug_coords[MAXD] = {0.6,0.6,0.42354};
            for (int j = 0; j < dim; ++j)
            {
                debug_coords[j] -= crx_coords[j];
            }
            debug_coords[2] = 0.0;
            if (Mag3d(debug_coords) < 3.0e-01)
            {
                printf("\nDARCY_DEBUG\n");
                printf("\nistart = %d nrad = %d n = %d\n",istart,nrad,n);
                printf("coords = %f %f %f\n",coords[0],coords[1],coords[2]);
                printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
                printf("coords_ghost = %f %f %f\n",coords_ghost[0],coords_ghost[1],coords_ghost[2]);
                printf("idir: %d \t nb: %d\n",idir,nb);
                printf("v_real = %f %f %f\n",v_real[0],v_real[1],v_real[2]);
                printf("rhol = %g pl = %g\n", rhol, pl);
                printf("rhor = %g pr = %g\n", rhor, pr);
                printf("rhor*pr - rhol*pl = %g \t sgn = %f\n",rhor*pr - rhol*pl,sgn);
                printf("Msqr = %g \t mdot = %g\n", Msqr, mdot);
                printf("nor = %f %f %f\n", nor[0],nor[1],nor[2]);
                printf("ghost_pres = %f\n", st_tmp_ghost.pres);
                printf("ghost_dens = %f\n", st_tmp_ghost.dens);
                printf("ghost_temp = %f\n", st_tmp_ghost.temp);
                printf("ghost_engy = %f\n", st_tmp_ghost.engy);
                printf("ghost_momn = %f %f %f\n",
                    st_tmp_ghost.momn[0],st_tmp_ghost.momn[1],st_tmp_ghost.momn[2]);
                printf("ghost_vel = %f %f %f\n",
                    st_tmp_ghost.vel[0],st_tmp_ghost.vel[1],st_tmp_ghost.vel[2]);
            }
        }
        ////////////////////////////////////////////////////////////////


	    if (std::isnan(st_tmp_ghost.momn[0]))
	    {
            //printf("negative energy! \n");
            printf("icoords = %d %d %d \n",icoords[0],icoords[1],icoords[2]);
            printf("coords = %f %f %f\n",coords[0],coords[1],coords[2]);
            printf("crx_coords = %f %f %f\n",crx_coords[0],crx_coords[1],crx_coords[2]);
            printf("coords_ghost = %f %f %f\n",coords_ghost[0],coords_ghost[1],coords_ghost[2]);
            printf("idir: %d \t nb: %d\n",idir,nb);
            printf("%f %f %f %f %f %f \n",st_tmp_ghost.dens,st_tmp_ghost.momn[0],
                st_tmp_ghost.momn[1],st_tmp_ghost.momn[2],st_tmp_ghost.pres,
                st_tmp_ghost.engy);
            printf("st_tmp_ghost.dim = %d, idir = %d, nb = %d \n",
                st_tmp_ghost.dim,idir,nb);
            printf("gamma = %f, einf = %f, pinf = %f \n",st_tmp_ghost.eos->gamma,
                st_tmp_ghost.eos->einf,st_tmp_ghost.eos->pinf);
            printf("coords_ref = %f %f %f \n",coords_ref[0],coords_ref[1],
                            coords_ref[2]);
            printf("v_real = %f %f %f\n",v_real[0],v_real[1],v_real[2]);
            printf("rhol = %f \t pl = %f \t rhor = %f \t pr = %f\n",rhol,pl,rhor,pr);
            printf("rhor*pr - rhol*pl = %f \t sgn = %f\n",rhor*pr - rhol*pl,sgn);
            printf("Msqr = %f \t mdot = %f\n", Msqr, mdot);
            printf("nor = %f %f %f\n", nor[0],nor[1],nor[2]);
            printf("ghost_pres = %f\n", st_tmp_ghost.pres);
            printf("ghost_dens = %f\n", st_tmp_ghost.dens);
            printf("ghost_temp = %f\n", st_tmp_ghost.temp);
            printf("ghost_engy = %f\n", st_tmp_ghost.engy);
            printf("ghost_momn = %f %f %f\n",
                st_tmp_ghost.momn[0],st_tmp_ghost.momn[1],st_tmp_ghost.momn[2]);
            printf("ghost_vel = %f %f %f\n",
                st_tmp_ghost.vel[0],st_tmp_ghost.vel[1],st_tmp_ghost.vel[2]);
            LOC(); clean_up(EXIT_FAILURE);
	    }
    }

    for (int i = istart; i <= nrad; ++i)
    {
        if (nb == 0)
        {
            vst->dens[nrad-i] = st_tmp_ghost.dens;
            vst->engy[nrad-i] = st_tmp_ghost.engy;
            vst->pres[nrad-i] = st_tmp_ghost.pres;
            
            for (int j = 0; j < 3; j++)
                vst->momn[j][nrad-i] = 0.0;

            if (dim == 3)
            {
                for (int j = 0; j < 3; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind3[idir][j]];
            }
            else if (dim == 2)
            {
                for (int j = 0; j < 2; j++)
                    vst->momn[j][nrad-i] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][nrad-i] = st_tmp_ghost.momn[0];
            }
        }
        else
        {
            vst->dens[n+nrad+i-1] = st_tmp_ghost.dens;
            vst->engy[n+nrad+i-1] = st_tmp_ghost.engy;
            vst->pres[n+nrad+i-1] = st_tmp_ghost.pres;
            
            for (int j = 0; j < 3; j++)
                vst->momn[j][n+nrad+i-1] = 0.0;
    
            if (dim == 3)
            {
                for (int j = 0; j < 3; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind3[idir][j]];
            }
            else if (dim == 2)
            {
                for (int j = 0; j < 2; j++)
                    vst->momn[j][n+nrad+i-1] = st_tmp_ghost.momn[ind2[idir][j]];
            }
            else
            {
                vst->momn[0][n+nrad+i-1] = st_tmp_ghost.momn[0];
            }
        }
    }

    /*
	if (debugging("trace"))
        (void) printf("Leaving setElasticStatesDarcy()\n\n");
    */
}	/* end setElasticStatesDarcy */


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
	    (void) printf("icoords = %d %d %d\n",icoords[0],icoords[1],icoords[2]);
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
                coords_ghost,NO_BOUNDARIES,crx_coords,intrp_a,
                &nearHse,&nearHs,5);
        /*status = FT_FindNearestIntfcPointInRange(front,comp_ghost,
                coords_ghost,INCLUDE_BOUNDARIES,crx_coords,intrp_a,
                &nearHse,&nearHs,5);*/
        
        
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

void CFABRIC_CARTESIAN::setViscousGhostState(
        int* icoords,
        COMPONENT comp,
        VSWEEP* vs,
        SWEEP* m_vst)
{
    double nip_coords[MAXD];
    double intrp_coeffs[MAXD];
    HYPER_SURF* hs;
    HYPER_SURF_ELEMENT* hse;
    
    int ghost_index = d_index(vs->icoords,top_gmax,dim);
    auto ghost_coords = cell_center[ghost_index].getCoords();
    COMPONENT ghost_comp = vs->comp;
    
    INTERFACE* intfc = front->grid_intfc;
    if (!gas_comp(ghost_comp))
    {
        intfc = front->interf;
    }


    bool nip_found = nearest_interface_point(&ghost_coords[0],
                comp,intfc,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);

    /*
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                ghost_comp,intfc,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    */


    /*
    bool nip_found = nearest_interface_point(&ghost_coords[0],
                ghost_comp,front->grid_intfc,NO_SUBDOMAIN,nullptr,
                nip_coords,intrp_coeffs,&hse,&hs);
    */
    

    if (debugging("viscous_ghost"))
    {
        printf("\nsetViscousGhostState() DEBUGGING\n");
        int index = d_index(icoords,top_gmax,dim);
        auto coords = cell_center[index].getCoords();
        fprint_general_vector(stdout,"coords",&coords[0],dim,"\n");
        fprint_general_vector(stdout,"ghost_coords",&ghost_coords[0],dim,"\n");
        printf("comp = %d ghost_comp = %d\n", comp,ghost_comp);
        fprint_general_vector(stdout,"coords_nip",nip_coords,dim,"\n");
        printf("wave_type(hs) = %d\n",wave_type(hs));
    }
    
    if (!nip_found)
    {
        printf("ERROR G_CARTESIAN::setViscousGhostState(): "
                "can't find nearest interface point\n");
        LOC(); clean_up(EXIT_FAILURE);
    }

    switch (wave_type(hs))
    {
        case DIRICHLET_BOUNDARY:
        {
            setDirichletViscousGhostState(vs,comp,intrp_coeffs,hse,hs);
            break;
        }
        case MOVABLE_BODY_BOUNDARY:
        case NEUMANN_BOUNDARY:
        {
            setNeumannViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                    nip_coords,comp,intrp_coeffs,hse,hs);
            break;
        }
        case ELASTIC_BOUNDARY: //TODO: DOES THE SOLVER EVER SEE THIS BOUNDARY?
        {
            if (eqn_params->porosity == 0 || !eqn_params->with_porosity)
            {
                setNeumannViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                        nip_coords,comp,intrp_coeffs,hse,hs);
            }
            else
            {
                setElasticViscousGhostState(icoords,m_vst,vs,&ghost_coords[0],
                        nip_coords,comp,intrp_coeffs,hse,hs);
            }
            break;
        }
        default:
        {
            printf("\n\nsetViscousGhostState() ERROR: unknown boundary type\n\n");
            LOC(); clean_up(EXIT_FAILURE);
        }
    }

}

void CFABRIC_CARTESIAN::setElasticViscousGhostState(
        int* icoords,
        SWEEP* m_vst,
        VSWEEP* vs,
        double* ghost_coords,
        double* crx_coords,
        COMPONENT comp,
        double* intrp_coeffs,
        HYPER_SURF_ELEMENT* hse,
        HYPER_SURF* hs)
{

    int ghost_index = d_index(vs->icoords,top_gmax,dim);
    COMPONENT ghost_comp = vs->comp;
    
    double nor[MAXD] = {0.0};
    double vel_intfc[MAXD] = {0.0};

    switch (dim)
    {
    case 2:
        {
            STATE* ss;
            STATE* se;

            if (gas_comp(negative_component(hs)))
            {
                ss = (STATE*)left_state(Bond_of_hse(hse)->start);
                se = (STATE*)left_state(Bond_of_hse(hse)->end);
            }
            else if (gas_comp(positive_component(hs)))
            {
                ss = (STATE*)right_state(Bond_of_hse(hse)->start);
                se = (STATE*)right_state(Bond_of_hse(hse)->end);
            }
            else
            {
                printf("setElasticViscousGhostState() ERROR: "
                        "no gas component on hypersurface\n");
                LOC(); clean_up(EXIT_FAILURE);
            }

            double ns[MAXD] = {0.0};
            double ne[MAXD] = {0.0};
            
            normal(Bond_of_hse(hse)->start,hse,hs,ns,front);
            normal(Bond_of_hse(hse)->end,hse,hs,ne,front);

            for (int i = 0; i < dim; ++i)
            {
                nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
            }
        }
        break;

    case 3:
        {
            STATE* st[3];
            TRI* nearTri = Tri_of_hse(hse);

            if (comp == negative_component(hs))
            {
                for (int j = 0; j < 3; ++j)
                    st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
            }
            else if (comp == positive_component(hs))
            {
                for (int j = 0; j < 3; ++j)
                    st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
            }
            else
            {
                printf("setElasticViscousGhostState() ERROR: "
                        "no gas component on hypersurface\n");
                LOC(); clean_up(EXIT_FAILURE);
            }

            //NOTE: Tri_normal() does not return a unit vector
            const double* tnor = Tri_normal(nearTri);

            for (int i = 0; i < dim; ++i)
            {
                nor[i] = tnor[i];

                vel_intfc[i] = 0.0;
                for (int j = 0; j < 3; ++j)
                    vel_intfc[i] += intrp_coeffs[j]*st[j]->vel[i];
            }
        }
        break;
    }

    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hs))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    
    double dist_reflect = FT_GridSizeInDir(nor,front);
    
    double dist_ghost = distance_between_positions(ghost_coords,crx_coords,dim);
    
    //The desired reflected point
    double coords_reflect[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    }

    int index = d_index(icoords,top_gmax,dim);
    //TODO: ADD DEBUG BLOCK HERE -- LIKE setSlipBoundaryNIP()


    //Interpolate Density and Momentum at the reflected point and compute the velocity.
    double dens_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->dens,
            getStateDens,&dens_reflect,&m_vst->dens[index]);

    double mom_reflect[MAXD];
    double vel_reflect[MAXD];
    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->momn[j],
                getStateMom[j],&mom_reflect[j],&m_vst->momn[j][index]);
        
        vel_reflect[j] = mom_reflect[j]/dens_reflect;
    }

    double pres_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->pres,
            getStatePres,&pres_reflect,&m_vst->pres[index]);
    
    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->mu,
            getStateMu,&mu_reflect,&m_vst->mu[index]);
    
    double mu_turb_reflect;
    //FT_IntrpStateVarAtCoords(front,comp,coords_reflect,m_vst->mu_turb,
      //      getStateMuTurb,&mu_turb_reflect,&m_vst->mu_turb[index]);
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field.mu_turb,
            getStateMuTurb,&mu_turb_reflect,&field.mu_turb[index]);

    double mu_total_reflect = mu_reflect + mu_turb_reflect;
    
    double dens_real = m_vst->dens[ghost_index];
    double v_real[MAXD];
    for (int j = 0; j < dim; j++)
    {
        v_real[j] = m_vst->momn[j][ghost_index]/dens_real;
    }
    
    double vel_rel_real[MAXD] = {0.0};
    double vn_real = 0.0;
    for (int j = 0; j < dim; j++)
    {
        vel_rel_real[j] = v_real[j] - vel_intfc[j];
        vn_real += (v_real[j] - vel_intfc[j])*nor[j];
    }

    double vel_rel_real_tan[MAXD] = {0.0};
    for (int j = 0; j < dim; j++)
    {
        vel_rel_real_tan[j] = vel_rel_real[j] - vn_real*nor[j];
    }

    
    double pl = pres_reflect;
    double rhol = dens_reflect;

    double pr = m_vst->pres[ghost_index];
    double rhor = m_vst->dens[ghost_index];
    
    EOS_PARAMS eos = eqn_params->eos[GAS_COMP2];
    double gamma = eos.gamma;

    double Msqr = gamma/(gamma + 1.0)*std::abs(rhor*pr - rhol*pl);

    double poro = eqn_params->porosity;
    double alpha = eqn_params->porous_coeff[0];
    double beta = eqn_params->porous_coeff[1];
        
        //TODO: OR SHOULD INTERPOLATE TEMPERATURES AND USE TO COMPUTE VISCOSITIES

    double mu_ghost =  m_vst->mu[ghost_index];
    double mu_total_ghost = m_vst->mu[ghost_index] + field.mu_turb[ghost_index];
    //double mu_total_ghost = m_vst->mu[ghost_index] + m_vst->mu_turb[ghost_index];
    
    double A = 0.5*(sqr(mu_ghost) - sqr(mu_reflect))*alpha;
    //double A = 0.5*(sqr(mu_total_ghost) - sqr(mu_total_reflect))*alpha;
        //double A = 0.5*(sqr(mu_total_real) - sqr(mu_total))/k_perm;
    //double B = rhol*beta;
        
    double sgn = (rhor*pr - rhol*pl >= 0) ? 1.0 : -1.0;

    double mdot = -2.0*sgn*Msqr/(A + std::sqrt(A*A + 4.0*beta*Msqr));

    double ghost_nor_vel = 0.0;
    if (std::abs(mdot) > MACH_EPS)
    {
        ghost_nor_vel = mdot/rhor;
    }

    //if (front->step < fsi_startstep)

    double velo[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        //TODO: NEED TO HANDLE TANGENTIAL JUMP?
        velo[j] = vel_rel_real_tan[j] + ghost_nor_vel*nor[j] + vel_intfc[j];

        vs->vel[j] = velo[j];
    }

    
    vs->mu = m_vst->mu[ghost_index];
    vs->dens = m_vst->dens[ghost_index];
    vs->temp = m_vst->temp[ghost_index];
        

    if (debugging("elastic_viscous_ghost"))
    {
        printf("setElasticViscousGhostState():\n");
        int index = d_index(icoords,top_gmax,dim);
        auto coords = cell_center[index].getCoords();
        printf("comp = %d ghost_comp = %d\n", comp,ghost_comp);
        printf("cc-comp = %d cc-ghost_comp = %d\n",
                cell_center[index].comp,cell_center[ghost_index].comp);
        fprint_general_vector(stdout,"coords",&coords[0],dim,"\n");
        fprint_general_vector(stdout,"ghost_coords",ghost_coords,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",crx_coords,dim,"\n");
        printf("\nsetElasticViscousGhostState() DEBUGGING\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_reflect = %f \t dist_ghost = %f\n",dist_reflect,dist_ghost);
        printf("pres_reflect = %f \t pres_real = %f\n",pl,pr);
        printf("dens_reflect = %f \t dens_real = %f\n",rhol,rhor);
        printf("Msqr = %g \t mdot = %g\n", Msqr, mdot);
        fprint_general_vector(stdout,"vel_ghost",velo,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
    }
}

