#include "iFluid.h"


//TODO: remove after verify not needed
extern int ifluid_find_state_at_dual_crossing(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	int comp,
	POINTER *state,
	HYPER_SURF **hs,
	double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->comp_grid_intfc;
	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);
	if (status == NO) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == NEUMANN_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == MOVABLE_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == ICE_PARTICLE_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY) 
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_dual_crossing */

extern int ifluid_find_state_at_cg_crossing(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	int comp,
	POINTER *state,
	HYPER_SURF **hs,
	double *crx_coords)
{
	boolean status;
	INTERFACE *grid_intfc = front->comp_grid_intfc;
	status = FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir,
				comp,state,hs,crx_coords);

	if (status == NO) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == FIRST_PHYSICS_WAVE_TYPE) 
	    return NO_PDE_BOUNDARY;
	if (wave_type(*hs) == NEUMANN_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == GROWING_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == MOVABLE_BODY_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == ICE_PARTICLE_BOUNDARY) 
	    return CONST_V_PDE_BOUNDARY;
	if (wave_type(*hs) == DIRICHLET_BOUNDARY) 
	{
	    if (boundary_state(*hs))
	    	return CONST_V_PDE_BOUNDARY;
	    else
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_crossing */

void Incompress_Solver_Smooth_Basis::setDualDomain()
{
	static boolean first = YES;
	INTERFACE *comp_grid_intfc;
	Table *T;
	int i,size,comp_size;

	comp_grid_intfc = front->comp_grid_intfc;
	ctop_grid = &topological_grid(comp_grid_intfc);
	ctop_gmax = ctop_grid->gmax;
	ctop_L = ctop_grid->L;
	ctop_U = ctop_grid->U;
	T = table_of_interface(comp_grid_intfc);
	ctop_comp = T->components;
	cimin = (lbuf[0] == 0) ? 1 : lbuf[0] + 1;
	cjmin = (lbuf[1] == 0) ? 1 : lbuf[1] + 1;
	ckmin = (lbuf[2] == 0) ? 1 : lbuf[2] + 1;
	cimax = (ubuf[0] == 0) ? ctop_gmax[0] - 1 : ctop_gmax[0] - ubuf[0];
	cjmax = (ubuf[1] == 0) ? ctop_gmax[1] - 1 : ctop_gmax[1] - ubuf[1];
	ckmax = (ubuf[2] == 0) ? ctop_gmax[2] - 1 : ctop_gmax[2] - ubuf[2];
	for (i = 0; i < dim; ++i)
	    offset[i] = (lbuf[i] == 0) ? 1 : 0;

	comp_size = ctop_gmax[0]+1;
        for (i = 1; i < dim; ++i)
            comp_size *= (ctop_gmax[i]+1);
	size = top_gmax[0]+1;
        for (i = 1; i < dim; ++i)
            size *= (top_gmax[i]+1);

	if (first)
	{
	    iFparams->field = field;
	    FT_VectorMemoryAlloc((POINTER*)&field->d_phi,comp_size,
					sizeof(double));
	    if (size < comp_size)
	    {
		FT_FreeThese(4,source,diff_coeff,field->div_U,array);
            	FT_VectorMemoryAlloc((POINTER*)&source,comp_size,
					sizeof(double));
            	FT_VectorMemoryAlloc((POINTER*)&diff_coeff,comp_size,
					sizeof(double));
            	FT_VectorMemoryAlloc((POINTER*)&field->div_U,comp_size,
                        		sizeof(double));
            	FT_VectorMemoryAlloc((POINTER*)&array,comp_size,
					sizeof(double));
	    }
	    first = NO;
	}
}	/* end setDualDomain */

void Incompress_Solver_Smooth_Basis::setDualGlobalIndex()
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&cn_dist,num_nodes,sizeof(int));
	}
	cNLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index1d(i,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		cNLblocks++;
	    }
	    break;
	case 2:
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index2d(i,j,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		cNLblocks++;
	    }
	    break;
	case 3:
	    for (k = ckmin; k <= ckmax; k++)
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index3d(i,j,k,ctop_gmax);
		if (ctop_comp[ic] == SOLID_COMP) continue;
		//if (domain_status[ic] != TO_SOLVE) continue;
		cNLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	cn_dist[myid] = cNLblocks;
	pp_global_imax(cn_dist,num_nodes);
	cilower = 0;
        ciupper = cn_dist[0];

        for (i = 1; i <= myid; i++)
        {
            cilower += cn_dist[i-1];
            ciupper += cn_dist[i];
        }	
}	/* setDualGlobalIndex */

void Incompress_Solver_Smooth_Basis::setDualIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int size = ciupper - cilower;
	static int old_size;

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&cij_to_I,ctop_gmax[0]+1,
				ctop_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&cijk_to_I,ctop_gmax[0]+1,
				ctop_gmax[1]+1,ctop_gmax[2]+1,INT);
	    	break;
	    }
	    old_size = size;
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] + 1: 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= ctop_gmax[1]; j++)
	    for (i = 0; i <= ctop_gmax[0]; i++)
		    cij_to_I[i][j] = -1;
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index2d(i,j,ctop_gmax);
                if (ctop_comp[ic] != SOLID_COMP)
                {
                    cij_to_I[i][j] = index + cilower;
                    index++;
                }
	    }
	    FT_ParallelExchCompGridCellIndex(front,llbuf,uubuf,
				(POINTER)cij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= ctop_gmax[2]; k++)
	    for (j = 0; j <= ctop_gmax[1]; j++)
	    for (i = 0; i <= ctop_gmax[0]; i++)
		    cijk_to_I[i][j][k] = -1;
	    for (k = ckmin; k <= ckmax; k++)
	    for (j = cjmin; j <= cjmax; j++)
	    for (i = cimin; i <= cimax; i++)
	    {
		ic = d_index3d(i,j,k,ctop_gmax);
		/*
		if (domain_status[ic] != TO_SOLVE)
		    continue;
		*/
                if (ctop_comp[ic] != SOLID_COMP)
		{
                    cijk_to_I[i][j][k] = index + cilower;
                    index++;
                }
	    }
	    FT_ParallelExchCompGridCellIndex(front,llbuf,uubuf,
				(POINTER)cijk_to_I);
	    break;
	}
}	/* end setDualIndexMap */

void Incompress_Solver_Smooth_Basis::computeDualFieldPointGrad(
        int *icoords,
        double *field,
        double *grad_field)
{
	int i,j,k,index,index_l,index_u;
	int dual_icl[MAXD],dual_icu[MAXD];
	double pu,pl;
	COMPONENT comp,dual_compl,dual_compu;
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int status;
	double denom;

	index = d_index(icoords,top_gmax,dim);
	comp = top_comp[index];
	switch (dim)
	{
	case 2:
	    for (i = 0; i < dim; ++i)
	    {
		denom = 0.0;
	    	pu = pl = 0.0;
	    	dual_icu[i] = icoords[i] - offset[i] + 1;
	    	dual_icl[i] = icoords[i] - offset[i];
	    	for (j = 0; j <= 1; ++j)
	    	{
		    dual_icu[(i+1)%dim] = icoords[(i+1)%dim] - 
				offset[(i+1)%dim] + j;
		    dual_icl[(i+1)%dim] = icoords[(i+1)%dim] - 
				offset[(i+1)%dim] + j;
		    index_u = d_index(dual_icu,ctop_gmax,dim);
		    index_l = d_index(dual_icl,ctop_gmax,dim);
		    dual_compu = ctop_comp[index_u];
		    dual_compl = ctop_comp[index_l];
		    if (dual_compu == comp && dual_compl == comp &&
                                !ifluid_comp(comp))
                        continue;	/*Common all solid case*/
		    else if (dual_compl != comp)
		    {
			status = (*findStateAtCGCrossing)(front,dual_icu,
				dir[i][0],comp,&intfc_state,&hs,crx_coords);
			if (status == CONST_P_PDE_BOUNDARY)
                        {
                            pu += field[index_u];
                            pl += getStatePhi(intfc_state);
                            denom += 1.0;
                        }
			else	 /*CONST_V_PDE_BOUNDARY and NO_PDE_BOUNDARY*/ 
			{
			    if (dual_compl != dual_compu)
                            {
                                denom += 1.0;
                                continue;
                            }
                            pu += field[index_u];
                            pl += field[index_l];
                            denom += 1.0;

			}
		    }
		    else if (dual_compu != comp)
		    {
			status = (*findStateAtCGCrossing)(front,dual_icl,
				dir[i][1],comp,&intfc_state,&hs,crx_coords);
			if (status == CONST_P_PDE_BOUNDARY)
                        {
                            pu += getStatePhi(intfc_state);
                            pl += field[index_l];
                            denom += 1.0;
                        }
			else	/*CONST_V_PDE_BOUNDARY and NO_PDE_BOUNDARY*/
			{
			    if (dual_compl != dual_compu)
                            {
                                denom += 1.0;
                                continue;
                            }
                            pu += field[index_u];
                            pl += field[index_l];
                            denom += 1.0;
			}
		    }
		    else	/*Common all fluid case*/
		    {
		    	pu += field[index_u];
		    	pl += field[index_l];
		    	denom += 1.0;
		    }
	    	}
		grad_field[i] = (denom == 0.0) ? 0.0:(pu - pl)/top_h[i]/denom;
	    }
	    break;
	case 3:
	    for (i = 0; i < dim; ++i)
	    {
		denom = 0.0;
	    	pu = pl = 0.0;
	    	dual_icu[i] = icoords[i] - offset[i] + 1;
	    	dual_icl[i] = icoords[i] - offset[i];
	    	for (j = 0; j <= 1; ++j)
	    	for (k = 0; k <= 1; ++k)
	    	{
		    dual_icl[(i+1)%dim] = dual_icu[(i+1)%dim] = 
			icoords[(i+1)%dim] - offset[(i+1)%dim] + j;
		    dual_icl[(i+2)%dim] = dual_icu[(i+2)%dim] = 
			icoords[(i+2)%dim] - offset[(i+2)%dim] + k;
		    index_u = d_index(dual_icu,ctop_gmax,dim);
		    index_l = d_index(dual_icl,ctop_gmax,dim);
		    dual_compu = ctop_comp[index_u];
                    dual_compl = ctop_comp[index_l];
                    if (dual_compu == comp && dual_compl == comp &&
				!ifluid_comp(comp))
                        continue;	/*Common all solid case*/
                    else if (dual_compl != comp)
                    {
                        status = (*findStateAtCGCrossing)(front,dual_icu,
                                dir[i][0],comp,&intfc_state,&hs,crx_coords);
			if (status == CONST_P_PDE_BOUNDARY)
                        {
                            pu += field[index_u];
                            pl += getStatePhi(intfc_state);
                            denom += 1.0;
                        }
			else	/*CONST_V_PDE_BOUNDARY and NO_PDE_BOUNDARY*/
			{
			    if (dual_compl != dual_compu)
                            {
                                denom += 1.0;
                                continue;
                            }
                            pu += field[index_u];
                            pl += field[index_l];
                            denom += 1.0;
			}
		    }
		    else if (dual_compu != comp)
                    {
                        status = (*findStateAtCGCrossing)(front,dual_icl,
                                dir[i][1],comp,&intfc_state,&hs,crx_coords);
			if (status == CONST_P_PDE_BOUNDARY)
                        {
                            pu += getStatePhi(intfc_state);
                            pl += field[index_l];
                            denom += 1.0;
                        }
			else	/*CONST_V_PDE_BOUNDARY and NO_PDE_BOUNDARY*/
			{
			    if (dual_compl != dual_compu)
                            {
                                denom += 1.0;
                                continue;
                            }
                            pu += field[index_u];
                            pl += field[index_l];
                            denom += 1.0;
			}
                    }
                    else	/*Common all fluid case*/
                    {
                        pu += field[index_u];
                        pl += field[index_l];
                        denom += 1.0;
                    }
	    	}
		grad_field[i] = (denom == 0.0) ? 0.0 : (pu - pl)/top_h[i]/denom;
	    }
	    break;
	}
}	/* end computeDualFieldPointGrad */

double Incompress_Solver_Smooth_Basis::computeDualFieldPointDiv(
        int *icoords,
        double **field)
{
        int i,j,k,index,index_u,index_l;
        COMPONENT comp_l,comp_u;
	double div;
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	int status;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int idir,nb;
	int ic[MAXD],icu[MAXD],icl[MAXD];
	double du;

	for (i = 0; i < dim; ++i)
	{
	    ic[i] = icoords[i] - 1 + offset[i];
	}
	div = 0.0;
	for (i = 0; i < dim; ++i)
        {
	    du = 0.0;
            icl[i] = ic[i];
            icu[i] = ic[i] + 1;
            for (j = 0; j <= 1; ++j)
	    {
		icl[(i+1)%dim] = ic[(i+1)%dim] + j;
		icu[(i+1)%dim] = ic[(i+1)%dim] + j;
		if (dim == 3)
		{
		    for (k = 0; k <= 1; ++k)
		    {
			icl[(i+2)%dim] = ic[(i+2)%dim] + k;
                	icu[(i+2)%dim] = ic[(i+2)%dim] + k;
                	index_l = d_index(icl,top_gmax,dim);
                	index_u = d_index(icu,top_gmax,dim);
		    	comp_l = ctop_comp[index_l];
		    	comp_u = ctop_comp[index_l];
		    	if (comp_l != comp_u)
			    printf("Component not equal!\n");
			du += (field[i][index_u] - field[i][index_l]);
		    }
		}
		else
		{
		    index_l = d_index(icl,top_gmax,dim);
		    index_u = d_index(icu,top_gmax,dim);
		    comp_l = ctop_comp[index_l];
		    comp_u = ctop_comp[index_l];
		    if (comp_l != comp_u)
			printf("Component not equal!\n");
		    du += (field[i][index_u] - field[i][index_l]);
		}
	    }
            du /= pow(2,dim-1);
	    div += du/top_h[i];
        }
        return div;
}       /* end computeDualFieldPointDiv */

double Incompress_Solver_Smooth_Basis::computeDualMu(
        int *icoords,
        double *mu)
{
	int icl[MAXD],icu[MAXD],ic[MAXD];
	int i,j,k,n;
	int index,dual_index;
	double ave_mu = 0.0;

	dual_index = d_index(icoords,ctop_gmax,dim);
	for (i = 0; i < dim; ++i)
	{
	    icl[i] = icoords[i] - 1 + offset[i];
	    icu[i] = icoords[i] + offset[i];
	}
	n = 0;
	switch (dim)
	{
	case 2:
	    for (ic[0] = icl[0]; ic[0] <= icu[0]; (ic[0])++)
	    for (ic[1] = icl[1]; ic[1] <= icu[1]; (ic[1])++)
	    {
		index = d_index(ic,top_gmax,dim);
		if (top_comp[index] != ctop_comp[dual_index])
		    continue;
		ave_mu += mu[index];
		n++;
	    }
	    break;
	case 3:
	    for (ic[0] = icl[0]; ic[0] <= icu[0]; (ic[0])++)
	    for (ic[1] = icl[1]; ic[1] <= icu[1]; (ic[1])++)
	    for (ic[2] = icl[2]; ic[2] <= icu[2]; (ic[2])++)
	    {
		index = d_index(ic,top_gmax,dim);
		if (top_comp[index] != ctop_comp[dual_index])
		    continue;
		ave_mu += mu[index];
		n++;
	    }
	    break;
	}
	ave_mu /= (double)n;
	return ave_mu;
}	/* end computeDualMu */

void Incompress_Solver_Smooth_2D_Cartesian::updateComponent(void)
{
	int i,j,l,icoords[MAXD];
	int index;

	/*update the component of pressure on dual grid*/
	for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index]))
            {
                int cl[MAXD], cu[MAXD];
                for (l = 0; l < dim; ++l)
                    cl[l] = icoords[l] - offset[l];
                for (int m = 0; m < 2; ++m)
                for (int n = 0; n < 2; ++n)
                {
                    cu[0] = cl[0] + m;
                    cu[1] = cl[1] + n;
		    if (cu[0]<cimin||cu[0]>cimax||cu[1]<cjmin||cu[1]>cjmax)
			continue;
                    int index_tmp = d_index(cu,ctop_gmax,dim);
		    if (!ifluid_comp(ctop_comp[index_tmp]))
                    	ctop_comp[index_tmp] = top_comp[index];
                }
            }
        }

	/*Set rho for boundary layer on computational grid*/
	if (field->rho == NULL)
	    return;
	for (j = jmin; j <= jmax; ++j)
	for (i = imin; i <= imax; ++i)
	{
	    icoords[0] = i;
            icoords[1] = j;
            index = d_index(icoords,top_gmax,dim);
            if (!ifluid_comp(top_comp[index]))
	    {
		int cl[MAXD],cu[MAXD],indexl,indexu;
		boolean VelSet = NO;
		for (l = 0; l < dim && !VelSet; ++l)
		{
		    cl[l] = icoords[l]-1;
		    cu[l] = icoords[l]+1;
		    for (int n = -1; n < 2 && !VelSet; ++n)
		    {
			cl[(l+1)%dim] = cu[(l+1)%dim] = icoords[(l+1)%dim]+n;
			indexl = d_index(cl,top_gmax,dim);
                        if (ifluid_comp(top_comp[indexl]))
                        {
                            field->rho[index] = field->rho[indexl];
                            VelSet = YES;
                            continue;
                        }
                        indexu = d_index(cu,top_gmax,dim);
                        if (ifluid_comp(top_comp[indexu]))
                        {
                            field->rho[index] = field->rho[indexu];
                            VelSet = YES;
                            continue;
                        }
		    }
		}
	    }
	}
}	/* end updateComponent */

void Incompress_Solver_Smooth_3D_Cartesian::updateComponent(void)
{
        int i,j,k,l,icoords[MAXD];
        int index;

	/*update the component of pressure on dual grid*/
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (ifluid_comp(top_comp[index]))
            {
                int cl[MAXD], cu[MAXD];
                for (l = 0; l < dim; ++l)
                    cl[l] = icoords[l] - offset[l];
                for (int m = 0; m < 2; ++m)
                for (int n = 0; n < 2; ++n)
                for (int r = 0; r < 2; ++r)
                {
                    cu[0] = cl[0] + m;
                    cu[1] = cl[1] + n;
                    cu[2] = cl[2] + r;
		    if (cu[0]<cimin || cu[0]>cimax || cu[1]<cjmin 
			|| cu[1]>cjmax || cu[2]<ckmin || cu[2]>ckmax)
			continue;
                    int index_tmp = d_index(cu,ctop_gmax,dim);
		    if (!ifluid_comp(ctop_comp[index_tmp]))
		    	ctop_comp[index_tmp] = top_comp[index];
                }
            }
        }

        /*Set rho for boundary layer on computational grid*/
	if (field->rho == NULL)
	    return;
        for (k = kmin; k <= kmax; ++k)
        for (j = jmin; j <= jmax; ++j)
        for (i = imin; i <= imax; ++i)
        {
            icoords[0] = i;
            icoords[1] = j;
            icoords[2] = k;
            index = d_index(icoords,top_gmax,dim);
            if (!ifluid_comp(top_comp[index])&&!InsideSolid(icoords))
            {
                int cl[MAXD],cu[MAXD],indexl,indexu;
		boolean VelSet = NO;
		for (l = 0; l < dim && !VelSet; ++l)
		{
		    cl[l] = icoords[l]-1;
		    cu[l] = icoords[l]+1;
		    for (int n = -1; n < 2 && !VelSet; ++n)
                    for (int r = -1; r < 2 && !VelSet; ++r)
		    {
			cl[(l+1)%dim] = cu[(l+1)%dim] = icoords[(l+1)%dim]+n;
			cl[(l+2)%dim] = cu[(l+2)%dim] = icoords[(l+2)%dim]+r;
			indexl = d_index(cl,top_gmax,dim);
			if (ifluid_comp(top_comp[indexl]))
			{
			    field->rho[index] = field->rho[indexl];
			    VelSet = YES;
			    continue;
			}
			indexu = d_index(cu,top_gmax,dim);
			if (ifluid_comp(top_comp[indexu]))
                        {
                            field->rho[index] = field->rho[indexu];
                            VelSet = YES;
                            continue;
                        }
		    }
		}
            }
        }
}	/* end updateComponent */


