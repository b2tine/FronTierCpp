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
 * 		         iFbasic.cpp
 *******************************************************************/
#include "iFluid.h"
#include "iFturb.h"
#include "keps.h"

#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>
#include <stack>
#include <set>

//----------------------------------------------------------------
//		IF_RECTANGLE
//----------------------------------------------------------------

static double (*getStateVel[3])(POINTER) =
    {getStateXvel,getStateYvel,getStateZvel};

static double (*getStateOldVel[3])(POINTER) =
    {getStateOldXvel,getStateOldYvel,getStateOldZvel};

static double (*getStateGradPhi[3])(POINTER) =
    {getStateGradPhiX,getStateGradPhiY,getStateGradPhiZ};

IF_RECTANGLE::IF_RECTANGLE()
    : comp(-1)
{}

void IF_RECTANGLE::setCoords(
	double* m_coords,
	int dim)
{
	for (int i = 0; i < dim; ++i)
        coords[i] = m_coords[i];
}

std::vector<double> IF_RECTANGLE::getCoords()
{
    std::vector<double> m_coords(coords,coords+3);
    return m_coords;
}
//--------------------------------------------------------------------------
//               Incompress_Solver_Basis
//               Pure virtual class	
//--------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_Basis
//------------------------------------------------------------------------------
Incompress_Solver_Smooth_Basis::Incompress_Solver_Smooth_Basis(Front &front)
    : front(&front)
{
	skip_neumann_solver = 0;
}


//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup cell_center
//---------------------------------------------------------------
void Incompress_Solver_Smooth_Basis::initMesh(void)
{
	if (debugging("trace"))
            (void) printf("Entering initMesh()\n");
    
    iFparams = (IF_PARAMS*)front->extra1;
    
    FT_MakeGridIntfc(front);
	setDomain();

	int num_cells = 1;
	for (int i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	
    cell_center.clear();
	IF_RECTANGLE rectangle;
	cell_center.insert(cell_center.end(),num_cells,rectangle);

    int nfaces = 2*dim;
	for (int l = 0; l < nfaces; ++l)
    {
        ghost_data[l].clear();
        GHOST_COMPUTATION gc;
        ghost_data[l].insert(ghost_data[l].end(),num_cells,gc);
    }
	
	
    int index;
    double coords[MAXD];

	switch (dim)
	{
	case 2:
	    for (int j = 0; j <= top_gmax[1]; j++)
	    for (int i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
		    index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (int k = 0; k <= top_gmax[2]; k++)
	    for (int j = 0; j <= top_gmax[1]; j++)
	    for (int i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
	    	coords[2] = top_L[2] + top_h[2]*k;
		    index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(coords,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}

	setComponent();
    FT_FreeGridIntfc(front);

    if (debugging("trace"))
            (void) printf("Leaving initMesh()\n");
}

void Incompress_Solver_Smooth_Basis::setComponent(void)
{
	static POINTER state;
    double coords[MAXD];
	int size = (int)cell_center.size();
	double **vel = field->vel;
	double *pres = field->pres;
	
	for (int i = 0; i < size; i++)
	{
        cell_center[i].comp = getComponent(cell_center[i].icoords);
    }

    if(state == NULL)
        FT_ScalarMemoryAlloc((POINTER*)&state,front->sizest);

	for (int i = 0; i < size; i++)
	{
        if (cell_center[i].comp != -1 &&
            cell_center[i].comp != top_comp[i])
        {
            getRectangleCenter(i, coords);
            if (!FrontNearestIntfcState(front,coords,top_comp[i],(POINTER)state))
            {
                (void) printf("In setComponent()\n");
                (void) printf("FrontNearestIntfcState() failed\n");
                (void) printf("old_comp = %d new_comp = %d\n",
                        cell_center[i].comp,top_comp[i]);
                LOC(); clean_up(EXIT_FAILURE);
            }
    
            for (int l = 0; l < dim; ++l)
            {
                vel[l][i] = getStateVel[l](state);
            }
            pres[i] = getStatePres(state);
        }
        cell_center[i].comp = top_comp[i];
	}
}

void Incompress_Solver_Smooth_Basis::setIndexMap(void)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];
	int size = iupper - ilower;
	static int old_size;

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,
				top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,
				top_gmax[1]+1,top_gmax[2]+1,INT);
	    	break;
	    }
	    old_size = size;
	}
	else if (old_size < size)
	{
	    old_size = size;
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
		    ij_to_I[i][j] = -1;
	    
        for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
            ic = d_index2d(i,j,top_gmax);
            if (domain_status[ic] != TO_SOLVE) continue;

            if (cell_center[ic].comp != SOLID_COMP)
            {
                ij_to_I[i][j] = index + ilower;
                index++;
            }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
		    ijk_to_I[i][j][k] = -1;

	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
            ic = d_index3d(i,j,k,top_gmax);
            if (domain_status[ic] != TO_SOLVE) continue;

            if (cell_center[ic].comp != SOLID_COMP)
            {
                ijk_to_I[i][j][k] = index + ilower;
                index++;
            }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
}	/* end setIndexMap */

void Incompress_Solver_Smooth_Basis::getVelocity(double *p, double *U)
{
	int i;
	double **vel = field->vel;

	if (!FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[0],getStateXvel,
					&U[0],NULL))
	{
	    for (i = 0; i < dim; ++i) U[i] = 0.0;
	    return;
	}
	FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[1],getStateYvel,&U[1],
					NULL);
	if (dim == 3)
	    FT_IntrpStateVarAtCoords(front,NO_COMP,p,vel[2],getStateZvel,
					&U[2],NULL);
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void Incompress_Solver_Smooth_Basis::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int Incompress_Solver_Smooth_Basis::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void Incompress_Solver_Smooth_Basis::getRectangleCenter(
	int index, 
	double *coords)
{
	for (int i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

double Incompress_Solver_Smooth_Basis::getDistance(double *c0, double *c1)
{
	int i;
	double dist;
	dist = 0.0;
	for (i = 0; i < dim; ++i)
	    dist += sqr(c0[i] - c1[i]);
	dist = sqrt(dist);
	return dist;
}


// input : p[]
// output: q[]

void Incompress_Solver_Smooth_Basis::getNearestInterfacePoint(
	COMPONENT comp,
	double *p,          // mesh point
	double *q,          // intfc point
	double *nor,		// Normal vector
	double *kappa)		// curvature
{
	INTERFACE *intfc = front->interf;
	int i,j,dim = front->rect_grid->dim;
	double t[MAXD],mag_nor;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	POINT *pts[MAXD];
	BOND *bond;
	TRI *tri;
	static double **pts_nor,*pts_kappa;

	if (pts_kappa == NULL)
	{
	    FT_VectorMemoryAlloc((POINTER*)&pts_kappa,MAXD,FLOAT);
	    FT_MatrixMemoryAlloc((POINTER*)&pts_nor,MAXD,MAXD,FLOAT);
	}

	nearest_interface_point(p,comp,intfc,NO_BOUNDARIES,NULL,q,t,&hse,&hs);
	
    if (hse != NULL)
	{
	    switch (dim)
        {
            case 2:

            bond = Bond_of_hse(hse);
            pts[0] = bond->start;
            pts[1] = bond->end;
            for (i = 0; i < dim; ++i)
            {
                GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
                GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
            }
            *kappa = (1.0 - t[0])*pts_kappa[0] + t[0]*pts_kappa[1];
            for (i = 0; i < dim; ++i)
                nor[i] = (1.0 - t[0])*pts_nor[0][i] + t[0]*pts_nor[1][i];
            mag_nor = mag_vector(nor,dim);
            for (i = 0; i < dim; ++i)
                nor[i] /= mag_nor;
            break;

            
            case 3:

            tri = Tri_of_hse(hse);	
            for (i = 0; i < dim; ++i)
            {
                pts[i] = Point_of_tri(tri)[i];
                GetFrontCurvature(pts[i],hse,hs,&pts_kappa[i],front);
                GetFrontNormal(pts[i],hse,hs,pts_nor[i],front);
            }
            for (i = 0; i < dim; ++i)
            {
                *kappa = t[i]*pts_kappa[i];
                for (j = 0; j < dim; ++j)
                nor[j] = t[i]*pts_nor[i][j];
            }
            mag_nor = mag_vector(nor,dim);
            for (i = 0; i < dim; ++i)
                nor[i] /= mag_nor;
        }
	}
	else
	{
	    *kappa = 0.0;
	    nor[0] = 1.0;
	    for (i = 1; i < dim; ++i) nor[i] = 0.0;
	}
}	/* end getNearestInterfacePoint */

int Incompress_Solver_Smooth_Basis::getComponent(
	double *p)
{
	return I_ComponentAtCoords(p,front->interf);
}

int Incompress_Solver_Smooth_Basis::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
	return NO_COMP;
}

//TODO: should be SaveAsTecplot_rect_grid_and_interface()
void Incompress_Solver_Smooth_Basis::save(char *filename)
{
	
	INTERFACE *intfc    = front->interf;
		
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		clean_up(EXIT_FAILURE);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

void Incompress_Solver_Smooth_Basis::setDomain()
{
	static int current_size = 0;
	INTERFACE *grid_intfc;
	Table *T;

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
    field = iFparams->field;

	hmin = top_h[0];
	int size = top_gmax[0] + 1;
    for (int i = 1; i < dim; ++i)
	{
        if (top_h[i] < hmin) hmin = top_h[i];
        size *= (top_gmax[i] + 1);
	}

	switch (dim)
	{
    case 2:
	    if (size > current_size)
        {
            if (field != NULL)
            {
                FT_FreeThese(18,array,source,diff_coeff,field->mu,
                    field->rho,field->pres,field->phi,field->grad_phi,
                    field->q,field->div_U,field->vort,field->vel,
                    field->vel_star,field->prev_vel,field->grad_q,field->f_surf,
                    field->adv_term,field->adv_term_old,domain_status);
                
                if (debugging("field_var"))
                    FT_FreeThese(1,field->old_var);
                
                FT_FreeThese(1,field);
            }

            FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
            iFparams->field = field;

            FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&source,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&diff_coeff,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->mu,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->rho,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->pres,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->q,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->grad_q,2,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->phi,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->grad_phi,2,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->vel,2,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->vel_star,2,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->prev_vel,2,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->vort,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->div_U,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->f_surf,2,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->adv_term,2,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->adv_term_old,2,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&domain_status,size,INT);

            if (debugging("field_var"))
            {
                FT_MatrixMemoryAlloc((POINTER*)&field->old_var,2,size,sizeof(double));
            }

            current_size = size;
        }
	    
        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    
        break;
	
    case 3:
	    if (size > current_size)
        {
            if (field != NULL)
            {
                FT_FreeThese(18,array,source,diff_coeff,field->mu,
                    field->rho,field->pres,field->phi,field->grad_phi,
                    field->q,field->div_U,field->vel,field->vel_star,
                    field->prev_vel,field->vorticity,field->grad_q,field->f_surf,
                    field->adv_term,field->adv_term_old,domain_status);
                
                if (debugging("field_var"))
                    FT_FreeThese(1,field->old_var);
                
                FT_FreeThese(1,field);
            }

            FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(IF_FIELD));
            iFparams->field = field;

            FT_VectorMemoryAlloc((POINTER*)&array,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&source,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&diff_coeff,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->mu,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->rho,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->pres,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->q,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->grad_q,3,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->phi,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->grad_phi,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->vel,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->vel_star,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->prev_vel,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->vorticity,3,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&field->div_U,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->f_surf,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->adv_term,3,size,sizeof(double));
            FT_MatrixMemoryAlloc((POINTER*)&field->adv_term_old,3,size,sizeof(double));
            FT_VectorMemoryAlloc((POINTER*)&domain_status,size,INT);
            
            if (debugging("field_var"))
            {
                FT_MatrixMemoryAlloc((POINTER*)&field->old_var,3,size,sizeof(double));
            }
        
            current_size = size;
        }
	    
        imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    
        break;
	}
    
    
    if (iFparams->num_scheme.ellip_method == DOUBLE_ELLIP)
    {
        setDoubleDomain();
    }
}

void Incompress_Solver_Smooth_Basis::setGlobalIndex()
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		if (domain_status[ic] != TO_SOLVE) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == SOLID_COMP) continue;
		if (domain_status[ic] != TO_SOLVE) continue;
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
    iupper = n_dist[0];

    for (i = 1; i <= myid; i++)
    {
        ilower += n_dist[i-1];
        iupper += n_dist[i];
    }	
}

void Incompress_Solver_Smooth_Basis::printFrontInteriorStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[200];
	FILE *outfile;
	double **vel = field->vel;

	sprintf(filename,"%s/state.ts%s",out_name,
			right_flush(front->step,7));
#if defined(HAVE_MPI)
	if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(HAVE_MPI) */
	sprintf(filename,"%s-ifluid",filename);
	outfile = fopen(filename,"w");

	fluid_print_front_states(outfile,front);
	
	fprintf(outfile,"\nInterior ifluid states:\n");
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",field->rho[index]);
	        fprintf(outfile,"%24.18g\n",field->pres[index]);
	        fprintf(outfile,"%24.18g\n",field->phi[index]);
	        fprintf(outfile,"%24.18g\n",field->mu[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",vel[l][index]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",field->rho[index]);
	        fprintf(outfile,"%24.18g\n",field->pres[index]);
	        fprintf(outfile,"%24.18g\n",field->phi[index]);
	        fprintf(outfile,"%24.18g\n",field->mu[index]);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%24.18g\n",vel[l][index]);
	    }
	}
	fclose(outfile);
}

void Incompress_Solver_Smooth_Basis::readFrontInteriorStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	char fname[200];
	double **vel = field->vel;
	double speed;


	m_rho[0] = iFparams->rho1;		
	m_rho[1] = iFparams->rho2;		
	m_mu[0] = iFparams->mu1;		
	m_mu[1] = iFparams->mu2;		
	m_sigma = iFparams->surf_tension;
	m_comp[0] = iFparams->m_comp1;
	m_comp[1] = iFparams->m_comp2;
	m_smoothing_radius = iFparams->smoothing_radius;
	mu_min = HUGE;
	for (i = 0; i < 2; i++)
	{
	    if (ifluid_comp(m_comp[i]))
	    {
		    mu_min = std::min(mu_min,m_mu[i]);
		    rho_min = std::min(rho_min,m_rho[i]);
	    }
	}

	sprintf(fname,"%s-ifluid",restart_name);
	infile = fopen(fname,"r");
	
	/* Initialize states at the interface */
	fluid_read_front_states(infile,front);

	FT_MakeGridIntfc(front);
	setDomain();

	next_output_line_containing_string(infile,"Interior ifluid states:");

	max_speed = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    vmin[i] = HUGE;
	    vmax[i] = -HUGE;
        abs_vmax[i] = -HUGE;
	}

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
            index = d_index2d(i,j,top_gmax);
            fscanf(infile,"%lf",&field->rho[index]);
            fscanf(infile,"%lf",&field->pres[index]);
            fscanf(infile,"%lf",&field->phi[index]);
            fscanf(infile,"%lf",&field->mu[index]);
            
            speed = 0.0;
            for (l = 0; l < dim; ++l)
            {
                fscanf(infile,"%lf",&vel[l][index]);
                speed += sqr(vel[l][index]); 
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmin[l]));
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmax[l]));
            }

            speed = sqrt(speed);
            if (max_speed < speed)
            {
                max_speed = speed;
                icrds_max[0] = i;
                icrds_max[1] = j;
            }
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
            index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&field->rho[index]);
	    	fscanf(infile,"%lf",&field->pres[index]);
	    	fscanf(infile,"%lf",&field->phi[index]);
	    	fscanf(infile,"%lf",&field->mu[index]);
		
            speed = 0.0;
            for (l = 0; l < dim; ++l)
            {
                fscanf(infile,"%lf",&vel[l][index]);
                speed += sqr(vel[l][index]);
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmin[l]));
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmax[l]));
            }
    
            speed = sqrt(speed);
            if (max_speed < speed)
            {
                max_speed = speed;
                icrds_max[0] = i;
                icrds_max[1] = j;
                icrds_max[2] = k;
            }
	    }
	}

	fclose(infile);
	computeGradientQ();
	computeVelDivergence();
	copyMeshStates();
	setAdvectionDt();
}


//TODO: Need to be using vmax array to compute:
//
//      convect_max_dt =
//          1.0/(abs_vmax[0]/top_h[0] + abs_vmax[1]/top_h[1] + abs_vmax[2]/top_h[2]); 
//
//      viscous_max_dt =
//          1.0/((mu_max/rho)*(2/sqr(top_h[0]) + 2/sqr(top_h[1]) +2/sqr(top_h[2])));
//
//      From Fedkiw Paper:
//          "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
//
void Incompress_Solver_Smooth_Basis::setAdvectionDt()
{
    max_dt = HUGE;
    pp_global_max(&max_speed,1);
	pp_global_max(abs_vmax,dim);

    if (max_speed > MACH_EPS)
    {
        max_dt = hmin/max_speed;
    }
        
        /*
        //TODO: From Fedkiw paper -- producing very small dt's
        double mesh_val = 0.0;
        for (int i = 0; i < dim; ++i)
            mesh_val += abs_vmax[i]/top_h[i];
        max_dt = 1.0/mesh_val;
        */
 

    /*
    //TODO: input file option for setting iFparams->min_speed
    //          -- What/how is the min_speed used? 
    if (iFparams->min_speed > MACH_EPS) //if (iFparams->min_speed != 0.0)
    {
        max_dt = FT_Min(max_dt,hmin/iFparams->min_speed);
    }
    */

    //viscous time step restriction
    visc_max_dt = HUGE;
    pp_global_max(&mu_max,1);
    pp_global_max(&rho_min,1);

    if (mu_max > MACH_EPS)
    {
        //TODO: Is this correct for our scheme?
        double mesh_val = 0.0;
        for (int i = 0; i < dim; ++i)
        {
            mesh_val += 2.0/sqr(top_h[i]);
        }

        visc_max_dt = 1.0/((mu_max/rho_min)*mesh_val);
            //visc_max_dt = 0.5*hmin*hmin/mu_max;
    }

    max_dt = std::min(max_dt,visc_max_dt);
	
    //TODO: Is this correct??? See accum_dt and min_dt for projection step.
    min_dt = 0.0000001*sqr(hmin)/mu_min;
	
    if (debugging("ifluid_dt"))
	{
	    if (max_dt == HUGE)
	    	(void) printf("In Incompress_Solver_Smooth_Basis::setAdvectionDt(): \n"
                    "max_dt = HUGE min_dt = %24.18g visc_max_dt = %24.18g\n",
                    min_dt, visc_max_dt);
	    else
	    	(void) printf("In Incompress_Solver_Smooth_Basis::setAdvectionDt():\n"
                    "max_dt = %24.18g min_dt = %24.18g visc_max_dt = %24.18g\n",
                    max_dt, min_dt, visc_max_dt);
	}
}	/* end setAdvectionDt */


void Incompress_Solver_Smooth_Basis::initMovieVariables()
{
	boolean set_bound = NO;
	FILE *infile = fopen(InName(front),"r");
	char string[100];
	double var_max,var_min;

	if (CursorAfterStringOpt(infile,"Type y to set movie bounds:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                set_bound = YES;
        }

	    /* Begin hdf movies */
	switch (dim)
	{
	case 2:
        //TODO: make optional CursorAfterStringOpt()
	    CursorAfterString(infile,"Type y to make movie of pressure:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max pressure:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres",0,field->pres,getStatePres,
				var_max,var_min);
	    }
	    CursorAfterString(infile,"Type y to make movie of vorticity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max vorticity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"vort",0,field->vort,getStateVort,
				var_max,var_min);
	    }
	    CursorAfterString(infile,"Type y to make movie of velocity:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
		if (set_bound)
		{
		    CursorAfterString(infile,"Enter min and max velocity:");
                    fscanf(infile,"%lf %lf",&var_min,&var_max);
                    (void) printf("%f %f\n",var_min,var_max);
		}
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"xvel",0,field->vel[0],getStateXvel,
				var_max,var_min);
		FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"yvel",0,field->vel[1],getStateYvel,
				var_max,var_min);
	    }
	    if (CursorAfterStringOpt(infile,
		"Type y to make movie of viscosity:"))
	    {
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
	    	{
		    if (set_bound)
		    {
		        CursorAfterString(infile,
				"Enter min and max viscosity:");
                        fscanf(infile,"%lf %lf",&var_min,&var_max);
                        (void) printf("%f %f\n",var_min,var_max);
		    }
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"visc",0,field->mu,getStateMu,
				var_max,var_min);
		            //FT_AddVtkScalarMovieVariable(front,"VISC",field->mu);
	    	}
	    }
	    break;
	case 3:
        //TODO: make optional CursorAfterStringOpt()
	    CursorAfterString(infile,"Type y to make yz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-yz",0,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-y",0,field->vel[1],getStateYvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-yz-z",0,field->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterString(infile,"Type y to make xz cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xz",1,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-x",1,field->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xz-z",1,field->vel[2],getStateZvel,
				0.0,0.0);
		}
	    }
	    CursorAfterString(infile,"Type y to make xy cross section movie:");
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
	    {
	    	CursorAfterString(infile,"Type y to make movie of pressure:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"pres-xy",2,field->pres,getStatePres,
				0.0,0.0);
	    	CursorAfterString(infile,"Type y to make movie of velocity:");
            	fscanf(infile,"%s",string);
            	(void) printf("%s\n",string);
            	if (string[0] == 'Y' || string[0] == 'y')
		{
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-x",2,field->vel[0],getStateXvel,
				0.0,0.0);
		    FT_AddHdfMovieVariable(front,set_bound,YES,SOLID_COMP,
				"velo-xy-y",2,field->vel[1],getStateYvel,
				0.0,0.0);
		}
	    }
	}

    //TODO: make movies optional with input file
	
    if (dim != 1)
    {
        if (CursorAfterStringOpt(infile,
                    "Type y to make scalar pressure field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"PRESSURE",field->pres);
        }
        if (CursorAfterStringOpt(infile,
                    "Type y to make scalar phi field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"PHI",field->phi);
        }
        if (CursorAfterStringOpt(infile,
                    "Type y to make scalar intermediate velocity divergence field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"DIV_USTAR",field->div_U);
        }
        if (CursorAfterStringOpt(infile,
                    "Type y to make viscosity field movie:"))
        {
            fscanf(infile,"%s",string);
            (void)printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkScalarMovieVariable(front,"VISC",field->mu);
        }
	    if (CursorAfterStringOpt(infile,
                    "Type y to make vector velocity field movie:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'Y' || string[0] == 'y')
                FT_AddVtkVectorMovieVariable(front,"VELOCITY",field->vel);
        }

        //TODO: get vorticity plot for 2d vtk
        if (dim == 3)
        {
            if (CursorAfterStringOpt(infile,
                        "Type y to make vector vorticity field movie:"))
            {
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
                if (string[0] == 'Y' || string[0] == 'y')
                    FT_AddVtkVectorMovieVariable(front,"VORTICITY",field->vorticity);
            }
        }
    }

	fclose(infile);
}	/* end initMovieVariables */

//TODO: should move into 2d class
void Incompress_Solver_Smooth_Basis::computeSubgridModel(void)
{
        int i,j,k,index,index0,index1,index2,index3,index4,size;  
        double *u, *v;
        double ulx,urx,vlx,vrx;
        double uly,ury,vly,vry;
        double ux,uy,vx,vy;
        double s, *s11, *s12, *s22;
        double *ss11, *ss12, *ss22;
        double *tau00, *tau01, *tau10, *tau11;
        double *vel_u, *vel_v, *vel_uu, *vel_uv, *vel_vv;
        double sum_vel_u,sum_vel_v,sum_vel_uu,sum_vel_uv,sum_vel_vv;
        double sum_s11,sum_s12,sum_s22,sum_ss11,sum_ss12,sum_ss22,sum_s;
        double *ma11, *ma12, *la11, *la12, *la22;
        double *cs, *cs_ave, *deno, *nume, *co_coords_y;
        double coords[2];
        int    *r, num_r;
        int    ii,jj,iii,jjj;
        const int nn = pp_numnodes();
        num_r = (int)(((top_U[1]-top_L[1])/top_h[1])+1);
	double **vel = field->vel;

        size = (top_gmax[0]+1)*(top_gmax[1]+1);
        FT_VectorMemoryAlloc((POINTER*)&u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&s22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ss22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau00,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau01,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau10,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&tau11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_u,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_v,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uu,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_uv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&vel_vv,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&co_coords_y,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&ma12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la11,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la12,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&la22,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&r,size,sizeof(int));
        FT_VectorMemoryAlloc((POINTER*)&cs,size,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&cs_ave,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&deno,num_r,sizeof(double));
        FT_VectorMemoryAlloc((POINTER*)&nume,num_r,sizeof(double));

        for (j = 0; j <= top_gmax[1]; j++)
        for (i = 0; i <= top_gmax[0]; i++)
        {
            index = d_index2d(i,j,top_gmax);
            u[index] = vel[0][index];
            v[index] = vel[1][index];
            getRectangleCenter(index, coords);
            co_coords_y[index] = coords[1] + (top_h[1]/2.0);
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin-2; i <= imax+2; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            index1  = d_index2d(i-1,j,top_gmax);
            ulx = u[index1];
            uly = v[index1];
            index2  = d_index2d(i+1,j,top_gmax);
            urx = u[index2];
            ury = v[index2];
            index3  = d_index2d(i,j-1,top_gmax);
            vlx = u[index3];
            vly = v[index3];
            index4  = d_index2d(i,j+1,top_gmax);
            vrx = u[index4];
            vry = v[index4];

            ux = (urx - ulx) / (2.0*top_h[0]);
            uy = (ury - uly) / (2.0*top_h[1]);
            vx = (vrx - vlx) / (2.0*top_h[0]);
            vy = (vry - vly) / (2.0*top_h[1]);
            s11[index0] = ux;
            s12[index0] = (uy + vx)/2;
            s22[index0] = vy;
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                         + (2*(s12[index0]*s12[index0]))
                         + (s22[index0]*s22[index0])));
            ss11[index0] = s*s11[index0];
            ss12[index0] = s*s12[index0];
            ss22[index0] = s*s22[index0];
            vel_u[index0] = u[index0];
            vel_v[index0] = v[index0];  
            vel_uu[index0] = u[index0]*u[index0]; 
            vel_uv[index0] = u[index0]*v[index0];  
            vel_vv[index0] = v[index0]*v[index0];      
        }

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                sum_vel_u = sum_vel_v = 0.0;
                sum_vel_uu = sum_vel_uv = sum_vel_vv = 0.0;
                sum_s11 = sum_s12 = sum_s22 = 0.0;
                sum_ss11 = sum_ss12 = sum_ss22 = sum_s = 0.0;
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0  = d_index2d(iii,jjj,top_gmax);
                    sum_vel_u += vel_u[index0];
                    sum_vel_v += vel_v[index0];
                    sum_vel_uu += vel_uu[index0];
                    sum_vel_uv += vel_uv[index0];
                    sum_vel_vv += vel_vv[index0];
                    sum_s11 += s11[index0];
                    sum_s12 += s12[index0];
                    sum_s22 += s22[index0];
                    sum_ss11 += ss11[index0];
                    sum_ss12 += ss12[index0];
                    sum_ss22 += ss22[index0];
                    sum_s += sqrt(2*( (s11[index0]*s11[index0]) 
                                  + (2*(s12[index0]*s12[index0])) 
                                  + (s22[index0]*s22[index0])));
                } 
                ma11[index] = (2.0*top_h[1]*top_h[1]*(sum_ss11/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s11/4.0));
                ma12[index] = (2.0*top_h[1]*top_h[1]*(sum_ss12/4.0))
                        - (2.0*4*top_h[1]*top_h[1]*(sum_s/4.0)*(sum_s12/4.0));
                la11[index] = (sum_vel_uu/4.0)-((sum_vel_u/4.0)*
			(sum_vel_u/4.0));
                la12[index] = (sum_vel_uv/4.0)-((sum_vel_u/4.0)*
			(sum_vel_v/4.0));
                la12[index] = (sum_vel_vv/4.0)-((sum_vel_v/4.0)*
			(sum_vel_v/4.0));
                r[index] = (int)(co_coords_y[index]/(2*top_h[1]));
            }
        }

        for (k = 0; k < num_r; k++)
        {
            deno[k] = 0.0;
            nume[k] = 0.0;
        }

        for (k = 0; k < num_r; k++)
        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index0 = d_index2d(ii,jj,top_gmax);
                if(k == r[index0])
                {
                    deno[k] += (ma11[index0]*ma11[index0]) + 
				(ma12[index0]*ma12[index0]);
                    nume[k] += (((la11[index0]/2.0)-(la22[index0]/2.0))*
				ma11[index0]) + (la12[index0]*ma12[index0]);
                }
            }
        }

        pp_gsync();
        
        if (nn > 1)
        {
           for (k = 0; k < num_r; k++)
           {
              pp_global_sum(&deno[k],1L);
              pp_global_sum(&nume[k],1L);
           }
        }

        for (k = 0; k < num_r; k++)
        {
            if(deno[k] < 10e-16)
                cs_ave[k] = 0.0;
            else
                cs_ave[k] = nume[k]/deno[k];
        }

        for (j = jmin; j <= (jmax/2); j++)
        {
            jj = (2*j)-1;
            for (i = imin-1; i <= (imax/2)+1; i++)
            {
                ii = (2*i)-4;
                index = d_index2d(ii,jj,top_gmax);
                for(jjj = jj; jjj < jj+2; jjj++)
                for(iii = ii; iii < ii+2; iii++)
                {
                    index0 = d_index2d(iii,jjj,top_gmax);
                    cs[index0] = cs_ave[r[index]];
                }
            }
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin-1; i <= imax+1; i++)
        {
            index0  = d_index2d(i,j,top_gmax);
            s = sqrt(2*( (s11[index0]*s11[index0]) 
                          + (2*(s12[index0]*s12[index0]))
                          + (s22[index0]*s22[index0])));
            tau00[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s11[index0]/2.0)-(s22[index0]/2.0));
            tau01[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau10[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*s*(s12[index0]);
            tau11[index0] = - 2.0*cs[index0]*top_h[0]*top_h[0]*
                              s*((s22[index0]/2.0)-(s11[index0]/2.0));
        }

        for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            index0 = d_index2d(i,j,top_gmax);
            index1 = d_index2d(i-1,j,top_gmax);
            index2 = d_index2d(i+1,j,top_gmax);
            index3 = d_index2d(i,j-1,top_gmax);
            index4 = d_index2d(i,j+1,top_gmax);

            vel[0][index0] += -m_dt*(
                              ((tau00[index2]-tau00[index1])/(2.0*top_h[0])) + 
                                ((tau01[index4]-tau01[index3])/(2.0*top_h[1])));
            vel[1][index0] += -m_dt*(
                              ((tau10[index2]-tau10[index1])/(2.0*top_h[0])) + 
                              ((tau11[index4]-tau11[index3])/(2.0*top_h[1])));
        }
        FT_FreeThese(2,u,v);
        FT_FreeThese(4,tau00,tau01,tau10,tau11);
        FT_FreeThese(6,s11,s12,s22,ss11,ss12,ss22);
        FT_FreeThese(5,vel_u,vel_v,vel_uu,vel_uv,vel_vv);
        FT_FreeThese(11,co_coords_y,ma11,ma12,la11,la12,la22,r,cs,cs_ave,
					deno,nume);
}       /* end compSGS */

//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_2D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0;
	else if (phi > m_smoothing_radius)
	    return 1;
	else
	    return 1.0/2 + phi/(2*m_smoothing_radius) + 
		   1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_2D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1])))) 
		    /(16.0*top_h[0]*top_h[1]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 2D */

double Incompress_Solver_Smooth_2D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : 0.5*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_2D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_2D_Basis::sampleVelocity()
{
        int i,j,index;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *line = sample->sample_coords;
	char *out_name = front->out_name;
        double coords[MAXD];
        double var1,var2,var;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l = -1;
        static double lambda;
	char dirname[256];
	double **vel = field->vel;
	static char **sample_color;

	if (sample_color == NULL)
	{
	    FT_MatrixMemoryAlloc((POINTER*)&sample_color,10,20,sizeof(char));
	    sprintf(sample_color[0],"red");
	    sprintf(sample_color[1],"blue");
	    sprintf(sample_color[2],"green");
	    sprintf(sample_color[3],"violet");
	    sprintf(sample_color[4],"orange");
	    sprintf(sample_color[5],"yellow");
	    sprintf(sample_color[6],"pink");
	    sprintf(sample_color[7],"cyan");
	    sprintf(sample_color[8],"light-gray");
	    sprintf(sample_color[9],"dark-gray");
	}
	if (pp_numnodes() > 1)
	    return;
	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
	sprintf(dirname, "%s/samples/sample-%d",out_name,step);
	if (!create_directory(dirname,NO))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
        switch (sample_type[0])
        {
        case 'x':
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index2d(l,0,top_gmax);
                    getRectangleCenter(index, coords);
                } while(line[0] >= coords[0]);
                --l;
                index = d_index2d(l,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index2d(l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda = (line[0] - x1) / (x2 - line[0]);
            }
            i = l;
	    sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[0][index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = vel[0][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[1][index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = vel[1][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->pres[index];
                index = d_index2d(i+1,j,top_gmax);
                var2 = field->pres[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],var);
            }
            fclose(sfile);
            break;
        case 'y':
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index2d(0,l,top_gmax);
                    getRectangleCenter(index, coords);
                } while (line[0] >= coords[1]);
                --l;
                index = d_index2d(0,l,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index2d(0,l+1,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
               lambda = (line[0] - y1) / (y2 - line[0]);
            }
            j = l;
	    sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[0][index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = vel[0][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = vel[1][index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = vel[1][index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
	    sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index2d(i,j,top_gmax);
                var1 = field->pres[index];
                index = d_index2d(i,j+1,top_gmax);
                var2 = field->pres[index];
                var = (var1 + lambda*var2) / (1.0 + lambda);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],var);
            }
            fclose(sfile);
            break;
        }
	count++;
}	/* end sampleVelocity2d */

//TODO: factor into separate components -- surf tension, turbulence etc.
//
//TODO: put clock() calls on LES turb models
void Incompress_Solver_Smooth_2D_Basis::setSmoothedProperties(void)
{
	boolean status;
	int i,j,l,index,sign;
	COMPONENT comp;
	double t[MAXD],force[MAXD];
	double center[MAXD],point[MAXD],H,D;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	double **f_surf = field->f_surf;
	double *mu = field->mu;
    double *pres = field->pres;
	double *rho = field->rho;
	double dist;
	int range = (int)(m_smoothing_radius+1);
	boolean first = YES;

    //zero f_surf array
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
        for (l = 0; l < dim; ++l)
            f_surf[l][index] = 0.0;
    }
    
    if (iFparams->use_eddy_visc == YES &&
        iFparams->eddy_visc_model == BALDWIN_LOMAX)
    {
	    range = FT_Max(range,(int)(5*iFparams->ymax/top_h[0]));
    }

    KE_PARAMS* ke_params;
    double* mu_t;
    double* tke;

    if (iFparams->use_eddy_visc == YES &&
        iFparams->eddy_visc_model == KEPSILON)
    {
        ke_params = computeMuOfKepsModel();
        mu_t = ke_params->field->mu_t;
        tke = ke_params->field->k;
    }

    mu_max = 0.0;
    rho_min = HUGE;

	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
        mu[index] = 0.0;

	    comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;
 
        getRectangleCenter(index,center);
	    
        int icoords[MAXD];
        icoords[0] = i;
        icoords[1] = j;

        status = FT_FindNearestIntfcPointInRange(front,
                comp,center,NO_BOUNDARIES,point,t,&hse,&hs,range);

        if (status == YES && 
                ifluid_comp(positive_component(hs)) &&
                ifluid_comp(negative_component(hs)) &&
                positive_component(hs) != negative_component(hs))
        {
            sign = (comp == m_comp[0]) ? -1 : 1;
            D = smoothedDeltaFunction(center,point);
            H = smoothedStepFunction(center,point,sign);
            mu[index] = m_mu[0] + (m_mu[1]-m_mu[0])*H;
            rho[index] = m_rho[0] + (m_rho[1]-m_rho[0])*H; 
        
            if (m_sigma != 0.0 && D != 0.0)
            {
                for (l = 0; l < dim; ++l) force[l] = 0.0;

                surfaceTension(center,hse,hs,force,m_sigma);
                for (l = 0; l < dim; ++l)
                {
                    force[l] /= -rho[index];
                    f_surf[l][index] += force[l];
                }
            }
        }
        else if (iFparams->use_eddy_visc == YES)
	    {
            switch (comp)
            {
                case LIQUID_COMP1:
                    mu[index] = m_mu[0];
                    rho[index] = m_rho[0];
                    break;
                case LIQUID_COMP2:
                    mu[index] = m_mu[1];
                    rho[index] = m_rho[1];
                    break;
            }

            switch (iFparams->eddy_visc_model)
            {
            case BALDWIN_LOMAX:
                
                /*status = FT_FindNearestIntfcPointInRange(front,
                        comp,center,NO_BOUNDARIES,point,t,&hse,&hs,range);*/

                if (status == YES &&
                        (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                         wave_type(hs) == ELASTIC_BOUNDARY))
                {
                    dist = distance_between_positions(center,point,dim);
                    mu[index] += computeMuOfBaldwinLomax(icoords,dist,first);
                    first = NO;
                }
                break;
            case KEPSILON:
                mu[index] += mu_t[index];
                    //pres[index] += 2.0/3.0*tke[index];
                break;
            case VREMAN:
                mu[index] += computeMuOfVremanModel(icoords);
                break;
            case SMAGORINSKY:
                mu[index] += computeMuofSmagorinskyModel(icoords); 
                break;
            default:
                (void) printf("Unknown eddy viscosity model!\n");
                clean_up(ERROR);
            } 
        }
        else
	    {
            switch (comp)
            {
            case LIQUID_COMP1:
                mu[index] = m_mu[0];
                rho[index] = m_rho[0];
                break;
            case LIQUID_COMP2:
                mu[index] = m_mu[1];
                rho[index] = m_rho[1];
                break;
            }
	    }

        //For computing viscous flux time step restriction
        if (mu_max < mu[index]) mu_max = mu[index];
        if (rho[index] < rho_min) rho_min = rho[index];
	}

	FT_ParallelExchGridArrayBuffer(mu,front,NULL);
        //FT_ParallelExchGridArrayBuffer(pres,front,NULL);
	FT_ParallelExchGridArrayBuffer(rho,front,NULL);
	FT_ParallelExchGridVectorArrayBuffer(f_surf,front);
}	/* end setSmoothedProperties2d */

/*
void Incompress_Solver_Smooth_2D_Basis::precompute()
{
	for (int j = jmin; j <= jmax; j++)
    for (int i = imin; i <= imax; i++)
	{
	    int index  = d_index2d(i,j,top_gmax);
	    COMPONENT comp = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;
 
        int icoords[MAXD] = {0};
        icoords[0] = i;
        icoords[1] = j;

        icoordsPrecompute(icoords);
    }
}
*/

//-------------------------------------------------------------------------------
//               Incompress_Solver_Smooth_3D_Basis
//------------------------------------------------------------------------------

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunction(double phi)
{
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
	    return 0.0;
	else if (phi > m_smoothing_radius)
	    return 1.0;
	else
	    return 1.0/2.0 + phi/(2.0*m_smoothing_radius) + 
		   1.0/(2.0*PI)*sin(PI*phi/m_smoothing_radius);
}

double Incompress_Solver_Smooth_3D_Basis::getSmoothingFunctionD(double *center, double *point)
{
        if (fabs(center[0]-point[0]) < 2.0*top_h[0] && 
		fabs(center[1]-point[1]) < 2*top_h[1] && 
		fabs(center[2]-point[2]) < 2*top_h[2])
	    return ((1.0 + cos((PI*(center[0]-point[0]))/(2.0*top_h[0]))) 
		    *(1.0 + cos((PI*(center[1]-point[1]))/(2.0*top_h[1]))) 
		    *(1.0 + cos((PI*(center[2]-point[2]))/(2.0*top_h[2])))) 
		/(64.0*top_h[0]*top_h[1]*top_h[2]);
	else
            return 0.0;
}	/* end getSmoothingFunctionD in 3D */

double Incompress_Solver_Smooth_3D_Basis::smoothedDeltaFunction(double *p, double *center)
{
	int i;
	double len,d[MAXD],r,delta;

	for (i = 0; i < dim; ++i) d[i] = p[i] - center[i];
	len = mag_vector(d,dim);
	r = len/hmin/m_smoothing_radius;	
	delta = (fabs(r) > 1.0) ? 0.0 : 0.5*(1.0 + cos(PI*r));
	return delta;
}	/* end smoothedDeltaFunction */

double Incompress_Solver_Smooth_3D_Basis::smoothedStepFunction(double *p, double *center, int sign)
{
	int i;
	double dist,dp[MAXD],x,H;

	for (i = 0; i < dim; ++i) dp[i] = p[i] - center[i];
	dist = mag_vector(dp,dim);
	x = sign*dist/hmin/m_smoothing_radius;	
	if (x < -1.0) H = 0.0;
	else if (x > 1.0) H = 1.0;
	else H = 0.5*(1.0 + x) + 0.5/PI*sin(PI*x);
	return H;
}	/* end smoothedStepFunction */

void Incompress_Solver_Smooth_3D_Basis::sampleVelocity()
{
        int i,j,k,index;
        double coords[MAXD];
        double velo1,velo2,velo_tmp1,velo_tmp2,velo;
        FILE *sfile;
        char sname[100];
        static int count = 0;
        static int step = 0;
        static int l=-1,m=-1;
        static double lambda1,lambda2;
	SAMPLE *sample = front->sample;
	char *sample_type = sample->sample_type;
	double *sample_line = sample->sample_coords;
	char *out_name = front-> out_name;
	char dirname[256];
	double **vel = field->vel;
	static char **sample_color;
	int sample_dir;

	if (sample_color == NULL)
	{
	    FT_MatrixMemoryAlloc((POINTER*)&sample_color,10,20,sizeof(char));
	    sprintf(sample_color[0],"red");
	    sprintf(sample_color[1],"blue");
	    sprintf(sample_color[2],"green");
	    sprintf(sample_color[3],"violet");
	    sprintf(sample_color[4],"orange");
	    sprintf(sample_color[5],"yellow");
	    sprintf(sample_color[6],"pink");
	    sprintf(sample_color[7],"cyan");
	    sprintf(sample_color[8],"light-gray");
	    sprintf(sample_color[9],"dark-gray");
	}
	if (pp_numnodes() > 1)
	    return;
	if (front->step < sample->start_step || front->step > sample->end_step)
	    return;
	if ((front->step - sample->start_step)%sample->step_interval)
	    return;
        if (step != front->step)
        {
            step = front->step;
            count = 0;
        }
        switch (sample_type[0])
	{
	case 'x':
	    switch (sample_type[1])
	    {
	    case 'y':
		sample_dir = 2;
		break;
	    case 'z':
		sample_dir = 1;
		break;
	    }
	    break;
	case 'y':
	    switch (sample_type[1])
	    {
	    case 'x':
		sample_dir = 2;
		break;
	    case 'z':
		sample_dir = 0;
		break;
	    }
	    break;
	case 'z':
	    switch (sample_type[1])
	    {
	    case 'y':
		sample_dir = 0;
		break;
	    case 'x':
		sample_dir = 1;
		break;
	    }
	    break;
	case '0':
	    sample_dir = 0;
	    break;
	case '1':
	    sample_dir = 1;
	    break;
	case '2':
	    sample_dir = 2;
	    break;
	default:
	    (void) printf("Unknown sample type\n");
	    clean_up(ERROR);
	}

        sprintf(dirname, "%s/samples/sample-%d",out_name,step);
	if (!create_directory(dirname,NO))
        {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
        switch (sample_dir)
        {
	case 2:
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index3d(l,0,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[0]);
                --l;
                index = d_index3d(l,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index3d(l+1,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda1 = (sample_line[0] - x1) / (x2 - sample_line[0]);
            }
            if (m == -1)
            {
                double y1,y2;
                do
                {
                    ++m;
                    index = d_index3d(0,m,0,top_gmax);
                    getRectangleCenter(index,coords);
                }while(sample_line[1]>=coords[1]);
                --m;
                index = d_index3d(0,m,0,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index3d(0,m+1,0,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda2 = (sample_line[1] - y1)/(y2 - sample_line[1]);
            }
            i = l;
            j = m;
            sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[0][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j+1,k,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i+1,j+1,k,top_gmax);
                velo2 = vel[0][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
            }
            fclose(sfile);

            sprintf(sname,"%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[1][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j+1,k,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i+1,j+1,k,top_gmax);
                velo2 = vel[1][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
            }
            fclose(sfile);

            sprintf(sname,"%s/z-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[2][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j+1,k,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i+1,j+1,k,top_gmax);
                velo2 = vel[2][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
            }
            fclose(sfile);

            sprintf(sname,"%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = field->pres[index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j+1,k,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i+1,j+1,k,top_gmax);
                velo2 = field->pres[index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
            }
            fclose(sfile);

            sprintf(sname,"%s/mu-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (k = kmin; k <= kmax; ++k)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = field->mu[index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = field->mu[index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j+1,k,top_gmax);
                velo1 = field->mu[index];
                index = d_index3d(i+1,j+1,k,top_gmax);
                velo2 = field->mu[index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[2],velo);
            }
            fclose(sfile);
            printf("sample line: x = %20.14f, y = %20.14f\n",coords[0],
                coords[1]);

            break;
	case 1:
            if (l == -1)
            {
                double x1,x2;
                do
                {
                    ++l;
                    index = d_index3d(l,0,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[0]);
                --l;
                index = d_index3d(l,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x1 = coords[0];
                index = d_index3d(l+1,0,0,top_gmax);
                getRectangleCenter(index,coords);
                x2 = coords[0];
                lambda1 = (sample_line[0] - x1) / (x2 - sample_line[0]);
            }
            if (m == -1)
            {
                double z1,z2;
                do
                {
                    ++m;
                    index = d_index3d(0,0,m,top_gmax);
                    getRectangleCenter(index,coords);
                }while(sample_line[1]>=coords[2]);
                --m;
                index = d_index3d(0,0,m,top_gmax);
                getRectangleCenter(index,coords);
                z1 = coords[2];
                index = d_index3d(0,0,m+1,top_gmax);
                getRectangleCenter(index,coords);
                z2 = coords[2];
                lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
            }
            i = l;
            k = m;
            sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[0][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i+1,j,k+1,top_gmax);
                velo2 = vel[0][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[1][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i+1,j,k+1,top_gmax);
                velo2 = vel[1][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/z-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = vel[2][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i+1,j,k+1,top_gmax);
                velo2 = vel[2][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (j = jmin; j <= jmax; ++j)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i+1,j,k,top_gmax);
                velo2 = field->pres[index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i+1,j,k+1,top_gmax);
                velo2 = field->pres[index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[1],velo);
            }
            fclose(sfile);

            printf("sample line: x = %20.14f, z = %20.14f\n",coords[0],
                coords[2]);
            break;
	case 0:
            if (l == -1)
            {
                double y1,y2;
                do
                {
                    ++l;
                    index = d_index3d(0,l,0,top_gmax);
                    getRectangleCenter(index, coords);
                }while(sample_line[0]>=coords[1]);
                --l;
                index = d_index3d(0,l,0,top_gmax);
                getRectangleCenter(index,coords);
                y1 = coords[1];
                index = d_index3d(0,l+1,0,top_gmax);
                getRectangleCenter(index,coords);
                y2 = coords[1];
                lambda1 = (sample_line[0] - y1)/(y2 - sample_line[0]);
            }
            if (m == -1)
            {
                double z1,z2;
                do
                {
                    ++m;
                    index = d_index3d(0,0,m,top_gmax);
                    getRectangleCenter(index,coords);
                }while(sample_line[1]>=coords[2]);
                --m;
                index = d_index3d(0,0,m,top_gmax);
                getRectangleCenter(index,coords);
                z1 = coords[2];
                index = d_index3d(0,0,m+1,top_gmax);
                getRectangleCenter(index,coords);
                z2 = coords[2];
                lambda2 = (sample_line[1] - z1)/(z2 - sample_line[1]);
            }
            j = l;
            k = m;
            sprintf(sname, "%s/x-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i,j+1,k,top_gmax);
                velo2 = vel[0][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[0][index];
                index = d_index3d(i,j+1,k+1,top_gmax);
                velo2 = vel[0][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/y-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i,j+1,k,top_gmax);
                velo2 = vel[1][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[1][index];
                index = d_index3d(i,j+1,k+1,top_gmax);
                velo2 = vel[1][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/z-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
	    fprintf(sfile,"Next\n");
	    fprintf(sfile,"color=%s\n",sample_color[count]);
	    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i,j+1,k,top_gmax);
                velo2 = vel[2][index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = vel[2][index];
                index = d_index3d(i,j+1,k+1,top_gmax);
                velo2 = vel[2][index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);

            sprintf(sname, "%s/p-%d.xg",dirname,count);
            sfile = fopen(sname,"w");
		    fprintf(sfile,"Next\n");
		    fprintf(sfile,"color=%s\n",sample_color[count]);
		    fprintf(sfile,"thickness=1.5\n");
            for (i = imin; i <= imax; ++i)
            {
                index = d_index3d(i,j,k,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i,j+1,k,top_gmax);
                velo2 = field->pres[index];
                velo_tmp1 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                index = d_index3d(i,j,k+1,top_gmax);
                velo1 = field->pres[index];
                index = d_index3d(i,j+1,k+1,top_gmax);
                velo2 = field->pres[index];
                velo_tmp2 = (velo1 + lambda1*velo2)/(1.0 + lambda1);

                velo = (velo_tmp1 + lambda2*velo_tmp2)/(1.0 + lambda2);
                getRectangleCenter(index,coords);
                fprintf(sfile,"%20.14f   %20.14f\n",coords[0],velo);
            }
            fclose(sfile);

            printf("sample line: y = %20.14f, z = %20.14f\n",coords[1],
                coords[2]);
	    break;
        default:
            (void) printf("Incorrect input for sample velocity!\n");
            break;
        }
	count++;
}	/* end sampleVelocity in 3D */

void Incompress_Solver_Smooth_Basis::initSampleVelocity(char *in_name)
{
        FILE *infile;
	static SAMPLE *sample;
	char *sample_type;
	double *sample_line;

	infile = fopen(in_name,"r");
	FT_ScalarMemoryAlloc((POINTER*)&sample,sizeof(SAMPLE));
	sample_type = sample->sample_type;
	sample_line = sample->sample_coords;

	if (dim == 2)
	{
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf",sample_line);
            (void) printf(" %f\n",sample_line[0]);
	}
	else if (dim == 3)
        {
            CursorAfterString(infile,"Enter the sample line type:");
            fscanf(infile,"%s",sample_type);
            (void) printf(" %s\n",sample_type);
            CursorAfterString(infile,"Enter the sample line coordinate:");
            fscanf(infile,"%lf %lf",sample_line,sample_line+1);
            (void) printf(" %f %f\n",sample_line[0],sample_line[1]);
        }
        CursorAfterString(infile,"Enter the start step for sample: ");
        fscanf(infile,"%d",&sample->start_step);
        (void) printf("%d\n",sample->start_step);
        CursorAfterString(infile,"Enter the end step for sample: ");
        fscanf(infile,"%d",&sample->end_step);
        (void) printf("%d\n",sample->end_step);
        CursorAfterString(infile,"Enter the step interval for sample: ");
        fscanf(infile,"%d",&sample->step_interval);
        (void) printf("%d\n",sample->step_interval);
	front->sample = sample;
        fclose(infile);
}	/* end initSampleVelocity */

//TODO: factor into separate components -- surf tension, turbulence etc.
//
//TODO: put clock() calls on LES turb models
void Incompress_Solver_Smooth_3D_Basis::setSmoothedProperties(void)
{
	boolean status;
	int i,j,k,l,index,sign; 
	COMPONENT comp;
    double t[MAXD],force[MAXD];
	double center[MAXD],point[MAXD],H,D;
	HYPER_SURF_ELEMENT *hse;
    HYPER_SURF *hs;
	double **f_surf = field->f_surf;
	double *mu = field->mu;
    double *pres = field->pres;
	double *rho = field->rho;
	double dist;
	int range = (int)(m_smoothing_radius+1);
	boolean first = YES;

    //zero f_surf array
    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index = d_index3d(i,j,k,top_gmax);
        for (l = 0; l < dim; ++l)
            f_surf[l][index] = 0.0;
    }
    
    if (iFparams->use_eddy_visc == YES &&
        iFparams->eddy_visc_model == BALDWIN_LOMAX)
    {
	    range = FT_Max(range,(int)(5*iFparams->ymax/top_h[0]));
    }

    KE_PARAMS* ke_params;
    double* mu_t;
    double* tke;
    
    if (iFparams->use_eddy_visc == YES &&
        iFparams->eddy_visc_model == KEPSILON)
    {
        ke_params = computeMuOfKepsModel();
        mu_t = ke_params->field->mu_t;
        tke = ke_params->field->k;
    }

    mu_max = 0.0;
    rho_min = HUGE;

    for (k = kmin; k <= kmax; k++)
	for (j = jmin; j <= jmax; j++)
    for (i = imin; i <= imax; i++)
	{
	    index  = d_index3d(i,j,k,top_gmax);			
        mu[index] = 0.0;
	    
        comp  = cell_center[index].comp;
	    if (!ifluid_comp(comp)) continue;

	    getRectangleCenter(index,center);
        
        int icoords[MAXD];
        icoords[0] = i;
        icoords[1] = j;
        icoords[2] = k;

        status = FT_FindNearestIntfcPointInRange(front,
                comp,center,NO_BOUNDARIES,point,t,&hse,&hs,range);

        if (status == YES &&
                ifluid_comp(positive_component(hs)) &&
                ifluid_comp(negative_component(hs)) && 
                positive_component(hs) != negative_component(hs))
        {//SURFACE TENSION
            sign = (comp == m_comp[0]) ? -1 : 1;
            D = smoothedDeltaFunction(center,point);
            H = smoothedStepFunction(center,point,sign);
            mu[index] = m_mu[0] + (m_mu[1]-m_mu[0])*H;
            rho[index] = m_rho[0] + (m_rho[1]-m_rho[0])*H;

            if (m_sigma != 0.0 && D != 0.0)
            {
                for (l = 0; l < dim; ++l) force[l] = 0.0;
                surfaceTension(center,hse,hs,force,m_sigma);
                for (l = 0; l < dim; ++l)
                {
                    force[l] /= -rho[index];
                    f_surf[l][index] += force[l];
                }
            }
        }
        else if (iFparams->use_eddy_visc == YES)
        {//EDDY VISCOSITY
            switch (comp)
            {
                case LIQUID_COMP1:
                    mu[index] = m_mu[0];
                    rho[index] = m_rho[0];
                    break;
                case LIQUID_COMP2:
                    mu[index] = m_mu[1];
                    rho[index] = m_rho[1];
                    break;
            }

            switch (iFparams->eddy_visc_model)
            {
            case BALDWIN_LOMAX:
        
                /*status = FT_FindNearestIntfcPointInRange(front,
                        comp,center,NO_BOUNDARIES,point,t,&hse,&hs,range);*/

                if (status == YES &&
                        (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY ||
                         wave_type(hs) == ELASTIC_BOUNDARY))
                {
                    dist = distance_between_positions(center,point,dim);
                    mu[index] += computeMuOfBaldwinLomax(icoords,dist,first);
                    first = NO;
                }
                break;
            case KEPSILON:
                mu[index] += mu_t[index];
                break;
            case VREMAN:
                mu[index] += computeMuOfVremanModel(icoords);
                break;
            case SMAGORINSKY:
                mu[index] += computeMuofSmagorinskyModel(icoords);
                break;
            default:
                (void) printf("Unknown eddy viscosity model!\n");
                clean_up(ERROR);
            } 
        }
	    else
	    {
            switch (comp)
            {
            case LIQUID_COMP1:
                mu[index] = m_mu[0];
                rho[index] = m_rho[0];
                break;
            case LIQUID_COMP2:
                mu[index] = m_mu[1];
                rho[index] = m_rho[1];
                break;
            }
	    }

        //For computing viscous flux time step restriction
        if (mu_max < mu[index]) mu_max = mu[index];
        if (rho[index] < rho_min) rho_min = rho[index];
	}

	FT_ParallelExchGridArrayBuffer(mu,front,NULL);
	FT_ParallelExchGridArrayBuffer(rho,front,NULL);
	FT_ParallelExchGridVectorArrayBuffer(f_surf,front);
}	/* end setSmoothedProperties in 3D */

// Flux of Riemann solution of Burgers equation u_t + uu_x = 0
double burger_flux(	
	double ul,
	double um,
	double ur)
{
	double u_Rl,u_Rr;
	if (ul < um)
	{
	    if (ul > 0.0) u_Rl = ul;
	    else if (um < 0.0) u_Rl = um;
	    else u_Rl = 0.0;
	}
	else
	{
	    if (ul + um > 0.0) u_Rl = ul;
	    else u_Rl = um;
	}

	if (um < ur)
	{
	    if (um > 0.0) u_Rr = um;
	    else if (ur < 0.0) u_Rr = ur;
	    else u_Rr = 0.0;
	}
	else
	{
	    if (um + ur > 0.0) u_Rr = um;
	    else u_Rr = ur;
	}
	return 0.5*(u_Rr*u_Rr - u_Rl*u_Rl);
}	/* end flux */

// Flux of Riemann solution of linear equation u_t + au_x = 0

double linear_flux(	
	double a,
	double ul,
	double um,
	double ur)
{
	if (a > 0.0)
	    return a*(um - ul);
	else
	    return a*(ur - um);
}	/* end net_uwind_flux */


void Incompress_Solver_Smooth_Basis::paintAllGridPoint(int status)
{
	int i,j,k,ic;	
	switch(dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		domain_status[ic] = status;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		domain_status[ic] = status;
	    }
	    break;
	}
}	/* end paintAllGridPoint */

void Incompress_Solver_Smooth_Basis::paintSolvedGridPoint()
{
	int i,j,k,ic;
	//following is for debugging	
	static int round = 0;
	round ++;
	static std::vector<std::vector<int> > map(jmax-jmin+1,
					std::vector<int>(imax-imin+1,0));
	switch(dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (domain_status[ic] == TO_SOLVE) {
		    domain_status[ic] = SOLVED;
		    map[j-jmin][i-imin] = round;
		}
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (domain_status[ic] == TO_SOLVE)
		    domain_status[ic] = SOLVED;
	    }
	    break;
	}
}	/* end paintSolvedGridPoin */

void Incompress_Solver_Smooth_Basis::paintConnectedRegion(
	int ic, 
	int color) 
{
	//depth first search (DFS) to find connected region
	//from seed ic and paint it with given color
	std::stack<int> stk;
	int next_index, next_icrds[MAXD];
	int dim = FT_Dimension();
	GRID_DIRECTION  dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	int smin[] = {imin,jmin,kmin};
        int smax[] = {imax,jmax,kmax};
	int top_gmin[] = {0, 0, 0};
        DOMAIN_METHOD paint_method = BY_CROSSING;
	int ipn[MAXD];

	if (domain_status[ic] != NOT_SOLVED)
	    return;
	int stk_count = 0;
	stk.push(ic);
        while (!stk.empty()) 
	{
            next_index = stk.top();
            domain_status[next_index] = TO_SOLVE;
            color_map[next_index] = color;
            stk.pop();
            invert_to_icoords(next_index,next_icrds,top_gmax,dim);
            for (int idir = 0; idir < dim*2; ++idir) 
	    {
                if (nextConnectedPoint(next_icrds,dir[idir],ipn,paint_method,
                    		top_gmin,top_gmax))
                {
                    stk.push(d_index(ipn,top_gmax,dim));
                }
            }
       	}
}	/* end paintConnectedRegion */

boolean Incompress_Solver_Smooth_Basis::paintToSolveGridPoint2(int c) 
{
	int i,j,k,ic;
	//fetch a new color and remove it
	//note that the first one is always
	//zero for solid color

        switch(dim)
        {
        case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                ic = d_index2d(i,j,top_gmax);
		if (color_map[ic] == c)
                    domain_status[ic] = TO_SOLVE;
            }
            break;
        case 3:
            for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
            {
                ic = d_index3d(i,j,k,top_gmax);
		if (color_map[ic] == c)
                    domain_status[ic] = TO_SOLVE;
            }
            break;
        }
	return YES;
}

boolean Incompress_Solver_Smooth_Basis::paintToSolveGridPoint()
{
	//test the new algorithm
	//comment it out to use old method

	GRID_DIRECTION  dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
	int idir,i,j,k,n,ic,ic1,*ip,ipn[MAXD];
	static int **ips,ip_size;
	boolean seed_found;
	int paint_method;
	int nv0,nv1;

	paint_method = BY_CROSSING;

	start_clock("paint_color");
	smin[0] = imin;	smin[1] = jmin;	smin[2] = kmin;
	smax[0] = imax;	smax[1] = jmax;	smax[2] = kmax;
	if (ips == NULL)
	{
	    ip_size = 1;
	    for (i = 0; i < dim; ++i)
	    {
	    	ip_size *= (smax[i] - smin[i] + 1); 
	    }
	    FT_MatrixMemoryAlloc((POINTER*)&ips,ip_size,dim,sizeof(int));
	}

	seed_found = NO;
	n = 0;
	switch (dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
            {
                for (i = imin; i <= imax; i++)
                {
                    ic = d_index2d(i,j,top_gmax);
                    if (domain_status[ic] == NOT_SOLVED)
                    {
			ips[n][0] = i;
			ips[n][1] = j;
			ip = ips[n];
                        domain_status[ic] = TO_SOLVE;
                        seed_found = YES;
			n++;
                        break;
                    }
                }
                if (seed_found) break;
            }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    {
	    	for (j = jmin; j <= jmax; j++)
	    	{
		    for (i = imin; i <= imax; i++)
		    {
		    	ic = d_index3d(i,j,k,top_gmax);	
		    	if (domain_status[ic] == NOT_SOLVED)
		    	{
			    ips[n][0] = i;
			    ips[n][1] = j;
			    ips[n][2] = k;
			    ip = ips[n];
			    domain_status[ic] = TO_SOLVE;
			    seed_found = YES;
			    n++;
			    break;
		    	}
		    }
		    if (seed_found) break;
	    	}
	    	if (seed_found) break;
	    }
	    break;
	}
	seed_found = pp_max_status(seed_found);

	if (!seed_found)
	{
	    stop_clock("paint_color");
	    return NO;
	}

	/* Start traversing through connected neighbors */
	std::stack<int> stk;
	int next_index, next_icrds[MAXD];
	stk.push(ic);
	while (!stk.empty()) {
	    next_index = stk.top();
	    domain_status[next_index] = TO_SOLVE;
	    stk.pop();
	    invert_to_icoords(next_index,next_icrds,top_gmax,3);
	    for (idir = 0; idir < dim*2; ++idir) {
		if (nextConnectedPoint(next_icrds,dir[idir],ipn,paint_method,
			smin,smax))
		{
		    stk.push(d_index(ipn,top_gmax,dim));
		}
	    }
	}
	FT_ParallelExchGridIntArrayBuffer(domain_status,front);

	stop_clock("paint_color");
	return YES;
}	/* end paintToSolveGridPoint */

boolean Incompress_Solver_Smooth_Basis::nextConnectedPoint(
	int *ip,
	GRID_DIRECTION dir,
	int *ipn,
	int paint_method,
	int *smin,
	int *smax)
{
	int id,idn;
	double crx_coords[MAXD];
	HYPER_SURF *hs;
	POINTER intfc_state;

	if (!next_ip_in_dir(ip,dir,ipn,smin,smax))
	    return NO;
	idn = d_index(ipn,top_gmax,dim);
	if (domain_status[idn] != NOT_SOLVED)
	    return NO;

	switch (paint_method)
	{
	case BY_COMPONENT:
	    id = d_index(ip,top_gmax,dim); 
	    idn = d_index(ipn,top_gmax,dim); 
	    if (top_comp[id] == top_comp[idn]) 
		return YES;
	    else 
		return NO;
	case BY_CROSSING:
	    id = d_index(ip,top_gmax,dim); 
	    if ((*findStateAtCrossing)(front,ip,dir,top_comp[id],
                                &intfc_state,&hs,crx_coords))
		return NO;
	    else
		return YES;
	case BY_WALL:
	    idn = d_index(ipn,top_gmax,dim); 
	    if (top_comp[idn] == SOLID_COMP) 
		return NO;
	    else 
		return YES;
	}
}	/* end connectedPoint */

void Incompress_Solver_Smooth_Basis::checkVelocityDiv(
	const char *mesg)
{
	double **vel = field->vel;
	double div_tmp,denom;
    double L1_norm,L2_norm,Li_norm;
    double max_div;	
    int i,j,k,l,index,N;
	int icoords[MAXD], icoords_max[MAXD];

	(void) printf("\nCheck divergence at %s\n",mesg);
        L1_norm = L2_norm = Li_norm = max_div = 0.0;
	switch (dim)
	{
	case 2:
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
                icoords[1] = j;
		index = d_index2d(i,j,top_gmax);
		if (!ifluid_comp(top_comp[index]))
		    continue;
		div_tmp = computeFieldPointDiv(icoords,vel);
		if (Li_norm < fabs(div_tmp)) Li_norm = fabs(div_tmp);
                L1_norm += fabs(div_tmp);
                L2_norm += sqr(div_tmp);
                if (max_div < fabs(div_tmp))
                {
                    max_div = fabs(div_tmp);
                    icoords_max[0] = i;
                    icoords_max[1] = j;
                }
	    }
            N = (imax - imin + 1)*(jmax - jmin + 1);
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
            for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {
		icoords[0] = i;
                icoords[1] = j;
                icoords[2] = k;
		index = d_index3d(i,j,k,top_gmax);
		if (!ifluid_comp(top_comp[index]))
		    continue;
		div_tmp = computeFieldPointDiv(icoords,vel);
		if (Li_norm < fabs(div_tmp)) Li_norm = fabs(div_tmp);
                L1_norm += fabs(div_tmp);
                L2_norm += sqr(div_tmp);
                if (max_div < fabs(div_tmp))
                {
                    max_div = fabs(div_tmp);
                    icoords_max[0] = i;
                    icoords_max[1] = j;
                    icoords_max[2] = k;
                }
	    }
            N = (imax - imin + 1)*(jmax - jmin + 1)*(kmax - kmin + 1);
	    break;
	}
        L1_norm /= N;
        L2_norm = sqrt(L2_norm)/N;

	denom = 0.0;
	for (l = 0; l < dim; ++l)
	    denom += fabs((vmax[l] - vmin[l])/(top_U[l] - top_L[l]));
	if (denom == 0.0) denom = 1.0;
	(void) printf("Absolute: L1 = %5.3g  L2 = %5.3g  Li =  %5.3g\n",
			L1_norm,L2_norm,Li_norm);
	(void) printf("Relative: L1 = %5.3g  L2 = %5.3g  Li =  %5.3g\n",
			L1_norm/denom,L2_norm/denom,Li_norm/denom);
        (void) printf("Max divergence: %6.4g occured at: ");
	for (l = 0; l < dim; ++l)
            printf("%d ",icoords_max[l]);
        printf("\n");
        fflush(stdout);
}	/* end checkVelocityDiv */

double Incompress_Solver_Smooth_Basis::computeFieldPointDiv(
        int *icoords,
        double **field)
{
        switch (iFparams->num_scheme.ellip_method)
        {
        case SIMPLE_ELLIP:
            return computeFieldPointDivSimple(icoords,field);
        case DOUBLE_ELLIP:
            return computeFieldPointDivDouble(icoords,field);
        default:
            printf("Elliptic Method Not Implemented\n");
            clean_up(1);
        }
}       /* end computeFieldPointDiv */

double Incompress_Solver_Smooth_Basis::computeFieldPointDivSimple(
        int *icoords,
        double **field_array)
{
	int icnb[MAXD], icnb_opp1[MAXD], icnb_opp2[MAXD];
        int index,index_nb,index_oppnb1,index_oppnb2;
        COMPONENT comp;
	double u_edge[3][2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	int status;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int idir,nb;
	double u0,u_ref;

	index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];

        if (!ifluid_comp(comp)) return 0.0;

        for (idir = 0; idir < dim; idir++)
        {
            u0 = field_array[idir][index];
            for (int j = 0; j < dim; ++j)
            {
                icnb[j] = icoords[j];
                icnb_opp1[j] = icoords[j];
                icnb_opp2[j] = icoords[j];
            }

            for (nb = 0; nb < 2; nb++)
            {
                u_edge[idir][nb] = 0.0;
                
                icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                index_nb = d_index(icnb,top_gmax,dim);
                
                status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);
                
                if (status == NO_PDE_BOUNDARY)
                {
                    //NOTE: Includes ELASTIC_BOUNDARY when using af_findcrossing
                    u_edge[idir][nb] = field_array[idir][index_nb];
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
                    //OUTLET
                    u_edge[idir][nb] = getStateVel[idir](intfc_state);
                }
                else if (status == CONST_V_PDE_BOUNDARY)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        //INLET
                        u_edge[idir][nb] = getStateVel[idir](intfc_state);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                            wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                    {
                        //Use one sided 3pt derivative
                        icnb_opp1[idir] = (nb == 0) ? icoords[idir] + 1 : icoords[idir] - 1;
                        index_oppnb1 = d_index(icnb_opp1,top_gmax,dim);
                        double u_oppnb1 = field_array[idir][index_oppnb1];
                        
                        icnb_opp2[idir] = (nb == 0) ? icoords[idir] + 2 : icoords[idir] - 2;
                        index_oppnb2 = d_index(icnb_opp2,top_gmax,dim);
                        double u_oppnb2 = field_array[idir][index_oppnb2];

                        u_edge[idir][nb] = 3.0*u0 - 3.0*u_oppnb1 + u_oppnb2;
                    }
                }
            }
        }

    double div = 0.0;
    for (int i = 0; i < dim; ++i)
        div += 0.5*(u_edge[i][1] - u_edge[i][0])/top_h[i];

    return div;
}       /* end computeFieldPointDivSimple */

double Incompress_Solver_Smooth_Basis::computeFieldPointDivDouble(
        int *icoords,
        double **field)
{
	int icnb[MAXD];
        int i,j,index,index_nb;
        COMPONENT comp;
	double div,u_edge[3][2];
        double crx_coords[MAXD];
        POINTER intfc_state;
        HYPER_SURF *hs;
	int status;
        GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int idir,nb;
	double u0,u_ref;

	index = d_index(icoords,top_gmax,dim);
        comp = top_comp[index];

        if (iFparams->num_scheme.ellip_method != DOUBLE_ELLIP &&
            !ifluid_comp(comp)) return 0.0;

        for (idir = 0; idir < dim; idir++)
        {
            u0 = field[idir][index];
            for (j = 0; j < dim; ++j)
                icnb[j] = icoords[j];
            for (nb = 0; nb < 2; nb++)
            {
                u_edge[idir][nb] = 0.0;
                icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                index_nb = d_index(icnb,top_gmax,dim);
                status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                                comp,&intfc_state,&hs,crx_coords);
                if (status == NO_PDE_BOUNDARY)
                {
                    u_edge[idir][nb] = field[idir][index_nb];
                }
                else if (status == CONST_P_PDE_BOUNDARY)
                {
                    if (iFparams->num_scheme.ellip_method == DOUBLE_ELLIP)
                        return 0.0;
                    else
                        u_edge[idir][nb] = u0; 
                }
                else if (status == CONST_V_PDE_BOUNDARY)
                {
                    if (wave_type(hs) == DIRICHLET_BOUNDARY)
                    {
                        if (iFparams->num_scheme.ellip_method == DOUBLE_ELLIP)
                            return 0.0;
                        u_edge[idir][nb] = getStateVel[idir](intfc_state);
                    }
                    else if (wave_type(hs) == NEUMANN_BOUNDARY)
                    {
                        if (iFparams->num_scheme.ellip_method == DOUBLE_ELLIP)
                            u_edge[idir][nb] = -u0; 
                        else
                            u_edge[idir][nb] = u0;
                    }
                    else
                    {
                        //u_edge[idir][nb] = getStateVel[idir](intfc_state);
                        u_edge[idir][nb] = u0;
                    }
                }
            }
        }

	div = 0.0;
	for (i = 0; i < dim; ++i)
	    div += 0.5*(u_edge[i][1] - u_edge[i][0])/top_h[i];
        return div;
}       /* end computeFieldPointDivDouble */

//TODO: function has diverged from the original purpose and now is concerned
//      with computing the gradient of phi -- should rename function.      
void Incompress_Solver_Smooth_Basis::computeFieldPointGrad(
        int *icoords,
        double *field_array,
        double *grad_field)
{
    int index, index_nb, index_nb_opp;
    int icnb[MAXD], icnb_opp[MAXD];
    COMPONENT comp;
    int i,j,idir,nb;
	double p_edge[3][2],p0;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int status;
	boolean refl_side[2];
    bool flowthrough[3] = {false,false,false};

    double** vel_star = field->vel_star;
    double** prev_vel= field->prev_vel;
    double* rho = field->rho;
    double* mu = field->mu;

	index = d_index(icoords,top_gmax,dim);
    comp = top_comp[index];

    if (!ifluid_comp(comp))
    {
        for (i = 0; i < dim; ++i)
            grad_field[i] = 0.0;
        return;
    }
	
    p0 = field_array[index];

	for (idir = 0; idir < dim; idir++)
	{
	    for (j = 0; j < dim; ++j)
        {
	    	icnb[j] = icoords[j];
	    	icnb_opp[j] = icoords[j];
        }

	    for (nb = 0; nb < 2; nb++)
	    {
            refl_side[nb] = NO;
	    	icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
	    	icnb_opp[idir] = (nb == 0) ? icoords[idir] + 1 : icoords[idir] - 1;
	    	
            index_nb = d_index(icnb,top_gmax,dim);
            index_nb_opp = d_index(icnb_opp,top_gmax,dim);
	
            status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                    comp,&intfc_state,&hs,crx_coords);

	    	if (status == NO_PDE_BOUNDARY)
            {
		        p_edge[idir][nb] = field_array[index_nb];
            }
            else if(wave_type(hs) == DIRICHLET_BOUNDARY)
            {
	    	    if (status == CONST_P_PDE_BOUNDARY)
                {
                    //OUTLET
                    p_edge[idir][nb] = getStatePhi(intfc_state);

                    /*
                    //////////////////////////////////////////////////////////////////////
                    //      n dot grad(phi^n+1) = 0.5*mu * n dot (grad^2(u^{*} + u^{n})
                    //      t dot grad(phi^n+1) = 0.5*mu * t dot (grad^2(u^{*} + u^{n})

                    std::vector<double> vector_laplacian(dim,0.0);
                    for (int ii = 0; ii < dim; ++ii)
                    {
                        //compute laplacian i-th component of vel and prev_vel
                        for (int m = 0; m < dim; ++m)
                        {
                            double laplacian = vel_star[m][index_nb_opp]
                                - 2.0*vel_star[m][index] + getStateVel[m](intfc_state);
                            laplacian /= top_h[m]*top_h[m];

                            double prev_laplacian = prev_vel[m][index_nb_opp]
                                - 2.0*prev_vel[m][index] + getStateOldVel[m](intfc_state);
                            prev_laplacian /= top_h[m]*top_h[m];

                            vector_laplacian[ii] += laplacian + prev_laplacian;
                        }
                    }

                    double nor[MAXD];
                    FT_NormalAtGridCrossing(front,icoords,
                            dir[idir][nb],comp,nor,&hs,crx_coords);

                    double laplace_n = 0.0;
                    for (int ii = 0; ii < dim; ++ii)
                        laplace_n += vector_laplacian[ii]*nor[ii];
                    
                    //TODO: tangential components

	                grad_field[idir] = 0.5*mu[index]*laplace_n;
                    flowthrough[idir] = true;
                    continue;
                    //////////////////////////////////////////////////////////////////////
                    */
                }
                else 
                {
                    //INLET
                    p_edge[idir][nb] = p0;//conforms with do-nothing boundary
                }
            }
            else if (!is_bdry_hs(hs) && 
                     (wave_type(hs) == NEUMANN_BOUNDARY ||
                      wave_type(hs) == MOVABLE_BODY_BOUNDARY))
            {
                //grad(phi) dot normal = 0
                int icoords_ghost[MAXD];
                for (int m = 0; m < dim; ++m)
                    icoords_ghost[m] = icoords[m];
                
                icoords_ghost[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
                
                double coords_ghost[MAXD];
                double coords_reflect[MAXD];
                
                ////////////////////////////////////////////////////////////////////////
                ///  matches Incompress_Solver_Smooth_Basis::setSlipBoundaryGNOR()  ///
                //////////////////////////////////////////////////////////////////////

                for (int m = 0; m < dim; ++m)
                {
                    coords_ghost[m] = top_L[m] + icoords_ghost[m]*top_h[m];
                    coords_reflect[m] = coords_ghost[m];
                }

                double nor[MAXD];
                FT_NormalAtGridCrossing(front,icoords,
                        dir[idir][nb],comp,nor,&hs,crx_coords);
                        
                //Reflect the ghost point through intfc-mirror at crossing.
                //first reflect across the grid line containing intfc crossing,
                coords_reflect[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];
                
                //Reflect the displacement vector across the line
                //containing the intfc normal vector
                double v[MAXD];
                double vn = 0.0;

                for (int m = 0; m < dim; ++m)
                {
                    v[m] = coords_reflect[m] - crx_coords[m];
                    vn += v[m]*nor[m];
                }

                for (int m = 0; m < dim; ++m)
                    v[m] = 2.0*vn*nor[m] - v[m];

                //The desired reflected point
                for (int m = 0; m < dim; ++m)
                    coords_reflect[m] = crx_coords[m] + v[m];
                ///////////////////////////////////////////////////////////////////////

                //Interpolate phi at the reflected point,
                double phi_reflect;
                FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field_array,
                        getStatePhi,&phi_reflect,&field_array[index]);
                
                //FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field_array,
                //        getStatePhi,&phi_reflect,nullptr);//default_ans uses nearest intfc state

                p_edge[idir][nb] = phi_reflect;
            }
            else if (is_bdry_hs(hs) && wave_type(hs) == NEUMANN_BOUNDARY)
            {
                p_edge[idir][nb] = p0;
            }
            else
            {
                printf("ERROR in computeFieldPointGrad() : Unknown Boundary Type!\n");
                LOC(); clean_up(EXIT_FAILURE);
            }
	    
        }
	}

	for (i = 0; i < dim; ++i)
    {
        //if (flowthrough[i]) continue;
	    grad_field[i] = 0.5*(p_edge[i][1] - p_edge[i][0])/top_h[i];
    }
}      /* end computeFieldPointGrad */

void Incompress_Solver_Smooth_Basis::computeFieldPointGradQ(
        int *icoords,
        double *field_array,
        double *grad_field)
{
    COMPONENT comp;
    int i,j,idir,nb;
    double crx_coords[MAXD];
    POINTER intfc_state;
    HYPER_SURF *hs;
    GRID_DIRECTION dir[3][2] = {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};
	int status;
	boolean refl_side[2];

	int icnb[MAXD], icnb_opp[MAXD];
    int index_nb,index_nb_opp;
    int icnb_opp1[3], icnb_opp2[3];
    int index_oppnb1, index_oppnb2;


    for (i = 0; i < dim; ++i)
        grad_field[i] = 0.0;

    if (iFparams->num_scheme.projc_method == PMIII ||
        iFparams->num_scheme.projc_method == SIMPLE) return;

    int index = d_index(icoords,top_gmax,dim);
    comp = top_comp[index];
    if (!ifluid_comp(comp)) return;

    double q0 = field_array[index];
	double q_edge[3][2];

	for (idir = 0; idir < dim; idir++)
	{
	    for (j = 0; j < dim; ++j)
        {
	    	icnb[j] = icoords[j];
            icnb_opp[j] = icoords[j];
            icnb_opp1[j] = icoords[j];
            icnb_opp2[j] = icoords[j];
        }

	    for (nb = 0; nb < 2; nb++)
	    {
            refl_side[nb] = NO;
	    	icnb[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
	    	icnb_opp[idir] = (nb == 0) ? icoords[idir] + 1 : icoords[idir] - 1;
	    	index_nb = d_index(icnb,top_gmax,dim);
	    	index_nb_opp = d_index(icnb_opp,top_gmax,dim);
	
            status = (*findStateAtCrossing)(front,icoords,dir[idir][nb],
                    comp,&intfc_state,&hs,crx_coords);

	    	if (status == NO_PDE_BOUNDARY)
            {
		        q_edge[idir][nb] = field_array[index_nb];
            }
	    	else if (status ==CONST_P_PDE_BOUNDARY)
            {
                //OUTLET
                q_edge[idir][nb] = getStateQ(intfc_state);
            }
	    	else if (status == CONST_V_PDE_BOUNDARY)
            {
                if(wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    //INLET
                    q_edge[idir][nb] = q0;
                }
                else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                        wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                {
                    /*
                    //Use one sided 3pt derivative
                    icnb_opp1[idir] = (nb == 0) ? icoords[idir] + 1 : icoords[idir] - 1;
                    index_oppnb1 = d_index(icnb_opp1,top_gmax,dim);
                    double q_oppnb1 = field_array[index_oppnb1];

                    icnb_opp2[idir] = (nb == 0) ? icoords[idir] + 2 : icoords[idir] - 2;
                    index_oppnb2 = d_index(icnb_opp2,top_gmax,dim);
                    double q_oppnb2 = field_array[index_oppnb2];

                    q_edge[idir][nb] = 3.0*q0 - 3.0*q_oppnb1 + q_oppnb2;
                    */

                    //grad(q) dot normal = 0   //TODO: is this correct bdry condition????? 
                    int icoords_ghost[MAXD];
                    for (int m = 0; m < dim; ++m)
                        icoords_ghost[m] = icoords[m];

                    icoords_ghost[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;

                    double coords_ghost[MAXD];
                    double coords_reflect[MAXD];

                    ////////////////////////////////////////////////////////////////////////
                    ///  matches Incompress_Solver_Smooth_Basis::setSlipBoundaryGNOR()  ///
                    //////////////////////////////////////////////////////////////////////

                    for (int m = 0; m < dim; ++m)
                    {
                        coords_ghost[m] = top_L[m] + icoords_ghost[m]*top_h[m];
                        coords_reflect[m] = coords_ghost[m];
                    }

                    double nor[MAXD];
                    FT_NormalAtGridCrossing(front,icoords,
                            dir[idir][nb],comp,nor,&hs,crx_coords);

                    //Reflect the ghost point through intfc-mirror at crossing.
                    //first reflect across the grid line containing intfc crossing,
                    coords_reflect[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];

                    //Reflect the displacement vector across the line
                    //containing the intfc normal vector
                    double v[MAXD];
                    double vn = 0.0;

                    for (int m = 0; m < dim; ++m)
                    {
                        v[m] = coords_reflect[m] - crx_coords[m];
                        vn += v[m]*nor[m];
                    }

                    for (int m = 0; m < dim; ++m)
                        v[m] = 2.0*vn*nor[m] - v[m];

                    //The desired reflected point
                    for (int m = 0; m < dim; ++m)
                        coords_reflect[m] = crx_coords[m] + v[m];
                    ///////////////////////////////////////////////////////////////////////

                    //Interpolate phi at the reflected point,
                    double q_reflect;
                    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field_array,
                            getStateQ,&q_reflect,&field_array[index]);
                        //FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field_array,
                        //        getStateQ,&q_reflect,nullptr);//default_ans uses nearest intfc state

                    q_edge[idir][nb] = q_reflect;
                }
                else if (is_bdry_hs(hs) && wave_type(hs) == NEUMANN_BOUNDARY)
                {
                    q_edge[idir][nb] = q0;
                }
                else
                {
                    printf("ERROR in computeFieldPointGradQ() : Unknown Boundary Type!\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

            }
	    }
	}

	for (i = 0; i < dim; ++i)
	    grad_field[i] = 0.5*(q_edge[i][1] - q_edge[i][0])/top_h[i];
}      /* end computeFieldPointGradQ */

//TODO: See NOTE below
void Incompress_Solver_Smooth_Basis::setReferencePressure()
{
	int i,j,k,index;
	double *pres = iFparams->field->pres;
        double diff_pres;

	switch(dim)
        {
	case 2:
	    index = d_index2d(imin,jmin,top_gmax);
        diff_pres = pres[index] - iFparams->ref_pres;
	    for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	    {
            index = d_index2d(i,j,top_gmax);
            if (!ifluid_comp(top_comp[index]))
                continue;
            pres[index] -= diff_pres;
            //NOTE:
            //pres[index] = pres[index] - diff_pres 
            //pres[index] = pres[index] - (pres[index] - iFparams-ref_pres)
            //pres[index] = iFparams->ref_pres
	    }
	    break;
	case 3:
        index = d_index3d(imin,jmin,kmin,top_gmax);
        diff_pres = pres[index] - iFparams->ref_pres;
	    for (k = 0; k <= top_gmax[2]; ++k)
	    for (j = 0; j <= top_gmax[1]; ++j)
        for (i = 0; i <= top_gmax[0]; ++i)
	    {
            index = d_index3d(i,j,k,top_gmax);
            if (!ifluid_comp(top_comp[index]))
                continue;
            pres[index] -= diff_pres;
            //NOTE:
            //pres[index] = pres[index] - diff_pres 
            //pres[index] = pres[index] - (pres[index] - iFparams-ref_pres)
            //pres[index] = iFparams->ref_pres
	    }
	    break;
	}
}	/* end setReferencePressure */

extern int ifluid_find_state_at_crossing(
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
	    else if (boundary_state_function(*hs) &&
		strcmp(boundary_state_function_name(*hs),
		"iF_splitBoundaryState") == 0)
	    	return CONST_V_PDE_BOUNDARY;
	    else //flow through bdry
	    	return CONST_P_PDE_BOUNDARY;
	}
}	/* ifluid_find_state_at_crossing */

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

extern boolean neumann_type_bdry(int w_type)
{
	switch (w_type)
	{
	case NEUMANN_BOUNDARY:
	case MOVABLE_BODY_BOUNDARY:
	case GROWING_BODY_BOUNDARY:
	case ICE_PARTICLE_BOUNDARY:
	    return YES;
	default:
	    return NO;
	}
}	/* end neumann_type_bdry */

void Incompress_Solver_Smooth_Basis::applicationSetComponent(void)
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

void Incompress_Solver_Smooth_Basis::applicationSetStates(void)
{
	double coords[MAXD];
	int *icoords;
	int j;
	int id;
	STATE state;
	int ave_comp;
	double p_intfc[MAXD],t[MAXD];
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	double dist;
	double **vel = field->vel;
	double *pres = field->pres;
	double *phi = field->phi;
	
	setDomain();

    int size = (int)cell_center.size();
	for (int i = 0; i < size; i++)
    {
        icoords = cell_center[i].icoords;
        if (cell_center[i].comp != -1 &&
            cell_center[i].comp != top_comp[i])
        {
            for (int j = 0; j < dim; ++j)
                coords[j] = top_L[j] + icoords[j]*top_h[j];
            id = d_index(icoords,top_gmax,dim);

            if (fabs(cell_center[i].comp - top_comp[i]) != 2) continue;

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
            
            if (!FT_FindNearestIntfcPointInRange(front,ave_comp,coords,
                    NO_BOUNDARIES,p_intfc,t,&hse,&hs,2)) continue;

            dist = 0.0;
            for (int j = 0; j < dim; ++j)
                dist += sqr(coords[j] - p_intfc[j]);
            dist = sqrt(dist);
            
            if (debugging("set_crossed_state"))
            {
                printf("coords  = %f %f %f\n",coords[0],coords[1],
                        coords[2]);
                printf("p_intfc = %f %f %f\n",p_intfc[0],p_intfc[1],
                        p_intfc[2]);
            }
            
            //TODO: should use hmin
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
                printf("Old pressure   = %f  Old phi   = %f\n",
                    pres[id],phi[id]);
                printf("Intfc pressure = %f  Intfc phi = %f\n",
                    state.pres,state.phi);
            }

            for (int j = 0; j < dim; ++j)
                vel[j][id] = state.vel[j];
            //double speed = sqrt(sqr(vel[0][id]) + sqr(vel[1][id]) + sqr(vel[2][id]));
        }
    }

	FT_FreeGridIntfc(front);
	FT_MakeGridIntfc(front);
}	/* end applicationSetStates */

static void initTestParams(Front *front)
{
	FILE *infile = fopen(InName(front),"r");
	IF_PARAMS *params = (IF_PARAMS*)front->extra1;

	CursorAfterString(infile,"Enter base directory name:");
            fscanf(infile,"%s",params->base_dir_name);
	(void) printf("%s\n",params->base_dir_name);
	CursorAfterString(infile,"Enter index of comaparing step:");
            fscanf(infile,"%d",&params->base_step);
	(void) printf("%d\n",params->base_step);
	fclose(infile);

}

void Incompress_Solver_Smooth_Basis::compareWithBaseSoln()
{
	int ii,jj,kk,index,k,l,DD;
        double P[MAXD],L1,L2,L_inf,err;
        double uex,*u;
        double u_max,u_min;
        IF_PARAMS *params = (IF_PARAMS*)front->extra1;
        IF_PARAMS *base_params;
        IF_FIELD *base_field;
        RECT_GRID *base_grid;
        int *base_gmax;
        int *base_comp;
        int base_ic[MAXD];
        Table *base_T;
        int N,base_index;

        L1 = L2 = L_inf = 0.0;

        u_max = -HUGE;
        u_min = HUGE;
        u = field->vel[0];

 	initTestParams(front);	
        readBaseFront(params,0);
        base_params = (IF_PARAMS*)base_front->extra1;
        base_field = base_params->field;
        base_grid = &topological_grid(base_front->grid_intfc);
        base_gmax = base_grid->gmax;
        base_T = table_of_interface(base_front->grid_intfc);
        base_comp = base_T->components;

        N = 0;
        switch (dim)
        {
        case 2:
            for (ii = imin; ii <= imax; ++ii)
            for (jj = jmin; jj <= jmax; ++jj)
            {
                double R;
                index = d_index2d(ii,jj,top_gmax);
                getRectangleCenter(index,P);
                rect_in_which(P,base_ic,base_grid);
                base_index = d_index(base_ic,base_gmax,2);
                if (top_comp[index] != base_comp[base_index])
                    continue;

                DD = (top_comp[index] == 2) ? 1 : -1;
                R = sqrt(sqr(P[0]) + sqr(P[1]));
                FT_IntrpStateVarAtCoords(base_front,top_comp[index],
                        P,base_field->vel[0],getStateXvel,&uex,NULL);
                if (u_max < u[index]) u_max = u[index];
                if (u_min > u[index]) u_min = u[index];
                err = fabs(u[index] - uex);
                L1 += err;
                L2 += sqr(err);
                if (err > L_inf)
                {
                    L_inf = err;
                }
                N++;
            }
            break;
        }
        L1 /= N;
        L2 /= N;
        L2 = sqrt(L2);
        printf("L1 = %18.16f  L2 = %18.16f  L_inf = %18.16f\n",
                                L1,L2,L_inf);
        printf("u_max = %f  u_min = %f\n",u_max,u_min);
}

void Incompress_Solver_Smooth_Basis::readBaseFront(
        IF_PARAMS *iF_params,
        int i)
{
        char *dir_name = iF_params->base_dir_name;
        F_BASIC_DATA *f_basic;
	int RestartStep = iF_params->base_step;
	int j;

        FT_ScalarMemoryAlloc((POINTER*)&base_front,sizeof(Front));
        FT_ScalarMemoryAlloc((POINTER*)&f_basic,sizeof(F_BASIC_DATA));
	
	f_basic->RestartRun = YES;
	f_basic->dim = dim;
        f_basic->size_of_intfc_state = sizeof(STATE);
	for (j = 0; j < dim; j++)
	    f_basic->subdomains[j] = front->pp_grid->gmax[j];

	FT_ReadComparisonDomain(InName(front),f_basic);

        sprintf(f_basic->restart_name,"%s/intfc-ts%s",dir_name,
                        right_flush(RestartStep,7));
        printf("restart_name = %s\n",f_basic->restart_name);

        FT_StartUp(base_front,f_basic);


        sprintf(f_basic->restart_state_name,"%s/state.ts%s",dir_name,
                        right_flush(RestartStep,7));
        readBaseStates(f_basic->restart_state_name);
	
}       /* end readBaseFront */

void Incompress_Solver_Smooth_Basis::readBaseStates(
        char *restart_name)
{
        FILE *infile;
        int i,j,k,index;
        char fname[100];
        double *u;
	double rho, pres, phi, mu, v; //dummy variable
        int *base_gmax;
        RECT_GRID *base_grid;
        static IF_PARAMS params;

	printf("Entering readBaseStates()\n");
        sprintf(fname,"%s-ifluid",restart_name);
        infile = fopen(fname,"r");

        /* Initialize states at interface and in the interior regions */
        //fluid_read_front_states(infile,base_front);

        FT_MakeGridIntfc(base_front);
        base_grid = &topological_grid(base_front->grid_intfc);
        base_gmax = base_grid->gmax;
        FT_ScalarMemoryAlloc((POINTER*)&params.field,sizeof(IF_FIELD));
        FT_VectorMemoryAlloc((POINTER*)&(params.field->vel),dim,sizeof(POINTER));

        next_output_line_containing_string(infile,"Interior ifluid states:");

        switch (dim)
        {
        case 1:
            FT_VectorMemoryAlloc((POINTER*)&u,base_gmax[0]+1,FLOAT);
            for (i = 0; i <= base_gmax[0]; ++i)
            {
                index = d_index1d(i,base_gmax);
		fscanf(infile,"%lf",&rho);
		fscanf(infile,"%lf",&pres);
		fscanf(infile,"%lf",&phi);
		fscanf(infile,"%lf",&mu);
                fscanf(infile,"%lf",&u[index]);
            }
            break;
        case 2:
            FT_VectorMemoryAlloc((POINTER*)&u,
                        (base_gmax[0]+1)*(base_gmax[1]+1),FLOAT);
            for (i = 0; i <= base_gmax[0]; ++i)
            for (j = 0; j <= base_gmax[1]; ++j)
            {
                index = d_index2d(i,j,base_gmax);
		fscanf(infile,"%lf",&rho);
		fscanf(infile,"%lf",&pres);
		fscanf(infile,"%lf",&phi);
		fscanf(infile,"%lf",&mu);
                fscanf(infile,"%lf",&u[index]);
                fscanf(infile,"%lf",&v);
            }
            break;
        case 3:
            FT_VectorMemoryAlloc((POINTER*)&u,(base_gmax[0]+1)*
                        (base_gmax[1]+1)*(base_gmax[2]+1),FLOAT);
            for (i = 0; i <= base_gmax[0]; ++i)
            for (j = 0; j <= base_gmax[1]; ++j)
            for (k = 0; k <= base_gmax[2]; ++k)
            {
                index = d_index3d(i,j,k,base_gmax);
		fscanf(infile,"%lf",&rho);
                fscanf(infile,"%lf",&pres);
                fscanf(infile,"%lf",&phi);
                fscanf(infile,"%lf",&mu);
                fscanf(infile,"%lf",&u[index]);
                fscanf(infile,"%lf",&v);
                fscanf(infile,"%lf",&v);
            }
            break;
        }
        params.field->vel[0] = u;
        base_front->extra1 = (POINTER)&params;
        fclose(infile);
}       /* end readBaseStates */

void Incompress_Solver_Smooth_Basis::recordVelocity()
{
    double** vel = field->vel;
    double** prev_vel = field->prev_vel;

	switch (dim)
	{
	case 2:
	    for (int j = 0; j <= top_gmax[1]; ++j)
	    for (int i = 0; i <= top_gmax[0]; ++i)
	    {
		    int index = d_index2d(i,j,top_gmax);
	        for (int l = 0; l < dim; ++l)
                prev_vel[l][index] = vel[l][index];
	    }
	    break;
	case 3:
	    for (int k = 0; k <= top_gmax[2]; ++k)
	    for (int j = 0; j <= top_gmax[1]; ++j)
	    for (int i = 0; i <= top_gmax[0]; ++i)
	    {
		    int index = d_index3d(i,j,k,top_gmax);
	        for (int l = 0; l < dim; ++l)
                prev_vel[l][index] = vel[l][index];
	    }
	}

    FT_ParallelExchGridVectorArrayBuffer(prev_vel,front);
}

void Incompress_Solver_Smooth_Basis::computeMaxSpeed(void)
{
	double speed;
	int i,j,k,l,index;
	double **vel = field->vel;

	max_speed = 0.0;
	for (i = 0; i < dim; ++i)
    {
        vmin[i] = HUGE;
        vmax[i] = -HUGE;
        abs_vmax[i] = -HUGE;
    }
	
    switch (dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	    {
            index = d_index2d(i,j,top_gmax);
            
            speed = 0.0;
            for (l = 0; l < dim; ++l)
            {
                speed += sqr(vel[l][index]);
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmin[l]));
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmax[l]));
            }
            speed = sqrt(speed);
		
            if (max_speed < speed) 
            {
                max_speed = speed;
                icrds_max[0] = i;
                icrds_max[1] = j;
            }
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	    {
            index = d_index3d(i,j,k,top_gmax);
            
            speed = 0.0;
            for (l = 0; l < dim; ++l)
            {
                speed += sqr(vel[l][index]);
                if (vmin[l] > vel[l][index]) vmin[l] = vel[l][index];
                if (vmax[l] < vel[l][index]) vmax[l] = vel[l][index];
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmin[l]));
                abs_vmax[l] = std::max(abs_vmax[l], std::abs(vmax[l]));
            }
            speed = sqrt(speed);

            if (max_speed < speed) 
            {
                max_speed = speed;
                icrds_max[0] = i;
                icrds_max[1] = j;
                icrds_max[2] = k;
            }
	    }
	    break;
	}

	pp_global_max(&max_speed,1);
	pp_global_max(abs_vmax,dim);
	pp_global_max(vmax,dim);
	pp_global_min(vmin,dim);
}	/* end computeMaxSpeed */

/*Compute pressure jump due to fabric/canopy porosity*/
/*Ronald Fedkiw, A Boundary Condition Capturing Methods 
  for Poisson Equation on Irregular Domains, JCP, 1999*/
double Incompress_Solver_Smooth_Basis::computeFieldPointPressureJump(
        int *icoords,
        double alpha,
        double beta)
{
        if (!iFparams->with_porosity) return 0.0;

        int nb, index_nb, idir, i;
        int top_gmin[MAXD], ipn[MAXD];
        double crx_coords[MAXD], coords[MAXD], grad_phi[MAXD];
        double vel_intfc[MAXD], vel_rel[MAXD], nor[MAXD],vec[MAXD];
        double Un = 0.0, ans = 0.0, side = 0.0, d_p = 0.0;
        
        double *phi = field->phi;
        double* mu = field->mu;
        double* div_U = field->div_U;
        double* q = field->q;
        
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        POINTER intfc_state;
        HYPER_SURF *hs;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *rgr = computational_grid(front->interf);
        POINTER state;
        bool is_intfc = false;
        
        top_gmin[0] = top_gmin[1] = top_gmin[2] = 0;
        for (i = 0; i < dim; i++)
            coords[i] = top_L[i] + icoords[i]*top_h[i];
        
        
        int index = d_index(icoords,top_gmax,dim);
        double mudiv = 0.5*mu[index]*div_U[index]; 
        double q0 = q[index];

        double rho = (field->rho[index] == 0) ?
            iFparams->rho2 : field->rho[index];

        int max_nb = (dim == 2) ? 4 : 6;
        for (nb = 0; nb < max_nb; nb++)
        {
            is_intfc = FT_NormalAtGridCrossing(front,icoords,dir[nb],
                                  top_comp[index],nor,&hs,crx_coords);
            
            if (is_intfc && dim == 2 && is_bdry(Curve_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 3 && is_bdry(Surface_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 2 &&
                     hsbdry_type(Curve_of_hs(hs)) == STRING_HSBDRY)
                continue;
            
            /*Ghost Fluid Method(GFM): p_+ - p_- = alpha*Un*nor + beta*(Un*nor)^2*/
            /*p+ and p- is defined by normal vector*/
            if (is_intfc && (wave_type(hs) == ELASTIC_BOUNDARY))
            {
                next_ip_in_dir(icoords,dir[nb],ipn,top_gmin,top_gmax);
                index_nb = d_index(ipn,top_gmax,dim);

                double mudiv_nb = 0.5*mu[index_nb]*div_U[index_nb]; 
                double q_nb = q[index_nb];

                /*get relative velocity*/
                FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir[nb],
                                top_comp[index],&state,&hs,crx_coords);
                
                double vel_fluid[MAXD] = {0.0};
                for (idir = 0 ; idir < dim; idir++)
                {
                    if (state)
                        vel_intfc[idir] = (*getStateVel[idir])(state);
                    else
                        vel_intfc[idir] = 0.0;

                    FT_IntrpStateVarAtCoords(front,NO_COMP,crx_coords,
                                             field->vel[idir],getStateVel[idir],
                                             &vel_fluid[idir],NULL);

                    vel_rel[idir] = vel_fluid[idir] - vel_intfc[idir];
                }

                /*project to normal direction*/
                Un = (dim == 2) ? Dot2d(nor,vel_rel) : Dot3d(nor,vel_rel);
        
                //NOTE: alpha and beta include thickness factor
                d_p = (alpha + fabs(Un)*beta)*Un;

                for (i = 0; i < dim; i++)
                    vec[i] = coords[i] - crx_coords[i];

                side = (dim == 2) ? Dot2d(nor,vec) : Dot3d(nor,vec);

                /////////////////////////////////////////////////////////////////
                //       d_p += jump_mudiv - jump_q;
                
                double jump_mudiv;
                double jump_q;

                if (side <= 0)
                {
                    jump_mudiv = mudiv_nb - mudiv;
                    jump_q = q_nb - q0;
                }
                else
                {
                    jump_mudiv = mudiv - mudiv_nb;
                    jump_q = q0 - q_nb;
                }

                d_p += jump_mudiv;

                if (iFparams->num_scheme.projc_method == PMI ||
                    iFparams->num_scheme.projc_method == PMII)
                {
                    d_p -= jump_q;
                }
                /////////////////////////////////////////////////////////////////

                if (side <= 0)
                {
                    ans -= d_p/sqr(top_h[nb/2])/rho;
                }
                else
                {
                    ans += d_p/(sqr(top_h[nb/2]))/rho;
                }
            
                if (debugging("pressure_drop"))
                {
                    printf("\ncomputeFieldPointPressureJump()\n");
                    printf("crds = [%f %f %f], crx = [%f %f %f],"
                            " side = %f, nb = %d\n",coords[0],coords[1],coords[2],
                            crx_coords[0],crx_coords[1],crx_coords[2],side,nb);
                    printf("vel_rel = [%f %f %f]",vel_rel[0],vel_rel[1],vel_rel[2]);
                    printf("d_p = %f, Un = %f, jump_mudiv = %f, jump_q = %f\n",
                            d_p, Un, jump_mudiv, jump_q);
                }
            }
        }
        
        /*return source term for Poisson equation due to jump condition*/
        return ans;
}

double Incompress_Solver_Smooth_Basis::computeFieldPointPressureJumpQ(
        int *icoords,
        double alpha,
        double beta)
{
        if (!iFparams->with_porosity) return 0.0;

        int nb, index_nb, idir, i;
        int top_gmin[MAXD], ipn[MAXD];
        double crx_coords[MAXD], coords[MAXD], grad_phi[MAXD];
        double vel_intfc[MAXD], vel_rel[MAXD], nor[MAXD],vec[MAXD];
        double Un = 0.0, ans = 0.0, side = 0.0, d_p = 0.0;
        
        double *phi = field->phi;
        double* mu = field->mu;
        double* div_U = field->div_U;
        double* q = field->q;
        
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        POINTER intfc_state;
        HYPER_SURF *hs;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *rgr = computational_grid(front->interf);
        POINTER state;
        bool is_intfc = false;
        
        top_gmin[0] = top_gmin[1] = top_gmin[2] = 0;
        for (i = 0; i < dim; i++)
            coords[i] = top_L[i] + icoords[i]*top_h[i];
        
        
        int index = d_index(icoords,top_gmax,dim);
        double mudiv = 0.5*mu[index]*div_U[index]; 
        double q0 = q[index];

        double rho = (field->rho[index] == 0) ?
            iFparams->rho2 : field->rho[index];

        int max_nb = (dim == 2) ? 4 : 6;
        for (nb = 0; nb < max_nb; nb++)
        {
            is_intfc = FT_NormalAtGridCrossing(front,icoords,dir[nb],
                                  top_comp[index],nor,&hs,crx_coords);
            
            if (is_intfc && dim == 2 && is_bdry(Curve_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 3 && is_bdry(Surface_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 2 &&
                     hsbdry_type(Curve_of_hs(hs)) == STRING_HSBDRY)
                continue;
            
            /*Ghost Fluid Method(GFM): p_+ - p_- = alpha*Un*nor + beta*(Un*nor)^2*/
            /*p+ and p- is defined by normal vector*/
            if(is_intfc && (wave_type(hs) == ELASTIC_BOUNDARY))
            {
                next_ip_in_dir(icoords,dir[nb],ipn,top_gmin,top_gmax);
                index_nb = d_index(ipn,top_gmax,dim);

                double mudiv_nb = 0.5*mu[index_nb]*div_U[index_nb]; 
                double q_nb = q[index_nb];

                /*get relative velocity*/
                FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir[nb],
                                top_comp[index],&state,&hs,crx_coords);
                
                double vel_fluid[MAXD] = {0.0};
                for (idir = 0 ; idir < dim; idir++)
                {
                    if (state)
                        vel_intfc[idir] = (*getStateVel[idir])(state);
                    else
                        vel_intfc[idir] = 0.0;

                    FT_IntrpStateVarAtCoords(front,NO_COMP,crx_coords,
                                             field->vel[idir],getStateVel[idir],
                                             &vel_fluid[idir],NULL);

                    vel_rel[idir] = vel_fluid[idir] - vel_intfc[idir];
                }

                /*project to normal direction*/
                Un = (dim == 2) ? Dot2d(nor,vel_rel) : Dot3d(nor,vel_rel);
        
                //NOTE: alpha and beta include thickness factor
                d_p = (alpha + fabs(Un)*beta)*Un;

                for (i = 0; i < dim; i++)
                    vec[i] = coords[i] - crx_coords[i];

                side = (dim == 2) ? Dot2d(nor,vec) : Dot3d(nor,vec);
                if (side <= 0)
                {
                    ans -= d_p/sqr(top_h[nb/2])/rho;
                }
                else
                {
                    ans += d_p/(sqr(top_h[nb/2]))/rho;
                }
            
                if (debugging("pressure_drop"))
                {
                    printf("\ncomputeFieldPointPressureJump()\n");
                    printf("crds = [%f %f %f], crx = [%f %f %f],"
                            " side = %f, nb = %d\n",coords[0],coords[1],coords[2],
                            crx_coords[0],crx_coords[1],crx_coords[2],side,nb);
                    printf("vel_rel = [%f %f %f]",vel_rel[0],vel_rel[1],vel_rel[2]);
                    printf("d_p = %f, Un = %f\n",d_p, Un);
                }
            }
        }
        
        /*return source term for Poisson equation due to jump condition*/
        return ans;
}

void Incompress_Solver_Smooth_Basis::computeFieldPointGradJump(
        int* icoords,
        double* var,
        double* grad_var)
{
        computeFieldPointGrad(icoords,var,grad_var);

        if (!iFparams->with_porosity) return;

        int index_nb, i;
        int top_gmin[MAXD], ipn[MAXD], nb;
        bool is_intfc = false;
        double side,coords[MAXD],crx_coords[MAXD],nor[MAXD],vec[MAXD];
        double vel_rel[MAXD], Un = 0.0;
        double vel_intfc[MAXD];
        double alpha = iFparams->porous_coeff[0];
        double beta = iFparams->porous_coeff[1];

        double* mu = field->mu;
        double* div_U = field->div_U;
        double* q = field->q;
        
        POINTER intfc_state;
        HYPER_SURF *hs;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *rgr = computational_grid(front->interf);
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        
        top_gmin[0] = top_gmin[1] = top_gmin[2] = 0;
        for (i = 0; i < dim; i++)
            coords[i] = top_L[i] + icoords[i]*top_h[i];

        
        int index = d_index(icoords,top_gmax,dim);
        double mudiv = 0.5*mu[index]*div_U[index]; 
        double q0 = q[index];


        double max_nb = (dim == 2) ? 4 : 6;
        for (nb = 0; nb < max_nb; nb++)
        {
            is_intfc = FT_NormalAtGridCrossing(front,icoords,dir[nb],
                                  top_comp[index],nor,&hs,crx_coords);
            
            if (is_intfc && dim == 2 && is_bdry(Curve_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 3 && is_bdry(Surface_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 2 &&
                     hsbdry_type(Curve_of_hs(hs)) == STRING_HSBDRY)
                continue;
            
            if (is_intfc && (wave_type(hs) == ELASTIC_BOUNDARY))
            {
                next_ip_in_dir(icoords,dir[nb],ipn,top_gmin,top_gmax);
                index_nb = d_index(ipn,top_gmax,dim);

                double mudiv_nb = 0.5*mu[index_nb]*div_U[index_nb]; 
                double q_nb = q[index_nb];

                FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir[nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
                
                double vel_fluid[MAXD] = {0.0};
                for (int idir = 0 ; idir < dim; idir++)
                {
                    vel_intfc[idir] = (*getStateVel[idir])(intfc_state);
                    FT_IntrpStateVarAtCoords(front,NO_COMP,crx_coords,
                                             field->vel[idir],getStateVel[idir],
                                             &vel_fluid[idir],NULL);
                    vel_rel[idir] = vel_fluid[idir] - vel_intfc[idir];
                }
                    
                //project to relative fluid velocity to intfc normal direction
                Un = (dim == 2) ? Dot2d(nor,vel_rel) : Dot3d(nor,vel_rel);
                
                //Compute interface pressure jump
                //  
                //  d_p = p^{+} - p^{-}
                //
                //given by the Ergun Equation                
                //
                //  d_p/d_r = alpha*Un + beta*|Un|*Un
                
                //NOTE: alpha and beta include thickness factor d_r
                double d_p = (alpha + fabs(Un)*beta)*Un;
		 
                for (i = 0; i < dim; i++)
                    vec[i] = coords[i] - crx_coords[i];
                
                side = (dim == 2) ? Dot2d(nor,vec) : Dot3d(nor,vec);

                /////////////////////////////////////////////////////////////////
                //                 d_p += jump_mudiv - jump_q;
                double jump_mudiv;
                double jump_q;

                if (side <= 0)
                {
                    jump_mudiv = mudiv_nb - mudiv;
                    jump_q = q_nb - q0;
                }
                else
                {
                    jump_mudiv = mudiv - mudiv_nb;
                    jump_q = q0 - q_nb;
                }
                
                d_p += jump_mudiv;

                if (iFparams->num_scheme.projc_method == PMI ||
                    iFparams->num_scheme.projc_method == PMII)
                {
                    d_p -= jump_q;
                }
                /////////////////////////////////////////////////////////////////
                
                // modify pressure gradient
                if ((side <= 0 && nb%2 == 0) || (side > 0 && nb%2 == 1))
                {
                    grad_var[nb/2] -= 0.5*d_p/top_h[nb/2];
                }
                else if ((side <= 0 && nb%2 == 1) || (side > 0 && nb%2 == 0))
                {
                    grad_var[nb/2] += 0.5*d_p/top_h[nb/2];
                }

                if (debugging("pressure_drop"))
                {
                    printf("\ncomputeFieldPointGradJump()\n");
                    printf("d_p = %f, vel_rel = [%f %f %f], Un = %f, "
                            "jump_mudiv = %f, jump_q = %f\n",
                            d_p,vel_rel[0],vel_rel[1],vel_rel[2],Un,jump_mudiv,jump_q);
                    printf("crds = [%f %f %f], crx = [%f %f %f], side = %f, nb = %d\n",
                            coords[0],coords[1],coords[2],
                            crx_coords[0],crx_coords[1],crx_coords[2],side,nb);
                }

            }
        }
}

void Incompress_Solver_Smooth_Basis::computeFieldPointGradJumpQ(
        int* icoords,
        double* var,
        double* grad_var)
{
        computeFieldPointGradQ(icoords,var,grad_var);

        if (!iFparams->with_porosity) return;

        int index_nb, i;
        int top_gmin[MAXD], ipn[MAXD], nb;
        bool is_intfc = false;
        double side,coords[MAXD],crx_coords[MAXD],nor[MAXD],vec[MAXD];
        double vel_rel[MAXD], Un = 0.0;
        double vel_intfc[MAXD];
        double alpha = iFparams->porous_coeff[0];
        double beta = iFparams->porous_coeff[1];

        double* mu = field->mu;
        double* div_U = field->div_U;
        double* q = field->q;
        
        POINTER intfc_state;
        HYPER_SURF *hs;
        INTERFACE *grid_intfc = front->grid_intfc;
        RECT_GRID *rgr = computational_grid(front->interf);
        GRID_DIRECTION dir[6] = {WEST,EAST,SOUTH,NORTH,LOWER,UPPER};
        
        top_gmin[0] = top_gmin[1] = top_gmin[2] = 0;
        for (i = 0; i < dim; i++)
            coords[i] = top_L[i] + icoords[i]*top_h[i];

        
        int index = d_index(icoords,top_gmax,dim);
        double mudiv = 0.5*mu[index]*div_U[index]; 
        double q0 = q[index];


        double max_nb = (dim == 2) ? 4 : 6;
        for (nb = 0; nb < max_nb; nb++)
        {
            is_intfc = FT_NormalAtGridCrossing(front,icoords,dir[nb],
                                  top_comp[index],nor,&hs,crx_coords);
            
            if (is_intfc && dim == 2 && is_bdry(Curve_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 3 && is_bdry(Surface_of_hs(hs)))
                continue;
            else if (is_intfc && dim == 2 &&
                     hsbdry_type(Curve_of_hs(hs)) == STRING_HSBDRY)
                continue;
            
            if (is_intfc && (wave_type(hs) == ELASTIC_BOUNDARY))
            {
                next_ip_in_dir(icoords,dir[nb],ipn,top_gmin,top_gmax);
                index_nb = d_index(ipn,top_gmax,dim);

                double mudiv_nb = 0.5*mu[index_nb]*div_U[index_nb]; 
                double q_nb = q[index_nb];

                FT_StateStructAtGridCrossing(front,grid_intfc,icoords,dir[nb],
                                top_comp[index],&intfc_state,&hs,crx_coords);
                
                double vel_fluid[MAXD] = {0.0};
                for (int idir = 0 ; idir < dim; idir++)
                {
                    vel_intfc[idir] = (*getStateVel[idir])(intfc_state);
                    FT_IntrpStateVarAtCoords(front,NO_COMP,crx_coords,
                                             field->vel[idir],getStateVel[idir],
                                             &vel_fluid[idir],NULL);
                    vel_rel[idir] = vel_fluid[idir] - vel_intfc[idir];
                }
                    
                //project to relative fluid velocity to intfc normal direction
                Un = (dim == 2) ? Dot2d(nor,vel_rel) : Dot3d(nor,vel_rel);
                
                //Compute interface pressure jump
                //  
                //  d_p = p^{+} - p^{-}
                //
                //given by the Ergun Equation                
                //
                //  d_p/d_r = alpha*Un + beta*|Un|*Un
                
                //NOTE: alpha and beta include thickness factor d_r
                double d_p = (alpha + fabs(Un)*beta)*Un;
		 
                for (i = 0; i < dim; i++)
                    vec[i] = coords[i] - crx_coords[i];
                
                side = (dim == 2) ? Dot2d(nor,vec) : Dot3d(nor,vec);

                // modify pressure gradient
                if ((side <= 0 && nb%2 == 0) || (side > 0 && nb%2 == 1))
                {
                    grad_var[nb/2] -= 0.5*d_p/top_h[nb/2];
                }
                else if ((side <= 0 && nb%2 == 1) || (side > 0 && nb%2 == 0))
                {
                    grad_var[nb/2] += 0.5*d_p/top_h[nb/2];
                }

                if (debugging("pressure_drop"))
                {
                    printf("\ncomputeFieldPointGradJumpQ()\n");
                    printf("d_p = %f, vel_rel = [%f %f %f], Un = %f\n",
                            d_p,vel_rel[0],vel_rel[1],vel_rel[2],Un);
                    printf("crds = [%f %f %f], crx = [%f %f %f], side = %f, nb = %d\n",
                            coords[0],coords[1],coords[2],
                            crx_coords[0],crx_coords[1],crx_coords[2],side,nb);
                }

            }
        }
}
int Incompress_Solver_Smooth_Basis::drawColorMap() {
	int num_colors = 0;
	int size;

	size = 1;
	for (int i = 0; i < dim; ++i)
	    size *= (top_gmax[i]+1);

	color_map.resize(size,0);
	color_map_copy.resize(size,0);
	std::fill(color_map.begin(), color_map.end(), 0);
        std::fill(color_map_copy.begin(), color_map_copy.end(), 0);

	//find connected regions locally
	paintInterior(num_colors);

	//modify colors to be unique globally
	makeGlobalColorMap(num_colors);

	return num_colors;
}	/* end drawColorMap */

void Incompress_Solver_Smooth_Basis::paintInterior(int &num_colors) 
{
	int i,ic,size;

   	paintAllGridPoint(NOT_SOLVED);
	//count number of colors in local domain
	//excluding solid regions labelled by 0
    	num_colors = 0; 
	size = top_gmax[0]+1;
	for (i = 1; i < dim; ++i)
	    size *= (top_gmax[i]+1);

        for (ic = 0; ic < size; ++ic) 
	{
	    //skip visited cells and solid cells
            if (domain_status[ic] != NOT_SOLVED ||
		top_comp[ic] == SOLID_COMP)
                continue;
	    paintConnectedRegion(ic,num_colors++);
	}
}	/* end paintInterior */

void Incompress_Solver_Smooth_Basis::makeGlobalColorMap(int &num_colors) 
{
	//adjust the local colors and make it unique globally
	std::vector<int> color_count(pp_numnodes(),0);
	std::vector<int> color_convert;
	int i,j;
	int start = 1; //reserved 0 for solid
	int no_change;

	if (debugging("paint_color"))
	    (void) printf(" Local num_colors = %d\n",num_colors);
	std::fill(color_count.begin(),color_count.end(),0);
	color_count[pp_mynode()] = num_colors;
	pp_global_imax(color_count.data(),pp_numnodes());
	num_colors = std::accumulate(color_count.begin(),color_count.end(),1);
	if (debugging("paint_color"))
	    (void) printf("Global num_colors = %d\n",num_colors);
	
	for (i = 0; i < pp_mynode(); ++i) 
	{
	    start += color_count[i];
	}
	for (i = 0; i < color_map.size(); ++i)
	{
	    if (top_comp[i] != SOLID_COMP)
	    {
            	color_map[i] += start;
	    	color_map_copy[i] = color_map[i];
	    }
	}

	color_convert.resize(num_colors,0);
	std::fill(color_convert.begin(),color_convert.end(),0);

	for (i = 0; i < num_colors; ++i)
	    color_convert[i] = i;

	FT_ParallelExchGridIntArrayBuffer(color_map.data(),front);
	for (i = 0; i < color_map.size(); ++i)
	{
	    if (color_map[i] != color_map_copy[i])
	    {
		if (color_map[i] >= color_map_copy[i]) continue;
		j = color_map_copy[i];
		if (color_convert[j] > color_map[i])
		    color_convert[j] = color_map[i];
	    }
	}
	no_change = NO;
	while (!no_change)
	{
	    pp_global_imin(color_convert.data(),num_colors);
	    no_change = YES;
	    for (i = 0; i < num_colors; ++i)
	    {
		j = color_convert[i];
		if (color_convert[i] > color_convert[j])
		{
		    no_change = NO;
		    color_convert[i] = color_convert[j];
		}
	    }
	    pp_global_imin(&no_change,1);
	}
	if (debugging("paint_color"))
	{
	    for (i = 0; i < num_colors; ++i)
	    	(void) printf("color_convert[%d] = %d\n",i,color_convert[i]);
	}
	for (i = 0; i < num_colors; ++i)
	{
	    if (i != color_convert[i])	// vacancy
	    {
		for (j = num_colors-1; j > i; --j)
		{
		    if (j == color_convert[j])
		    {
			color_convert[j] = i;
			break;
		    }
		}
	    }
	}
	if (debugging("paint_color"))
	{
	    (void) printf("After compressing the index:\n");
	    for (i = 0; i < num_colors; ++i)
	    	(void) printf("color_convert[%d] = %d\n",i,color_convert[i]);
	}
	for (i = 0; i < color_map.size(); ++i)
	{
	    j = color_map[i];
	    color_map[i] = color_convert[j];
	}
	num_colors = 0;
	for (i = 0; i < color_convert.size(); ++i)
	{
	    if (num_colors < color_convert[i])
		num_colors = color_convert[i];
	}
	num_colors++;
	if (debugging("paint_color"))
	    (void) printf("Final number of colors = %d\n",num_colors);
}	/* end makeGlobalColorMap */

void Incompress_Solver_Smooth_Basis::setDoubleDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i,j,k,ic,id,size,icoords[MAXD];
        COMPONENT copy_comp;

        grid_intfc = front->grid_intfc;
	if (first)
	{
            dim = grid_intfc->dim;
            D_extension = 2;

            for (i = 0; i < dim; ++i)
            {
                ext_gmax[i] = top_gmax[i];
                ext_imin[i] = (lbuf[i] == 0) ? 1 : lbuf[i];
                ext_l[i] = ext_u[i] = 0;
                if (grid_intfc->rect_bdry_type[i][0] == DIRICHLET_BOUNDARY)
                {
                    ext_l[i] = D_extension;
                    ext_gmax[i] += D_extension;
                }
                if (grid_intfc->rect_bdry_type[i][1] == DIRICHLET_BOUNDARY)
                {
                    ext_u[i] = D_extension;
                    ext_gmax[i] += D_extension;
                }
                ext_imax[i] = (ubuf[i] == 0) ? ext_gmax[i] - 1 : 
                            ext_gmax[i] - ubuf[i];
            }

	    size = ext_gmax[0]+1;
            for (i = 1; i < dim; ++i)
                size *= (ext_gmax[i]+1);

            FT_VectorMemoryAlloc((POINTER*)&ext_comp,size,sizeof(COMPONENT));
	    first = NO;
	}
        switch(dim)
        {
        case 2:
	    for (j = 0; j <= ext_gmax[1]; j++)
	    for (i = 0; i <= ext_gmax[0]; i++)
	    {
                id = d_index2d(i,j,ext_gmax);
                ext_comp[id] = FILL_COMP;
            }
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
                ic = d_index2d(i,j,top_gmax);
                id = d_index2d(i+ext_l[0],j+ext_l[1],ext_gmax);
                ext_comp[id] = top_comp[ic];
            }
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
            {
                icoords[0] = i;
                if (ext_l[1] != 0)
                {
                    icoords[1] = jmin + ext_l[1];
                    id = d_index(icoords,ext_gmax,2);
                    copy_comp = ext_comp[id];
                    for (j = 0; j <= ext_l[1]; ++j)
                    {
                        icoords[1] = j;
                        id = d_index(icoords,ext_gmax,2);
                        ext_comp[id] = copy_comp;
                    }
                }
                if (ext_u[1] != 0)
                {
                    icoords[1] = ext_gmax[1] - ext_u[1] - 1;
                    id = d_index(icoords,ext_gmax,2);
                    copy_comp = ext_comp[id];
                    for (j = 0; j <= ext_u[1]; ++j)
                    {
                        icoords[1] = ext_gmax[1] - j;
                        id = d_index(icoords,ext_gmax,2);
                        ext_comp[id] = copy_comp;
                    }
                }
            }
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
            {
                icoords[1] = j;
                if (ext_l[0] != 0)
                {
                    icoords[0] = imin + ext_l[0];
                    id = d_index(icoords,ext_gmax,2);
                    copy_comp = ext_comp[id];
                    for (i = 0; i <= ext_l[0]; ++i)
                    {
                        icoords[0] = i;
                        id = d_index(icoords,ext_gmax,2);
                        ext_comp[id] = copy_comp;
                    }
                }
                if (ext_u[0] != 0)
                {
                    icoords[0] = ext_gmax[0] - ext_u[0] - 1;
                    id = d_index(icoords,ext_gmax,2);
                    copy_comp = ext_comp[id];
                    for (i = 0; i <= ext_u[0]; ++i)
                    {
                        icoords[0] = ext_gmax[0] - i;
                        id = d_index(icoords,ext_gmax,2);
                        ext_comp[id] = copy_comp;
                    }
                }
            }
            break;
        case 3:
	    for (k = ext_imin[2]; k <= ext_imax[2]; k++)
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
                id = d_index3d(i,j,k,ext_gmax);
                ext_comp[id] = FILL_COMP;
            }
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
                ic = d_index3d(i,j,k,top_gmax);
                id = d_index3d(i+ext_l[0],j+ext_l[1],k+ext_l[2],ext_gmax);
                ext_comp[id] = top_comp[ic];
            }
            break;
        }
}	/* end setDoubleDomain */

void Incompress_Solver_Smooth_Basis::setDoubleGlobalIndex()
{
	int i,j,k,id;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&dn_dist,num_nodes,sizeof(int));
	}
	dNLblocks = 0;
	switch (dim)
	{
	case 2:
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
		id = d_index2d(i,j,ext_gmax);
		if (ext_comp[id] == SOLID_COMP) continue;
		dNLblocks++;
	    }
	    break;
	case 3:
	    for (k = ext_imin[2]; k <= ext_imax[2]; k++)
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
		id = d_index3d(i,j,k,ext_gmax);
		if (ext_comp[id] == SOLID_COMP) continue;
		dNLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	dn_dist[myid] = dNLblocks;
	pp_global_imax(dn_dist,num_nodes);
	eilower = 0;
        eiupper = dn_dist[0];

        for (i = 1; i <= myid; i++)
        {
            eilower += dn_dist[i-1];
            eiupper += dn_dist[i];
        }	
}	/* setDoubleGlobalIndex */

void Incompress_Solver_Smooth_Basis::setDoubleIndexMap(void)
{
	static boolean first = YES;
        static int com_gmax[MAXD];
        INTERFACE *grid_intfc = front->grid_intfc;
	int i,j,k,ic,index;
	int size = eiupper - eilower;
        int llbuf[MAXD],uubuf[MAXD];

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&dij_to_I,ext_gmax[0]+1,
				ext_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&dijk_to_I,ext_gmax[0]+1,
				ext_gmax[1]+1,ext_gmax[2]+1,INT);
	    	break;
	    }
            for (i = 0; i < dim; ++i)
            {
                com_gmax[i] = front->rect_grid->gmax[i];
                if (grid_intfc->rect_bdry_type[i][0] == DIRICHLET_BOUNDARY)
                    com_gmax[i] += D_extension;
                if (grid_intfc->rect_bdry_type[i][1] == DIRICHLET_BOUNDARY)
                    com_gmax[i] += D_extension;
            }
	}

	index = 0;
        for (i = 0; i < dim; ++i)
        {
            llbuf[i] = lbuf[i] != 0 ? lbuf[i] : 1;
            uubuf[i] = ubuf[i] != 0 ? ubuf[i] : 1;
        }
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= ext_gmax[1]; j++)
	    for (i = 0; i <= ext_gmax[0]; i++)
		    dij_to_I[i][j] = -1;
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
		ic = d_index2d(i,j,ext_gmax);
                if (ext_comp[ic] != SOLID_COMP)
                {
                    dij_to_I[i][j] = index + eilower;
                    index++;
                }
	    }
	    FT_ParallelExchExtendedCellIndex(front,llbuf,uubuf,
				com_gmax,(POINTER)dij_to_I);
	    break;
	case 3:
	    for (k = 0; k <= ext_gmax[2]; k++)
	    for (j = 0; j <= ext_gmax[1]; j++)
	    for (i = 0; i <= ext_gmax[0]; i++)
		    dijk_to_I[i][j][k] = -1;
	    for (k = ext_imin[2]; k <= ext_imax[2]; k++)
	    for (j = ext_imin[1]; j <= ext_imax[1]; j++)
	    for (i = ext_imin[0]; i <= ext_imax[0]; i++)
	    {
		ic = d_index3d(i,j,k,ext_gmax);
		/*
		if (domain_status[ic] != TO_SOLVE)
		    continue;
		*/
                if (ext_comp[ic] != SOLID_COMP)
		{
                    dijk_to_I[i][j][k] = index + eilower;
                    index++;
                }
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,
				(POINTER)dijk_to_I);
	    break;
	}
}	/* end setDoubleIndexMap */

void Incompress_Solver_Smooth_Basis::setSlipBoundary(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_slip)
{
    setSlipBoundaryNIP(icoords,idir,nb,comp,hs,state,vel,v_slip);
    
    //TODO: Write GNOR implementation and compare results.
    //
    //      setSlipBoundaryGNOR(icoords,idir,nb,comp,hs,state,vel,v_slip);
}

// Based on finding the nearest interface point to the ghost point
// computed using FT_FindNearestIntfcPointInRange()
void Incompress_Solver_Smooth_Basis::setSlipBoundaryNIP(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_slip)
{
    GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    int ghost_ic[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coords_reflect[MAXD], coords_ghost[MAXD];
    double nor[MAXD];
    
    double vel_intfc_gcrx[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        vel_intfc_gcrx[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ghost_ic[i] = icoords[i];
    }
    
    ghost_ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    int ghost_index = d_index(ghost_ic,top_gmax,dim);
    COMPONENT ghost_comp = top_comp[ghost_index];
    
    for (int j = 0; j < dim; ++j)
        coords_ghost[j] = top_L[j] + ghost_ic[j]*top_h[j];

    ////////////////////////////////////////////////////////////////////////
    double intrp_coeffs[MAXD] = {0.0};
    HYPER_SURF_ELEMENT* hsurf_elem;
    HYPER_SURF* hsurf;
    double range = 2;

    //TODO: Why does this fail for INCLUDE_BOUNDARIES and NO_SUBDOMAIN values?
    //      Conversely, why does it work with NO_BOUNDARIES in the backward facing
    //      step scenario -- to what degree is it working?
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    /*      
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,INCLUDE_BOUNDARIES,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    FT_FindNearestIntfcPointInRange(front,ghost_comp,coords_ghost,NO_SUBDOMAIN,
            crx_coords,intrp_coeffs,&hsurf_elem,&hsurf,range);
    */

    double dist_ghost = distance_between_positions(coords_ghost,crx_coords,dim);
    
    //compute the normal and velocity vectors at the interface point
    double vel_intfc[MAXD] = {0.0};
    switch (dim)
	{
        case 2:
            {
                double ns[MAXD] = {0.0};
                double ne[MAXD] = {0.0};
                
                normal(Bond_of_hse(hsurf_elem)->start,hsurf_elem,hsurf,ns,front);
                normal(Bond_of_hse(hsurf_elem)->end,hsurf_elem,hsurf,ne,front);

                /*
                STATE* ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                STATE* se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);
                */

                /////////////////////////////////////////////////////////////
                STATE* ss;
                STATE* se;

                if (ifluid_comp(negative_component(hsurf)))
                {
                    ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);
                }
                else if (ifluid_comp(positive_component(hsurf)))
                {
                    ss = (STATE*)right_state(Bond_of_hse(hsurf_elem)->start);
                    se = (STATE*)right_state(Bond_of_hse(hsurf_elem)->end);
                }
                else
                {
                    printf("setSlipBoundaryNIP() ERROR: "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }
                /////////////////////////////////////////////////////////////

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                    vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
                }
            }
            break;

        case 3:
            {
                TRI* nearTri = Tri_of_hse(hsurf_elem);
                const double* tnor = Tri_normal(nearTri);
                //NOTE: Tri_normal() does not return a unit vector
                
                STATE* st[3];

                if (ifluid_comp(negative_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);
                }
                else if (ifluid_comp(positive_component(hsurf)))
                {
                    for (int j = 0; j < 3; ++j)
                        st[j] = (STATE*)right_state(Point_of_tri(nearTri)[j]);
                }
                else
                {
                    printf("setSlipBoundaryNIP() ERROR: "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

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

    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hsurf))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    
    //NOTE: must use unit-length vectors with FT_GridSizeInDir()
    double dist_reflect = FT_GridSizeInDir(nor,front);

    //TODO: need to check if dist_reflect > dist_ghost ???
    
        /*
        // Compute dist_reflect as the diagonal length of rect grid blocks
        double dist_reflect = 0.0;
        for (int j = 0; j < 3; ++j)
             dist_reflect += sqr(top_h[j]);
        dist_reflect = sqrt(dist_reflect);
        */

    //The desired reflected point
    for (int j = 0; j < dim; ++j)
        coords_reflect[j] = crx_coords[j] + dist_reflect*nor[j];
    ////////////////////////////////////////////////////////////////////////
   

    ////////////////////////////////////////////////////////////////////////
    if (debugging("slip_boundary"))
    {
        printf("\nsetSlipBoundaryNIP() DEBUGGING\n");
        printf("idir = %d nb = %d\n",idir,nb);
        fprint_int_vector(stdout,"icoords",icoords,dim,", ");
        fprint_int_vector(stdout,"ghost_ic",ghost_ic,dim,"\n");
        fprint_general_vector(stdout,"coords",coords,dim,"\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",crx_coords,dim,"\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g\n",dist_ghost,dist_reflect);
        printf("dist_ghost/dist_reflect = %g  dist_reflect - dist_ghost = %g\n",
                dist_ghost/dist_reflect, dist_reflect - dist_ghost);
    }
    ////////////////////////////////////////////////////////////////////////

    // Interpolate the velocity at the reflected point
    int index = d_index(icoords,top_gmax,dim);
    double vel_reflect[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,vel[j],
                getStateVel[j],&vel_reflect[j],&vel[j][index]);
    }
 
    double vel_rel[MAXD] = {0.0};
    double vn = 0.0;

    for (int j = 0; j < dim; ++j)
    {
        vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        vn += vel_rel[j]*nor[j];
    }

    double vel_rel_tan[MAXD] = {0.0};
    double vel_rel_nor[MAXD] = {0.0};
    double vel_ghost_nor[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_rel_nor[j] = vn*nor[j];
	    vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }
    double mag_vtan = Magd(vel_rel_tan,dim);

    /////////////////////////////////////////////////////////////////////////
    if (iFparams->use_eddy_visc == NO)
    {
        for (int j = 0; j < dim; ++j)
            v_slip[j] = vel_reflect[j] + vel_ghost_nor[j];
        return;
    }
    /////////////////////////////////////////////////////////////////////////
    
    if (debugging("slip_boundary"))
    {
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_nor",vel_rel_nor,dim,"\n");
        printf("Magd(vel_rel_tan,dim) = %g\n",mag_vtan);
    }

    double mu_l;
    double rho_l;

    switch (comp)
    {
        case LIQUID_COMP1:
            mu_l = m_mu[0];
            rho_l = m_rho[0];
            break;
        case LIQUID_COMP2:
            mu_l = m_mu[1];
            rho_l = m_rho[1];
            break;
        default:
            printf("Unknown fluid COMPONENT: %d\n",comp);
            LOC(); clean_up(EXIT_FAILURE);
            break;
    }
    
    double tau_wall[MAXD] = {0.0};
    double mag_tau_wall = computeWallShearStress(mag_vtan,
                    dist_reflect,mu_l,rho_l,45.0);
    //NOTE: In all numerical experiments, Newton's method converged
    //      when the initial guess for the dimensionless wall velocity
    //      was in the range of 40-50.

    if (mag_vtan > MACH_EPS)
    {
        for (int j = 0; j < dim; ++j)
            tau_wall[j] = mag_tau_wall*vel_rel_tan[j]/mag_vtan;
    }

    // Interpolate the effective viscosity at the reflected point
    double mu_reflect;
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field->mu,
                getStateMu,&mu_reflect,&field->mu[index]);
    if (mu_reflect < MACH_EPS) mu_reflect = field->mu[index];
    
    double vel_ghost_tan[MAXD] = {0.0};
    double vel_ghost_rel[MAXD] = {0.0};
    
    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] = vel_rel_tan[j]
            - (dist_reflect - dist_ghost)/mu_reflect*tau_wall[j];

        vel_ghost_rel[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
        v_slip[j] = vel_ghost_rel[j] + vel_intfc[j];
    }


    if (debugging("slip_boundary"))
    {
        printf("mu_reflect = %g , mu_[%d] = %g\n",mu_reflect,index,field->mu[index]);
        printf("mag_tau_wall = %g\n",mag_tau_wall);
        fprint_general_vector(stdout,"tau_wall",tau_wall,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_tan",vel_ghost_tan,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_nor",vel_ghost_nor,dim,"\n");
        fprint_general_vector(stdout,"vel_ghost_rel",vel_ghost_rel,dim,"\n");
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
    }
    
    
    //store data to avoid recomputing values in the fluid solver
    int fid = face_index(idir,nb);
    for (int i = 0; i < dim; ++i)
    {
        ghost_data[fid][index].vel[i] = v_slip[i];
        ghost_data[fid][index].force[i] = tau_wall[i];
    }
}

// Based on the normal at the interface grid crossing
// computed by FT_NormalAtGridCrossing().
//
//TODO: Add the missing features that have already been
//      implemented in setSlipBoundaryNIP().
//      Compare the implementations when complete. 
void Incompress_Solver_Smooth_Basis::setSlipBoundaryGNOR(
	int *icoords,
	int idir,
	int nb,
	int comp,
	HYPER_SURF *hs,
	POINTER state,
	double** vel,
	double* v_slip)
{
    GRID_DIRECTION  ldir[3] = {WEST,SOUTH,LOWER};
    GRID_DIRECTION  rdir[3] = {EAST,NORTH,UPPER};

    int ghost_ic[MAXD];
    double coords[MAXD], crx_coords[MAXD];
    double coords_reflect[MAXD], coords_ghost[MAXD];
    double nor[MAXD];
    
    double vel_intfc_gcrx[MAXD];
    for (int i = 0; i < dim; ++i)
    {
        vel_intfc_gcrx[i] = (*getStateVel[i])(state);
        coords[i] = top_L[i] + icoords[i]*top_h[i];
        ghost_ic[i] = icoords[i];
    }
    
    ghost_ic[idir] = (nb == 0) ? icoords[idir] - 1 : icoords[idir] + 1;
    
    for (int j = 0; j < dim; ++j)
    {
        coords_ghost[j] = top_L[j] + ghost_ic[j]*top_h[j];
        coords_reflect[j] = coords_ghost[j];
    }

    // Reflect ghost point through intfc-mirror at crossing
    GRID_DIRECTION dir = (nb == 0) ? ldir[idir] : rdir[idir];
    FT_NormalAtGridCrossing(front,icoords,dir,comp,nor,&hs,crx_coords);
	
    // first reflect across the grid line containing intfc crossing
    coords_reflect[idir] = 2.0*crx_coords[idir] - coords_ghost[idir];
    
    double vec_reflect[MAXD] = {0.0};
    double vec_midpoint[MAXD] = {0.0};
    double vn = 0.0;
    
    for (int j = 0; j < dim; ++j)
    {
        vec_reflect[j] = coords_reflect[j] - crx_coords[j];
        vec_midpoint[j] = coords_reflect[j] - crx_coords[j];
        vn += vec_reflect[j]*nor[j];
    }

    for (int j = 0; j < dim; ++j)
        vec_reflect[j] = 2.0*vn*nor[j] - vec_reflect[j];

    double coords_midpoint[MAXD] = {0.0};

    //The desired reflected point, and the midpoint between it and the ghost point
    for (int j = 0; j < dim; ++j)
    {
        coords_reflect[j] = crx_coords[j] + vec_reflect[j];
        coords_midpoint[j] = coords_reflect[j] - vn*nor[j];
        //coords_midpoint should be identical to 1/2(coords_reflect + coords_ghost)
    }
    
    // compute dist_ghost and dist_reflect
    double coords_nip[MAXD] = {0.0};
    double intrp_coeffs[MAXD] = {0.0};
    HYPER_SURF_ELEMENT* hsurf_elem;
    HYPER_SURF* hsurf;
    double range = 2;

    FT_FindNearestIntfcPointInRange(front,comp,coords_midpoint,NO_BOUNDARIES,
            coords_nip,intrp_coeffs,&hsurf_elem,&hsurf,range);

    double dist_ghost = distance_between_positions(coords_ghost,coords_nip,dim);
    double dist_reflect = distance_between_positions(coords_reflect,coords_nip,dim);
    
    //TODO: Should we set dist_reflect to length of grid block diagonal?
    //      how does FT_GridSizeInDir() or the previous method differ? 
    
    /*
    // Compute dist_reflect as the diagonal length of rect grid blocks
    double dist_reflect = 0.0;
    for (int j = 0; j < 3; ++j)
         dist_reflect += sqr(top_h[j]);
    dist_reflect = sqrt(dist_reflect);
    */

    ////////////////////////////////////////////////////////////////////////
    //Temp debugging
    if (debugging("slip_boundary"))
    {
        printf("\nsetSlipBoundaryGNOR() DEBUGGING\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",coords_nip,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g , dist_ghost/dist_reflect = %g\n",
                dist_ghost, dist_reflect, dist_ghost/dist_reflect);
    }
    ////////////////////////////////////////////////////////////////////////
   
    
    ////////////////////////////////////////////////////////////////////////
    // Recompute the normal vector and interface velocity at coords_nip
    // before computing th the normal component of relative velocity with
    // respect to the intfc.
    
    double vel_intfc[MAXD] = {0.0};
    switch (dim)
	{
        case 2:
            {
                double ns[MAXD] = {0.0};
                double ne[MAXD] = {0.0};
                
                normal(Bond_of_hse(hsurf_elem)->start,hsurf_elem,hsurf,ns,front);
                normal(Bond_of_hse(hsurf_elem)->end,hsurf_elem,hsurf,ne,front);

                STATE* ss = (STATE*)left_state(Bond_of_hse(hsurf_elem)->start);
                STATE* se = (STATE*)left_state(Bond_of_hse(hsurf_elem)->end);

                for (int i = 0; i < dim; ++i)
                {
                    nor[i] = (1.0 - intrp_coeffs[0])*ns[i] + intrp_coeffs[0]*ne[i];
                    vel_intfc[i] = (1.0 - intrp_coeffs[0])*ss->vel[i] + intrp_coeffs[0]*se->vel[i];
                }
            }
            break;

        case 3:
            {
                TRI* nearTri = Tri_of_hse(hsurf_elem);
                const double* tnor = Tri_normal(nearTri);
                //NOTE: Tri_normal() does not normalize the normal vector.
                
                STATE* st[3];
                for (int j = 0; j < 3; ++j)
                    st[j] = (STATE*)left_state(Point_of_tri(nearTri)[j]);

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

    double mag_nor = Magd(nor,dim);
    for (int i = 0; i < dim; ++i)
        nor[i] /= mag_nor;
        
    if (comp == negative_component(hsurf))
	{
	    for (int i = 0; i < dim; ++i)
            nor[i] *= -1.0;
	}
    ////////////////////////////////////////////////////////////////////////

    double vel_reflect[MAXD] = {0.0};
    double mu_reflect;
    
    int index = d_index(icoords,top_gmax,dim);

    /* Interpolate the state at the reflected point */
    for (int j = 0; j < dim; ++j)
    {
        //TODO: Is vel[j][index] the best default value upon failure??
        FT_IntrpStateVarAtCoords(front,comp,coords_reflect,vel[j],
                getStateVel[j],&vel_reflect[j],&vel[j][index]);
    }
    FT_IntrpStateVarAtCoords(front,comp,coords_reflect,field->mu,
                getStateMu,&mu_reflect,&field->mu[index]);
 
    double vel_rel[MAXD] = {0.0};
    vn = 0.0;

    for (int j = 0; j < dim; ++j)
    {
        //vel_rel[j] = vel_reflect[j] - vel_intfc_gcrx[j];
        vel_rel[j] = vel_reflect[j] - vel_intfc[j];
        vn += vel_rel[j]*nor[j];
    }

    //////////////////////////////////////////////////////////////////////////////
    //TODO: This would be omitted if lower code implementation finished
    for (int j = 0; j < dim; ++j)
    {
        v_slip[j] = vel_reflect[j] - (dist_ghost/dist_reflect)*vn*nor[j];
            //v_slip[j] = vel_reflect[j] - vn*nor[j]; //NOTE: is just the tangential velocity
    }
    //////////////////////////////////////////////////////////////////////////////
    

    //TODO: CONTINUE WRITING IMPLEMENTATION BELOW
    //
    //NOTE: setSlipBoundaryNIP() development is ahead
    //      and should be copied/merged here.
    
    /*
    double vel_rel_tan[MAXD] = {0.0};
    double vel_ghost_nor[MAXD] = {0.0};

    for (int j = 0; j < dim; ++j)
    {
	    vel_rel_tan[j] = vel_rel[j] - vn*nor[j];
	    vel_ghost_nor[j] = -1.0*(dist_ghost/dist_reflect)*vn*nor[j];
    }
    double mag_vtan = Magd(vel_rel_tan,dim);

    //TODO: Need to modify tangential velocity also.
    //      First need to solve spalding's wall law for the
    //      friction velocity, u_tau, then compute the wall
    //      shear stress, tau_wall. Then must modify the ghost's
    //      relative tangential slip velocity according to
    //
    //      v_tan_ghost = v_tan_reflect
    //            - (dist_reflect + dist_ghost)/(mu_reflect + mu_t_reflect)*tau_wall
    //
    //      where v_tan_ghost and v_tan_reflect are relative velocities
    //      with respect to the interface.
    
    if (debugging("slip_boundary"))
    {
        printf("setSlipBoundaryGNOR() DEBUGGING\n");
        fprint_general_vector(stdout,"coords",coords,dim,"\n");
        fprint_general_vector(stdout,"coords_ghost",coords_ghost,dim,"\n");
        fprint_general_vector(stdout,"coords_nip",coords_nip,dim,"\n");
        fprint_general_vector(stdout,"normal",nor,dim,"\n");
        fprint_general_vector(stdout,"coords_reflect",coords_reflect,dim,"\n");
        printf("dist_ghost = %g , dist_reflect = %g , dist_ghost/dist_reflect = %g\n",
                dist_ghost, dist_reflect, dist_ghost/dist_reflect);
        fprint_general_vector(stdout,"vel_reflect",vel_reflect,dim,"\n");
        fprint_general_vector(stdout,"vel_intfc",vel_intfc,dim,"\n");
        fprint_general_vector(stdout,"vel_rel_tan",vel_rel_tan,dim,"\n");
    }

    double mul;
    double rhol;

    switch (comp)
    {
        case LIQUID_COMP1:
            mul = m_mu[0];
            rhol = m_rho[0];
            break;
        case LIQUID_COMP2:
            mul = m_mu[1];
            rhol = m_rho[1];
            break;
        default:
            printf("Unknown fluid COMPONENT: %d\n",comp);
            LOC(); clean_up(EXIT_FAILURE);
    }
    
    double tau_wall = computeWallShearStress(mag_vtan,dist_reflect,mul,rhol);

    double vel_ghost_tan[MAXD] = {0.0};
    for (int j = 0; j < dim; ++j)
    {
        vel_ghost_tan[j] =
            vel_rel_tan[j] - (dist_reflect - dist_ghost)/mu_reflect*tau_wall;

        //TODO: need to revert back to world frame by adding back vel_intfc
        v_slip[j] = vel_ghost_tan[j] + vel_ghost_nor[j] + vel_intfc[j];
            //v_slip[j] = vel_ghost_tan[j] + vel_ghost_nor[j];
            //v_slip[j] += vel_intfc[j];
    }

    if (debugging("slip_boundary"))
    {
        printf("tau_wall = %f\n",tau_wall);
        fprint_general_vector(stdout,"v_slip",v_slip,dim,"\n");
        printf("\n");
    }
    */
}

std::vector<double> Incompress_Solver_Smooth_Basis::computeGradPhiTangential(
        int* icoords,
        GRID_DIRECTION dir,
        COMPONENT comp,
        HYPER_SURF *hs,
        double* crx_coords)
{
    double** grad_phi = field->grad_phi;
    int index = d_index(icoords,top_gmax,dim);
    
    double nor[MAXD] = {0.0};
    FT_NormalAtGridCrossing(front,icoords,dir,
            comp,nor,&hs,crx_coords);

    double vn = 0.0;
    for (int j = 0; j < dim; ++j)
        vn += grad_phi[j][index]*nor[j];

    std::vector<double> grad_phi_tangent(3,0.0);
    for (int j = 0; j < dim; ++j) 
        grad_phi_tangent[j] = grad_phi[j][index] - vn*nor[j];

    return grad_phi_tangent;
}

//TODO: Currently only looking for CONSTANT state Dirichlet
//      boundaries corresponding to an INLET boundary.
void Incompress_Solver_Smooth_Basis::setFreeStreamVelocity()
{
    U_FreeStream = 0.0;

    HYPER_SURF* hs;
    for (int idir = 0; idir < dim; ++idir)
    for (int side = 0; side < 2; ++side)
    {
        hs = FT_RectBoundaryHypSurf(front->interf,
                DIRICHLET_BOUNDARY,idir,side);
        
        if (hs == nullptr) continue;
        if (boundary_state(hs) == nullptr) continue;

        STATE* bstate = (STATE*)boundary_state(hs);
        double* bdry_vel = bstate->vel;

        double bdry_speed = 0.0;
        for (int l = 0; l < dim; ++l)
        {
            bdry_speed += sqr(bdry_vel[l]);
        }
        bdry_speed = sqrt(bdry_speed);

        if (U_FreeStream < bdry_speed) 
        {
            U_FreeStream = bdry_speed;
        }
    }

	pp_global_max(&U_FreeStream,1);
}

