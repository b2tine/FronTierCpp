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


#include "fabric.h"
#include "bending.h"

static void naturalStressOfTri(TRI*,double);
static void singleCanopyModification(Front*);
static void bifurcateCanopyModification(Front*);
static void copyParachuteSet(ELASTIC_SET,ELASTIC_SET*);
static void rotateParachuteSet(ELASTIC_SET*,double*,double,double);

void printAfExtraData(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int i,dim = intfc->dim;
	FILE *outfile;
	char filename[200];
	SURFACE **s;
	CURVE **c;
	NODE **n;
	BOND *b;
    TRI *t;

	sprintf(filename,"%s/state.ts%s",out_name,
                        right_flush(front->step,7));
#if defined(HAVE_MPI)
    if (pp_numnodes() > 1)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(HAVE_MPI) */
    sprintf(filename,"%s-afdata",filename);
    outfile = fopen(filename,"w");

	fprintf(outfile,"\nAirfoil extra front state data:\n");

    //TODO: When done debugging can package up the functionality in here.
    //      FT_WriteFrontState(outfile,front);
    
	next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        /*
        if (wave_type(hs) != ELASTIC_BOUNDARY &&
            wave_type(hs) != ELASTIC_STRING) continue;
        */

        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->force[i]);
	    fprintf(outfile,"\n");

        FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
        fwrite(sl,1,front->sizest,outfile);
        fwrite(sr,1,front->sizest,outfile);
        fprintf(outfile,"\n");
        
        /*
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g %24.18g\n",sl->impulse[i],sr->impulse[i]);
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g %24.18g\n",sl->vel[i],sr->vel[i]);
        */
    }

	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;
        p = b->start;
        
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->force[i]);
	    fprintf(outfile,"\n");
        
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
        fwrite(sl,1,front->sizest,outfile);
        fwrite(sr,1,front->sizest,outfile);
	    fprintf(outfile,"\n");

        /*
        fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sl->vel[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sr->vel[i]);
	    fprintf(outfile,"\n");
        */
	    
        for (b = (*c)->first; b != NULL; b = b->next)
	    {
		    p = b->end;
	    	
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
	    	fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->force[i]);
	    	fprintf(outfile,"\n");
            
            sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            fwrite(sl,1,front->sizest,outfile);
            fwrite(sr,1,front->sizest,outfile);
            fprintf(outfile,"\n");
            
            /*
            fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    	fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    	fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->vel[i]);
	    	fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->vel[i]);
	    	fprintf(outfile,"\n");
            */
	    }
	}
	
    for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
        
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",p->force[i]);
	    fprintf(outfile,"\n");
	    
        sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
        fwrite(sl,1,front->sizest,outfile);
        fwrite(sr,1,front->sizest,outfile);
	    fprintf(outfile,"\n");
        
        /*
        fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sl->vel[i]);
	    fprintf(outfile,"\n");
        for (i = 0; i < dim; ++i)
            fprintf(outfile,"%24.18g ",sr->vel[i]);
	    fprintf(outfile,"\n");
        */
	}

	fprintf(outfile,"\nSurface extra data:\n");
    intfc_surface_loop(intfc,s) 
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_STRING) continue;

        int num_pts;
        REGISTERED_PTS *registered_pts;

        if ((*s)->extra == nullptr)
            num_pts = 0;
        else
        {
            registered_pts = (REGISTERED_PTS*)(*s)->extra;
            num_pts = registered_pts->num_pts;
        }
        fprintf(outfile,"number of registered points = %d\n",num_pts);
        for (i = 0; i < num_pts; ++i)
            fprintf(outfile,"%d\n",registered_pts->global_ids[i]);
    }
	
    fprintf(outfile,"\nCurve extra data:\n");
    intfc_curve_loop(intfc,c)
	{
		if (hsbdry_type(*c) != STRING_HSBDRY) continue;
        
        //TODO: C_PARAMS no longer being used, and should be
        //      replaced with the appropiate data structure for
        //      point mass runs (load_nodes)
	    C_PARAMS *c_params = (C_PARAMS*)(*c)->extra;
	    if (c_params == NULL)
        {
            fprintf(outfile,"curve extra: no\n");
        }
        else
	    {
            fprintf(outfile,"curve extra: yes\n");
            fprintf(outfile,"point_mass = %24.18g\n",c_params->point_mass);
            fprintf(outfile,"load_mass = %24.18g\n",c_params->load_mass);
            fprintf(outfile,"load_type = %d\n",c_params->load_type);
            fprintf(outfile,"dir = %d\n",c_params->dir);
	    }

        for (b = (*c)->first; b != (*c)->last; b = b->next)
        {
		    p = b->end;
            if (p->extra == nullptr)
            {
                fprintf(outfile,"string point extra: no\n");
            }
            else
            {
                fprintf(outfile,"string point extra: yes\n");
                BOND_BENDER* bond_bender = (BOND_BENDER*)p->extra;
                fwrite(bond_bender,1,sizeof(BOND_BENDER),outfile);
                fprintf(outfile,"\n");
            }
        }
	}
	
    fprintf(outfile,"\nNode extra data:\n");
    intfc_node_loop(intfc,n)
	{
	    AF_NODE_EXTRA *n_params = (AF_NODE_EXTRA*)(*n)->extra;
	    if (n_params == NULL)
                fprintf(outfile,"node extra: no\n");
	    else
	    {
                fprintf(outfile,"node extra: yes\n");
                fprintf(outfile,"af_node_type = %d\n",n_params->af_node_type);
	    }
	}
	
    fprintf(outfile,"\nGlobal index of points\n");
	next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
        fprintf(outfile,"%ld\n",Gindex(p));
    
    for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
            fprintf(outfile,"%ld\n",Gindex(p));
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
            	fprintf(outfile,"%ld\n",Gindex(p));
	    }
	}

    intfc_node_loop(intfc,n)
	{
	    p = (*n)->posn;
            fprintf(outfile,"%ld\n",Gindex(p));
	}

	fprintf(outfile,"\nGlobal index of triangles\n");
    intfc_surface_loop(intfc,s)
    {
        surf_tri_loop(*s,t)
        {
            fprintf(outfile,"%ld\n",Gindex(t));
        }
    }

	fprintf(outfile,"\nGlobal index of curves\n");
	for (c = intfc->curves; c && *c; ++c)
	    fprintf(outfile,"%d\n",Gindex(*c));

	fprintf(outfile,"\nGlobal index of surfaces\n");
	for (s = intfc->surfaces; s && *s; ++s)
	    fprintf(outfile,"%d\n",Gindex(*s));
	
    fprintf(outfile,"\nPoint periodic shift\n");
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
	{
            fprintf(outfile,"%24.18g %24.18g %24.18g",p->pshift[0],
				p->pshift[1],p->pshift[2]);
	}
	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
            fprintf(outfile,"%24.18g %24.18g %24.18g",p->pshift[0],
				p->pshift[1],p->pshift[2]);
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
            	fprintf(outfile,"%24.18g %24.18g %24.18g",p->pshift[0],
				p->pshift[1],p->pshift[2]);
	    }
	}
       
    intfc_node_loop(intfc,n)
	{
	    p = (*n)->posn;
            fprintf(outfile,"%24.18g %24.18g %24.18g",p->pshift[0],
				p->pshift[1],p->pshift[2]);
	}

    fclose(outfile);
}	/* end printAfExtraData */

void readAfExtraData(
	Front *front,
	char *restart_name)
{
	INTERFACE *intfc = front->interf;
        STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	int i,dim = intfc->dim;
	FILE *infile;
	char filename[200];
	SURFACE **s;
	CURVE **c;
	NODE **n;
	BOND *b;
    TRI *t;

	char string[100];
	long max_point_gindex = 0;
	long max_tri_gindex = 0;

        sprintf(filename,"%s-afdata",restart_name);
        infile = fopen(filename,"r");

	printf("filename = %s\n",filename);
	next_output_line_containing_string(infile,
		"Airfoil extra front state data:");

    //TODO: When done debugging can package up the functionality in here.
    //      FT_ReadFrontState(infile,front);
    
	next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
    {
        /*
        if (wave_type(hs) != ELASTIC_BOUNDARY &&
            wave_type(hs) != ELASTIC_STRING) continue;
        */

        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->vel[i]);
        fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->force[i]);
        fscanf(infile,"\n");

        FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
        fread(sl,1,front->sizest,infile);
        fread(sr,1,front->sizest,infile);
        fscanf(infile,"\n");

        /*
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf %lf\n",&sl->impulse[i],&sr->impulse[i]);
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf %lf\n",&sl->vel[i],&sr->vel[i]);
        */
    }

	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;
        p = b->start;
	    
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->vel[i]);
        fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->force[i]);
        fscanf(infile,"\n");
        
        sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
        fread(sl,1,front->sizest,infile);
        fread(sr,1,front->sizest,infile);
        fscanf(infile,"\n");
       
        /*
        fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sl->impulse[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sr->impulse[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sl->vel[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sr->vel[i]);
	    fscanf(infile,"\n");
        */
	    
        for (b = (*c)->first; b != NULL; b = b->next)
	    {
		    p = b->end;

            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
	    	fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->force[i]);
	    	fscanf(infile,"\n");

	    	sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            fread(sl,1,front->sizest,infile);
            fread(sr,1,front->sizest,infile);
            fscanf(infile,"\n");
            
            /*
            fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    	fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sl->impulse[i]);
	    	fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->impulse[i]);
	    	fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sl->vel[i]);
	    	fscanf(infile,"\n");
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->vel[i]);
	    	fscanf(infile,"\n");
            */
	    }
	}

	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;

        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->vel[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&p->force[i]);
	    fscanf(infile,"\n");
	    
        sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
        fread(sl,1,front->sizest,infile);
        fread(sr,1,front->sizest,infile);
        fscanf(infile,"\n");

        /*
        fscanf(infile,"%lf %lf",&sl->pres,&sr->pres);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sl->impulse[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sr->impulse[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sl->vel[i]);
	    fscanf(infile,"\n");
        for (i = 0; i < dim; ++i)
            fscanf(infile,"%lf ",&sr->vel[i]);
	    fscanf(infile,"\n");
        */
	}

	next_output_line_containing_string(infile,"Surface extra data:");
    intfc_surface_loop(intfc,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_STRING) continue;
    
        int num_pts;
        fgetstring(infile,"number of registered points = ");
        fscanf(infile,"%d",&num_pts);
        if (num_pts != 0)
        {
            static REGISTERED_PTS *registered_pts;
            FT_ScalarMemoryAlloc((POINTER*)&registered_pts,
                        sizeof(REGISTERED_PTS));
            FT_VectorMemoryAlloc((POINTER*)&registered_pts->global_ids,
                        num_pts,sizeof(int));
            (*s)->extra = (REGISTERED_PTS*)registered_pts;
            registered_pts->num_pts = num_pts;
            for (i = 0; i < num_pts; ++i)
                fscanf(infile,"%d",registered_pts->global_ids+i);
        }
    }

	next_output_line_containing_string(infile,"Curve extra data:");
	for (c = intfc->curves; c && *c; ++c)
	{
		if (hsbdry_type(*c) != STRING_HSBDRY) continue;

	    C_PARAMS *c_params;
	    fgetstring(infile,"curve extra:");
        fscanf(infile,"%s",string);
	    if (string[0] == 'y')
        {
            FT_ScalarMemoryAlloc((POINTER*)&c_params,sizeof(C_PARAMS));
            fgetstring(infile,"point_mass = ");
                fscanf(infile,"%lf",&c_params->point_mass);
            fgetstring(infile,"load_mass = ");
                fscanf(infile,"%lf",&c_params->load_mass);
            fgetstring(infile,"load_type = ");
                fscanf(infile,"%d",(int*)&c_params->load_type);
            fgetstring(infile,"dir = ");
                fscanf(infile,"%d",&c_params->dir);
            (*c)->extra = (POINTER)c_params;
        }

        for (b = (*c)->first; b != (*c)->last; b = b->next)
        {
            fgetstring(infile,"string point extra:");
            fscanf(infile,"%s",string);
	        if (string[0] == 'n') continue;
            
            BOND_BENDER* bond_bender;
            FT_ScalarMemoryAlloc((POINTER*)&bond_bender,sizeof(BOND_BENDER));
            fread(bond_bender,1,sizeof(BOND_BENDER),infile);
            b->end->extra = (POINTER)bond_bender;
            fscanf(infile,"\n");
        }
	}

	next_output_line_containing_string(infile,"Node extra data:");
	for (n = intfc->nodes; n && *n; ++n)
	{
	    AF_NODE_EXTRA *n_params;
	    fgetstring(infile,"node extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') 
	    {
	    	(*n)->extra = NULL;
		continue;
	    }
	    FT_ScalarMemoryAlloc((POINTER*)&n_params,sizeof(AF_NODE_EXTRA));
	    fgetstring(infile,"af_node_type =");
            fscanf(infile,"%d",(int*)&n_params->af_node_type);
	    (*n)->extra = (POINTER)n_params;
	    (*n)->size_of_extra = sizeof(AF_NODE_EXTRA);
	}
	
    if (fgetstring(infile,"Global index of points") == FUNCTION_FAILED)
	{
	    (void) printf("String \"Global index of points\" not found\n");
	    clean_up(ERROR);
	}
	
    next_point(intfc,NULL,NULL,NULL);
    while (next_point(intfc,&p,&hse,&hs))
	{
        fscanf(infile,"%ld",&Gindex(p));
	    if (max_point_gindex < Gindex(p))
            max_point_gindex = Gindex(p);
	}

	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
            fscanf(infile,"%ld",&Gindex(p));
	    if (max_point_gindex < Gindex(p))
		max_point_gindex = Gindex(p);
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
            	fscanf(infile,"%ld",&Gindex(p));
	    	if (max_point_gindex < Gindex(p))
		    max_point_gindex = Gindex(p);
	    }
	}
	
    for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
            fscanf(infile,"%ld",&Gindex(p));
	    if (max_point_gindex < Gindex(p))
		max_point_gindex = Gindex(p);
	}

	max_point_gindex++;
	pp_global_lmax(&max_point_gindex,1);
	intfc->max_point_gindex = max_point_gindex;

	if (fgetstring(infile,"Global index of triangles") == FUNCTION_FAILED)
	{
	    (void) printf("String \"Global index of triangles\" not found\n");
	    clean_up(ERROR);
	}
    
    intfc_surface_loop(intfc,s)
    {
        surf_tri_loop(*s,t)
        {
            fscanf(infile,"%ld",&Gindex(t));
            if (max_tri_gindex < Gindex(t))
                max_tri_gindex = Gindex(t);
        }
    }

	max_tri_gindex++;
	pp_global_lmax(&max_tri_gindex,1);
	intfc->max_tri_gindex = max_tri_gindex;

	if (fgetstring(infile,"Global index of curves") == FUNCTION_FAILED)
	{
	    (void) printf("String \"Global index of curves\" not found\n");
	    clean_up(ERROR);
	}
	for (c = intfc->curves; c && *c; ++c)
            fscanf(infile,"%d",&Gindex(*c));

	if (fgetstring(infile,"Global index of surfaces") == FUNCTION_FAILED)
	{
	    (void) printf("String \"Global index of surfaces\" not found\n");
	    clean_up(ERROR);
	}
	for (s = intfc->surfaces; s && *s; ++s)
            fscanf(infile,"%d",&Gindex(*s));
	if (fgetstring(infile,"Point periodic shift") == FUNCTION_FAILED)
	{
	    (void) printf("String \"Point periodic shift\" not found\n");
	    return;
	}
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
	{
            fscanf(infile,"%lf %lf %lf",p->pshift,p->pshift+1,p->pshift+2);
	}
	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
            fscanf(infile,"%lf %lf %lf",p->pshift,p->pshift+1,p->pshift+2);
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
            	fscanf(infile,"%lf %lf %lf",p->pshift,p->pshift+1,p->pshift+2);
	    }
	}
	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
            fscanf(infile,"%lf %lf %lf",p->pshift,p->pshift+1,p->pshift+2);
	}

    fclose(infile);
}	/* end readAfExtraData */

void printHyperSurfQuality(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	int dim = Dimension(intfc);
	CURVE **c,*curve;
	SURFACE **s,*surf;
	BOND *bond;
	TRI *tri;
	double max_area,min_area,max_length,min_length;
	double scaled_area,len[MAXD];
	int i;
	RECT_GRID *gr = &topological_grid(intfc);
	double *h = gr->h;

	switch (dim)
	{
	case 2:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (wave_type(*c) != ELASTIC_BOUNDARY &&
		    wave_type(*c) != ELASTIC_STRING)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic surface quality:\n");
	    (void) printf("min_length = %f\n",min_length);
	    (void) printf("max_length = %f\n",max_length);
	    (void) printf("\n\n");
	    break;
	case 3:
	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY)
		    continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\n\nElastic curve quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_length = 0.0;
	    min_length = HUGE;
	    for (c = intfc->curves; c && *c; ++c)
	    {
		if (hsbdry_type(*c) != STRING_HSBDRY) continue;
		curve = *c;
		for (bond = curve->first; bond != NULL; bond = bond->next)
		{
		    if (scaled_bond_length(bond,h,dim) > max_length)
			max_length = scaled_bond_length(bond,h,dim);
		    if (scaled_bond_length(bond,h,dim) < min_length)
			min_length = scaled_bond_length(bond,h,dim);
		}
	    }
	    (void) printf("\nElastic string quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_area = max_length = 0.0;
	    min_area = min_length = HUGE;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
		surf = *s;
		for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); 
				tri = tri->next)
		{
		    scaled_tri_params(tri,h,&scaled_area,len);
		    if (scaled_area > max_area) max_area = scaled_area;
		    if (scaled_area < min_area) min_area = scaled_area;
		    for (i = 0; i < 3; ++i)
		    {
			if (len[i] > max_length) 
			    max_length = len[i];
			if (len[i] < min_length) 
			    min_length = len[i];
		    }
		}
	    }
	    (void) printf("\nElastic surface quality:\n");
	    (void) printf("min_scaled_area = %14.10f\n",min_area);  
	    (void) printf("max_scaled_area = %14.10f\n",max_area); 
	    (void) printf("min_scaled_tri_side = %14.10f\n",sqrt(min_length));
	    (void) printf("max_scaled_tri_side = %14.10f\n",sqrt(max_length));
	    (void) printf("\n\n");
	    break;
	}
}	/* end printHyperSurfQuality */

extern void print_elastic_params(ELASTIC_SET geom_set)
{
	double *spfr;
	Front *fr = geom_set.front;

	spfr = Spfr(fr);
    for (int i = 0; i <= 3; ++i)
    {
        printf("Max front speed(%d) = %f\n",i,spfr[i]);
    }
    
    (void) printf("Input surface parameters:\n");
    (void) printf("ks = %f  m_s = %f  lambda_s = %f\n",
                    geom_set.ks,
                    geom_set.m_s,
                    geom_set.lambda_s);
    (void) printf("Input string parameters:\n");
    (void) printf("kl = %f  m_l = %f  lambda_l = %f\n",
                    geom_set.kl,
                    geom_set.m_l,
                    geom_set.lambda_l);
    (void) printf("Input gore parameters:\n");
    (void) printf("kg = %f  m_g = %f  lambda_g = %f\n",
                    geom_set.kg,
                    geom_set.m_g,
                    geom_set.lambda_g);
	(void) printf("\ndt_tol = %20.14f  dt = %20.14f\n",
                        geom_set.dt_tol,geom_set.dt);
}	/* end print_elastic_params */

void optimizeElasticMesh(
	Front *front)
{
	if (debugging("no_optimize")) return;
	if (FT_Dimension() != 3) return;

	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = computational_grid(intfc);
	boolean nothing_done;
	int i,status;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	SCALED_REDIST_PARAMS scaled_redist_params;
	int old_string_pts,new_string_pts,old_canopy_pts,new_canopy_pts;

	if (debugging("trace"))
	    (void) printf("Entering optimizeElasticMesh()\n");

    char gvdir[100];
	if (debugging("optimize_intfc"))
	{
	    (void) printf("Quality of mesh before optimization:\n");
	    printHyperSurfQuality(front);
	    (void) printf("Checking consistency of interface\n");
	    consistent_interface(front->interf);
	    (void) printf("Checking completed\n");
        sprintf(gvdir,"%s/gview-before-optimize",OutName(front));
	    gview_plot_interface(gvdir,intfc);
	}

    int num_opt_round = 0;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    if(af_params)
    {
        num_opt_round = af_params->num_opt_round;
    }
	
    scaled_redist_params.min_scaled_bond_length = 0.45/2.0;
	scaled_redist_params.max_scaled_bond_length = 1.05/2.0;

	scaled_redist_params.min_scaled_tri_area = 0.1083;
	scaled_redist_params.max_scaled_tri_area = 0.4330;
	scaled_redist_params.min_scaled_tri_area = 0.1083/2.0;
	scaled_redist_params.max_scaled_tri_area = 0.4330/2.0;
	scaled_redist_params.min_scaled_side_length = 0.45/2.0;
	scaled_redist_params.max_scaled_side_length = 1.05/2.0;
	scaled_redist_params.aspect_tol = 3.0;

	old_string_pts = old_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		old_canopy_pts += I_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		old_string_pts += I_NumOfCurvePoints(*c) - 2;

	printf("num_opt_round = %d\n\n",num_opt_round);
	
	for (i = 0; i < num_opt_round; ++i)
	{
	    status = YES;
	    if (debugging("optimize_intfc"))
		(void) printf("Optimization round %d\n",i);
	    for (c = intfc->curves; c && *c; ++c)
	    {
	    	if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != STRING_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY)
		    continue;
	    	curve = *c;
	    	nothing_done = FT_OptimizeCurveMesh(front,curve,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
	    	if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
	    	surf = *s;
	    	nothing_done = FT_OptimizeSurfMesh(front,surf,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
	    FT_ParallelExchIntfcBuffer(front);
	    if (debugging("optimize_intfc"))
	    {
		(void) printf("Quality of mesh after %d-th round:\n",i);
	    	printHyperSurfQuality(front);
		(void) printf("Checking consistency of interface\n");
		consistent_interface(front->interf);
		(void) printf("After checking\n");
	    }
	    if (status) break;
	}

	new_string_pts = new_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		new_canopy_pts += I_NumOfSurfPoints(*s);
	for (c = intfc->curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		new_string_pts += I_NumOfCurvePoints(*c) - 2;
	if (debugging("optimize_intfc"))
	{
        sprintf(gvdir,"%s/gview-after-optimize",OutName(front));
	    gview_plot_interface(gvdir,intfc);
	}
	if (debugging("trace"))
	    (void) printf("Leaving optimizeElasticMesh()\n");
}	/* end optimizeElasticMesh */

void modifyInitialization(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	boolean bifurcate_initialization;
	
	if (!CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
	{
	    fclose(infile);
	    return;
	}
	else
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
	    {
	    	fclose(infile);
		return;
	    }
        }
	bifurcate_initialization = NO;
	if (CursorAfterStringOpt(infile,
            "Enter yes to bifurcate initialization:"))
	{
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
		bifurcate_initialization = YES;
        }
    	fclose(infile);
	if (bifurcate_initialization)
	    bifurcateCanopyModification(front);
	else
	    singleCanopyModification(front);
}	/* end modifyInitialization */

void setStressColor(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	SURFACE **s;
	TRI *tri;
	POINT *p;
	double f[MAXD];
	int i,j;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double ks = af_params->ks;
	double max_color = -HUGE;
	double min_color = HUGE;
	int n,N;
	double *color;
	
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
	    {
		n++;
		naturalStressOfTri(tri,ks);
		if (max_color < tri->color)
		    max_color = tri->color;
		if (min_color > tri->color)
		    min_color = tri->color;
	    }
	}
	FT_VectorMemoryAlloc((POINTER*)&color,n,sizeof(double));
	n = 0;
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    surf_tri_loop(*s,tri)
		tri->color = log(tri->color-min_color+1);
	}
	/* Smoothing loop */
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    I_SmoothSurfColor(*s,3);
	}
}	/* end setStressColor */

void initMovieStress(
	char *inname,
	Front *front)
{
	char string[100];
	FILE *infile = fopen(inname,"r");
	front->print_gview_color = NO;
	if (CursorAfterStringOpt(infile,
            "Type y to plot surface stress: "))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		front->print_gview_color = YES;
		FT_AddVtkIntfcMovieVariable(front,"VMSTRESS");
		return;
	    }
        }
        fclose(infile);
}	/* end initMovieStress */

void poisson_ratio(
	Front *front)
{
        INTERFACE *intfc = front->interf;
        SURFACE **s;
	CURVE **c;
	BOND *b;
	double x_min, y_min, x_max, y_max;

        intfc_surface_loop(intfc,s)
        {
            if (Boundary(*s)) continue;
	    surf_pos_curve_loop(*s,c)
	    {
		b = (*c)->first;
		x_min = x_max = b->start->_coords[0];
		y_min = y_max = b->start->_coords[1];
		for (b = (*c)->first; b != NULL; b = b->next)
		{
		    if (x_min > b->end->_coords[0])
			x_min = b->end->_coords[0];

		    if (x_max < b->end->_coords[0])
			x_max = b->end->_coords[0];

		    if (y_min > b->end->_coords[1])
			y_min = b->end->_coords[1];

		    if (y_max < b->end->_coords[1])
			y_max = b->end->_coords[1];
		}
		printf("curve x-d boundary: x_min: %f\t x_max: %f\n",
					x_min,x_max);
		printf("curve y-d boundary: y_min: %f\t y_max: %f\n\n",
					y_min,y_max);
	    }
	}
}	/* poisson_ratio */

static void naturalStressOfTri(
	TRI *tri,
	double ks)
{
	double tau[3];
	double sigma[3];
	int i,j;
	double len,len0;
	double vec[3];
	double s[3],c[3];
	double b1,b2,arg,sigma1,sigma2;

	for (i = 0; i < 3; ++i)
	{
	    len0 = tri->side_length0[i];
	    for (j = 0; j < 3; ++j)
	    {
		vec[j] = Coords(Point_of_tri(tri)[(i+1)%3])[j] -
			Coords(Point_of_tri(tri)[i])[j];
	    }
	    len = Mag3d(vec);
	    tau[i] = ks*(len - len0);
	    c[i] = vec[0]/len;
	    s[i] = vec[1]/len;
	}
	if (Mag3d(tau) < 1.0e-8) // tolerance
	{
	    tri->color = 0.0;
	    return;
	}
	// Convert to Cartesian tensor
	for (i = 0; i < 3; ++i)
	{
	    sigma[i] = sqr(c[i])*tau[0] + sqr(s[i])*tau[1] + s[i]*c[i]*tau[2];
	}
	// Diagonalize the stress tensor for principal directions
	b1 = -(sigma[0] + sigma[1]);
	b2 = sigma[0]*sigma[1] - 0.25*sqr(sigma[2]);
	arg = sqr(b1) - 4.0*b2;
	sigma1 = 0.5*(-b1 + sqrt(arg));
	sigma2 = 0.5*(-b1 - sqrt(arg));
	// Use von Mises stress as a measure
	tri->color = sqrt(sqr(sigma1) + sqr(sigma2) - sigma1*sigma2);
}	/* end naturalStressOfTri */

static void singleCanopyModification(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	double disp[MAXD],center[MAXD];
	double phi,theta;
	INTERFACE *intfc = front->interf;
	double L[MAXD],U[MAXD];
	int i,dim,gmax[MAXD];

	dim = FT_Dimension();
	CursorAfterString(infile,
		"Enter yes for translation of interior interface:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] != 'y' || string[0] != 'Y')
	{
	    CursorAfterString(infile,"Enter displacement of translation:");
            fscanf(infile,"%lf %lf %lf",disp,disp+1,disp+2);
            (void) printf("%f %f %f\n",disp[0],disp[1],disp[2]);
	    I_TransInteriorIntfcPoints(intfc,disp);
	}
	CursorAfterString(infile,
		"Enter yes for rotation of interior interface:");
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
	if (string[0] != 'y' || string[0] != 'Y')
	{
	    CursorAfterString(infile,"Enter center of rotation:");
            fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
            (void) printf("%f %f %f\n",center[0],center[1],center[2]);
	    CursorAfterString(infile,"Enter azimuthal and polar angles:");
            fscanf(infile,"%lf %lf",&phi,&theta);
            (void) printf("%f %f\n",phi,theta);
	    theta *= PI/180.0;
	    phi *= PI/180.0;
	    I_SphericalRotateInteriorIntfcPoints(intfc,center,phi,theta);
	}
	if (CursorAfterStringOpt(infile,
            "Entering yes to modify computational grid:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		for (i = 0; i < dim; ++i)
        	{
	            sprintf(input_string,
				"New domain limit in %d-th dimension:",i);
	            CursorAfterString(infile,input_string);
	            fscanf(infile,"%lf %lf",&L[i],&U[i]);
	            (void) printf("%f %f\n",L[i],U[i]);
        	}
		CursorAfterString(infile,"New computational grid:");
        	for (i = 0; i < dim; ++i)
        	{
	            fscanf(infile,"%d",&gmax[i]);
		    (void) printf("%d ",gmax[i]);
        	}
        	(void) printf("\n");
		FT_ResetDomainAndGrid(front,L,U,gmax);
	    }
        }
	fclose(infile);
}	/* end singleCanopyModification */

static void bifurcateCanopyModification(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	ELASTIC_SET *parachute_set;
	int num_canopy;
	double disp[MAXD],center[MAXD];
	double phi,theta;
	INTERFACE *intfc = front->interf;
	double L[MAXD],U[MAXD];
	int i,dim,gmax[MAXD];

	dim = FT_Dimension();

	CursorAfterString(infile,
		"Enter total number of canopy sets:");
        fscanf(infile,"%d",&num_canopy);
        (void) printf("%d\n",num_canopy);
	FT_VectorMemoryAlloc((POINTER*)&parachute_set,num_canopy,
			sizeof(ELASTIC_SET));
	/* Get the original set */
	parachute_set[0].front = front;
	assembleParachuteSet(front->interf,&parachute_set[0]);
	for (i = 1; i < num_canopy; ++i)
	{
	    copyParachuteSet(parachute_set[0],&parachute_set[i]);
	}
	for (i = 0; i < num_canopy; ++i)
	{
	    sprintf(string,"For canopy set %d",i+1);
	    CursorAfterString(infile,string);
	    (void) printf("\n");

	    	CursorAfterString(infile,
		"Enter yes for translation of canopy set:");
            fscanf(infile,"%s",input_string);
            (void) printf("%s\n",input_string);
	    if (input_string[0] == 'y' || input_string[0] == 'Y')
	    {
	    	CursorAfterString(infile,"Enter displacement of translation:");
            	fscanf(infile,"%lf %lf %lf",disp,disp+1,disp+2);
            	(void) printf("%f %f %f\n",disp[0],disp[1],disp[2]);
	    }
	    CursorAfterString(infile,
		"Enter yes for rotation of canopy set:");
            fscanf(infile,"%s",input_string);
            (void) printf("%s\n",input_string);
	    if (input_string[0] == 'y' || input_string[0] == 'Y')
	    {
	    	CursorAfterString(infile,"Enter center of rotation:");
            	fscanf(infile,"%lf %lf %lf",center,center+1,center+2);
            	(void) printf("%f %f %f\n",center[0],center[1],center[2]);
	    	CursorAfterString(infile,"Enter azimuthal and polar angles:");
            	fscanf(infile,"%lf %lf",&phi,&theta);
            	(void) printf("%f %f\n",phi,theta);
	    	theta *= PI/180.0;
	    	phi *= PI/180.0;
	    	rotateParachuteSet(&parachute_set[i],center,phi,theta);
	    }
	}
	InstallNewLoadNode(front,num_canopy);

	if (CursorAfterStringOpt(infile,
            "Entering yes to modify computational grid:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
	    {
		for (i = 0; i < dim; ++i)
        	{
	            sprintf(input_string,
				"New domain limit in %d-th dimension:",i);
	            CursorAfterString(infile,input_string);
	            fscanf(infile,"%lf %lf",&L[i],&U[i]);
	            (void) printf("%f %f\n",L[i],U[i]);
        	}
		CursorAfterString(infile,"New computational grid:");
        	for (i = 0; i < dim; ++i)
        	{
	            fscanf(infile,"%d",&gmax[i]);
		    (void) printf("%d ",gmax[i]);
        	}
        	(void) printf("\n");
		FT_ResetDomainAndGrid(front,L,U,gmax);
	    }
        }
	fclose(infile);
}	/* end bifurcateCanopyModification */

static void copyParachuteSet(
	ELASTIC_SET orig_set,
	ELASTIC_SET *copy_set)
{
	int i,j,ns,nc,nn;
	INTERFACE *cur_intfc = current_interface();
	Front *front = orig_set.front;
	NODE *start,*end;
	CURVE **c,**pos_curves,**neg_curves;
	AF_NODE_EXTRA *extra;

	set_current_interface(front->interf);

	ns = copy_set->num_surfs = orig_set.num_surfs;
	nc = copy_set->num_curves = orig_set.num_curves;
	nn = copy_set->num_nodes = orig_set.num_nodes;

	for (i = 0; i < nn; ++i)
	{
	    copy_set->nodes[i] = copy_node(orig_set.nodes[i]);
	    if (is_load_node(orig_set.nodes[i]))
	    {
		copy_set->load_node = copy_set->nodes[i];
        	FT_ScalarMemoryAlloc((POINTER*)&extra,sizeof(AF_NODE_EXTRA));
		copy_set->load_node->extra = extra;
		copy_set->load_node->size_of_extra = sizeof(AF_NODE_EXTRA);
		extra->af_node_type = LOAD_NODE;
	    }
	}
	for (i = 0; i < nc; ++i)
	{
	    start = end = NULL;
	    for (j = 0; j < nn; ++j)
	    {
		if (orig_set.curves[i]->start == orig_set.nodes[j])
		    start = copy_set->nodes[j];
		if (orig_set.curves[i]->end == orig_set.nodes[j])
		    end = copy_set->nodes[j];
	    }
	    if (start == NULL || end == NULL)
	    {
		printf("Cannot find start or end node\n");
		clean_up(ERROR);
	    }
	    copy_set->curves[i] = copy_curve(orig_set.curves[i],start,end);
	}
	for (i = 0; i < ns; ++i)
	{
	    pos_curves = neg_curves = NULL;
	    for (j = 0; j < nc; ++j)
	    {
		surf_pos_curve_loop(orig_set.surfs[i],c)
		{
		    if (*c == orig_set.curves[j])
			unique_add_to_pointers(copy_set->curves[j],&pos_curves);
		}
		surf_neg_curve_loop(orig_set.surfs[i],c)
		{
		    if (*c == orig_set.curves[j])
			unique_add_to_pointers(copy_set->curves[j],&neg_curves);
		}
	    }
	    copy_set->surfs[i] = copy_surface(orig_set.surfs[i],pos_curves,
				neg_curves,YES);
	}
}	/* end copyParachuteSet */

static void rotateParachuteSet(
	ELASTIC_SET *parachute_set,
	double *center,
	double phi,
	double theta)
{
	int i;
	for (i = 0; i < parachute_set->num_surfs; ++i)
	{
	    I_SphericalRotateInteriorSurfPoints(parachute_set->surfs[i],
					center,phi,theta);
	}
	for (i = 0; i < parachute_set->num_curves; ++i)
	{
	    I_SphericalRotateInteriorCurvePoints(parachute_set->curves[i],
					center,phi,theta);
	}
	for (i = 0; i < parachute_set->num_nodes; ++i)
	{
	    I_SphericalRotatePoint(parachute_set->nodes[i]->posn,
					center,phi,theta,NO);
	}
}	/* end rotateParachuteSet */

extern int numOfMonoHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

extern int numOfGoreHsbdry(
	INTERFACE *intfc)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) nc++;
	} 
	return nc;
}	/* end numOfMonoBdry */

extern int numOfGoreNodes(
	INTERFACE *intfc)
{
	NODE **n;
	CURVE **c;
	int num_gore_nodes = 0;
	AF_NODE_EXTRA *extra;
	boolean is_string_node;

	intfc_node_loop(intfc,n)
	{
	    if ((*n)->extra == NULL)
		continue;
	    is_string_node = NO;
	    for (c = (*n)->in_curves; c && *c; ++c)
		if (hsbdry_type(*c) == STRING_HSBDRY)
		    is_string_node = YES;
	    for (c = (*n)->out_curves; c && *c; ++c)
		if (hsbdry_type(*c) == STRING_HSBDRY)
		    is_string_node = YES;
	    if (is_string_node) continue;
	    extra = (AF_NODE_EXTRA*)(*n)->extra;
	    if (extra->af_node_type == GORE_NODE)
		num_gore_nodes++;
	}
	return num_gore_nodes;
}	/* numOfGoreNodes */

extern int arrayOfMonoHsbdry(
	INTERFACE *intfc,
	CURVE **mono_curves)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY) 
	    {
		mono_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfMonoBdry */

extern int arrayOfGoreHsbdry(
	INTERFACE *intfc,
	CURVE **gore_curves)
{
	CURVE **c;
	int nc = 0;
	intfc_curve_loop(intfc,c)
	{
	    if (hsbdry_type(*c) == GORE_HSBDRY) 
	    {
		gore_curves[nc] = *c;
		nc++;
	    }
	} 
	return nc;
}	/* end arrayOfGoreBdry */

extern int getGoreNodes(
	INTERFACE *intfc,
	NODE **gore_nodes)
{
	NODE **n;
	int num_nodes = 0;

	intfc_node_loop(intfc,n)
	{
	    if (is_gore_node(*n))
		gore_nodes[num_nodes++] = *n;
	}
	return num_nodes;
}	/* getGoreNodes */

extern boolean is_bdry_node(
	NODE *node)
{
	CURVE **c;
	for (c = node->in_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == NEUMANN_HSBDRY ||
		hsbdry_type(*c) == DIRICHLET_HSBDRY ||
		hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
	    {
		return YES;
	    }
	} 
	for (c = node->out_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == NEUMANN_HSBDRY ||
		hsbdry_type(*c) == DIRICHLET_HSBDRY ||
		hsbdry_type(*c) == SUBDOMAIN_HSBDRY) 
	    {
		return YES;
	    }
	} 
	return NO;
}	/* is_bdry_node */

extern boolean is_gore_node(
	NODE *node)
{
	CURVE **c;
	AF_NODE_EXTRA *extra;

	if (node->extra == NULL)
	    return NO;
	for (c = node->in_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return NO;
	for (c = node->out_curves; c && *c; ++c)
	    if (hsbdry_type(*c) == STRING_HSBDRY)
		return NO;
	extra = (AF_NODE_EXTRA*)(node)->extra;
	if (extra->af_node_type == GORE_NODE)
	    return YES;
	else 
	    return NO;
}	/* end is_gore_node */

extern boolean is_string_node(NODE *n)
{
        AF_NODE_EXTRA *af_node_extra;
        if (n->extra == NULL) return NO;
        af_node_extra = (AF_NODE_EXTRA*)n->extra;
        if (af_node_extra->af_node_type == STRING_NODE) return YES;
        return NO;
}       /* end is_load_node */

extern boolean is_load_node(NODE *n)
{
        AF_NODE_EXTRA *af_node_extra;
        if (n->extra == NULL) return NO;
        af_node_extra = (AF_NODE_EXTRA*)n->extra;
        if (af_node_extra->af_node_type == LOAD_NODE) return YES;
        return NO;
}       /* end is_load_node */

extern boolean is_rg_string_node(NODE *n)
{
        AF_NODE_EXTRA *af_node_extra;
        if (n->extra == NULL) return NO;
        af_node_extra = (AF_NODE_EXTRA*)n->extra;
        if (af_node_extra->af_node_type == RG_STRING_NODE) return YES;
        return NO;
}       /* end is_rg_string_node */

extern boolean goreInIntfc(
	INTERFACE *intfc)
{
	NODE **n;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    if (is_gore_node(*n))
		return YES;
	}
	return NO;
}	/* end goreInIntfc */

extern double springCharTimeStep(
	Front *fr)
{
	AF_PARAMS *af_params = (AF_PARAMS*)fr->extra2;
	double dt_tol;
	dt_tol = sqrt((af_params->m_s)/(af_params->ks));
        if (af_params->m_l != 0.0 &&
            dt_tol > sqrt((af_params->m_l)/(af_params->kl)))
            dt_tol = sqrt((af_params->m_l)/(af_params->kl));
        if (af_params->m_g != 0.0 &&
            dt_tol > sqrt((af_params->m_g)/(af_params->kg)))
            dt_tol = sqrt((af_params->m_g)/(af_params->kg));
	return dt_tol;
}	/* end springCharTimeStep */


