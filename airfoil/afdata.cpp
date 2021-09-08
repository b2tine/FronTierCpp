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

#include "airfoil.h"
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

        /*
        //TODO: This all gets handled in one of the user_fprint_surface
        //      functions and can be removed -- keep for now as reference.
        //
        else if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
        {
            HYPER_SURF* hs = Hyper_surf(*s);
            fprintf(outfile,"body_index = %d\n",body_index(hs));
            fprintf(outfile,"total_mass = %g\n",total_mass(hs));
            fprintf(outfile,"moment_of_inertial = %g\n",mom_inertial(hs));
            fprintf(outfile,"angular_velo = %g\n",angular_velo(hs));
            fprintf(outfile,"motion_type = %d\n",motion_type(hs));
            fprintf(outfile,"xcom = %g %g %g\n", center_of_mass(hs)[0],
                    center_of_mass(hs)[1],center_of_mass(hs)[2]);
            fprintf(outfile,"vcom = %g %g %g\n",center_of_mass_velo(hs)[0],
                    center_of_mass_velo(hs)[1],center_of_mass_velo(hs)[2]);
            fprintf(outfile,"crot = %g %g %g\n",rotation_center(hs)[0],
                    rotation_center(hs)[1],rotation_center(hs)[2]);
            fprintf(outfile,"tdir = %g %g %g\n",translation_dir(hs)[0],
                    translation_dir(hs)[1],translation_dir(hs)[2]);
            fprintf(outfile,"rotdir = %g %g %g\n",rotation_direction(hs)[0],
                    rotation_direction(hs)[1],rotation_direction(hs)[2]);
            fprintf(outfile,"pmomi = %g %g %g\n",p_mom_inertial(hs)[0],
                    p_mom_inertial(hs)[1],p_mom_inertial(hs)[1]);
            fprintf(outfile,"pangv = %g %g %g\n",p_angular_velo(hs)[0],
                    p_angular_velo(hs)[1],p_angular_velo(hs)[2]);
            fprintf(outfile,"eulerp = %g %g %g %g\n",euler_params(hs)[0],
                    euler_params(hs)[1],euler_params(hs)[2],euler_params(hs)[3]);
        }
        */
    }

    //TODO: need FINITE_STRING also ... do we still even us these below? or is this all in node extra?
    fprintf(outfile,"\nCurve extra data:\n");
    intfc_curve_loop(intfc,c)
	{
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;

        /*
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
        */

        if ((*c)->extra == nullptr)
        {
            fprintf(outfile,"curve extra: no\n");
        }
        else
        {
            fprintf(outfile,"curve extra: yes\n");
            FINITE_STRING* s_params = (FINITE_STRING*)(*c)->extra;
            fprintf(outfile,"radius = %24.18g\n",s_params->radius);
            fprintf(outfile,"density = %24.18g\n",s_params->dens);
            fprintf(outfile,"c_drag = %24.18g\n",s_params->c_drag);
            fprintf(outfile,"ampFluidFactor = %24.18g\n",s_params->ampFluidFactor);
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
	printf("filename = %s\n",filename);
    infile = fopen(filename,"r");

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


    //////////////////////////////////////////////////////////
    /*
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY &&
		wave_type(hs) != ELASTIC_STRING) 
		continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            
            //TODO: Replace below with this
            //
            //      fread(sl,1,front->sizest,outfile);
            //      fread(sr,1,front->sizest,outfile);
            //
            
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf %lf\n",&sl->impulse[i],&sr->impulse[i]);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf %lf\n",&sl->vel[i],&sr->vel[i]);

            //TODO:Keep this though
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
                
            fscanf(infile,"\n");
        }
	
    for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
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
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
	    	sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            	for (i = 0; i < dim; ++i)
                    fscanf(infile,"%lf ",&p->vel[i]);
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
	    }
	}
	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&p->vel[i]);
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
	}
    */
    ////////////////////////////////////////////////////////////
	
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
        
        /*
        //TODO: This gets handled by the user_read_print_surface()
        //      function somewhere and can be removed -- keep for
        //      now as reference.

        else if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
        {
            HYPER_SURF* hs = Hyper_surf(*s);
	        fgetstring(infile,"body_index = "); fscanf(infile,"%d",&body_index(hs));
	        fgetstring(infile,"total_mass = "); fscanf(infile,"%g",&total_mass(hs));
	        fgetstring(infile,"moment_of_inertial = "); fscanf(infile,"%g",&mom_inertial(hs));
	        fgetstring(infile,"angular_velo = "); fscanf(infile,"%g",&angular_velo(hs));
	        fgetstring(infile,"motion_type = "); fscanf(infile,"%d",&motion_type(hs));
	        fgetstring(infile,"xcom = "); fscanf(infile,"%g %g %g",&center_of_mass(hs)[0],
                    &center_of_mass(hs)[1],&center_of_mass(hs)[2]);
	        fgetstring(infile,"vcom = "); fscanf(infile,"%g %g %g",&center_of_mass_velo(hs)[0],
                    &center_of_mass_velo(hs)[1],&center_of_mass_velo(hs)[2]);
	        fgetstring(infile,"crot = "); fscanf(infile,"%g %g %g",&rotation_center(hs)[0],
                    &rotation_center(hs)[1],&rotation_center(hs)[2]);
	        fgetstring(infile,"tdir = "); fscanf(infile,"%g %g %g",&translation_dir(hs)[0],
                    &translation_dir(hs)[1],&translation_dir(hs)[2]);
	        fgetstring(infile,"rotdir = "); fscanf(infile,"%g %g %g",&rotation_direction(hs)[0],
                    &rotation_direction(hs)[1],&rotation_direction(hs)[2]);
	        fgetstring(infile,"pmomi = "); fscanf(infile,"%g %g %g",&p_mom_inertial(hs)[0],
                    &p_mom_inertial(hs)[1],&p_mom_inertial(hs)[1]);
            fgetstring(infile,"pangv = "); fscanf(infile,"%g %g %g",&p_angular_velo(hs)[0],
                    &p_angular_velo(hs)[1],&p_angular_velo(hs)[2]);
            fgetstring(infile,"eulerp = "); fscanf(infile,"%g %g %g %g",&euler_params(hs)[0],
                    &euler_params(hs)[1],&euler_params(hs)[2],&euler_params(hs)[3]);
        }
        */
    }
    
    //TODO: C_PARAMS no longer in use, should be removed.
    next_output_line_containing_string(infile,"Curve extra data:");
	for (c = intfc->curves; c && *c; ++c)
	{
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;

        /*
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
        */

        FINITE_STRING* s_params;
	    fgetstring(infile,"curve extra:");
        fscanf(infile,"%s",string);
	    if (string[0] == 'y')
        {
            FT_ScalarMemoryAlloc((POINTER*)&s_params,sizeof(FINITE_STRING));
            fgetstring(infile,"radius = ");
            fscanf(infile,"%lf",&s_params->radius);
            fgetstring(infile,"density = ");
            fscanf(infile,"%lf",&s_params->dens);
            fgetstring(infile,"c_drag = ");
            fscanf(infile,"%lf",&s_params->c_drag);
            fgetstring(infile,"ampFluidFactor = ");
            fscanf(infile,"%lf",&s_params->ampFluidFactor);
            (*c)->extra = (POINTER)s_params;
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

void clearRegisteredPoints(Front* front)
{
    FILE *infile = fopen(InName(front),"r");
    char string[100];

    if (CursorAfterStringOpt(infile,"Enter yes to keep registered points: "))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
        {
            fclose(infile);
            return;
        }
    }
    fclose(infile);

    SURFACE** s;
    intfc_surface_loop(front->interf,s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_STRING) continue;
        
        if ((*s)->extra)
        {
            REGISTERED_PTS* registered_pts = (REGISTERED_PTS*)(*s)->extra;
            FT_FreeThese(2,registered_pts->global_ids,registered_pts);
            (*s)->extra = nullptr;
        }
    }
}

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
		if (hsbdry_type(*c) != STRING_HSBDRY)
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
	    (void) printf("\nElastic string quality:\n");
	    (void) printf("min_scaled_length = %14.10f\n",min_length);
	    (void) printf("max_scaled_length = %14.10f\n",max_length);

	    max_area = max_length = 0.0;
	    min_area = min_length = HUGE;
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
		if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
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
    fflush(stdout);
}	/* end printHyperSurfQuality */

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

void optimizeElasticStrings(
	Front *front)
{
	if (debugging("no_optimize")) return;
	if (FT_Dimension() != 3) return;

	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = computational_grid(intfc);
	boolean nothing_done;
	int i,status;
	CURVE **c,*curve;
	//SURFACE **s,*surf;
	SCALED_REDIST_PARAMS scaled_redist_params;
	int old_string_pts,new_string_pts;
    //int old_canopy_pts,new_canopy_pts;

	if (debugging("trace"))
	    (void) printf("Entering optimizeElasticCurves()\n");

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

    /*
	//scaled_redist_params.min_scaled_tri_area = 0.1083;
	//scaled_redist_params.max_scaled_tri_area = 0.4330;
	scaled_redist_params.min_scaled_tri_area = 0.1083/2.0;
	scaled_redist_params.max_scaled_tri_area = 0.4330/2.0;
	
    scaled_redist_params.min_scaled_side_length = 0.45/2.0;
	scaled_redist_params.max_scaled_side_length = 1.05/2.0;
    */

    scaled_redist_params.aspect_tol = 3.0;

    /*
    old_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		old_canopy_pts += I_NumOfSurfPoints(*s);
    */

	old_string_pts = 0;
	for (c = intfc->curves; c && *c; ++c)
    {
	    if (hsbdry_type(*c) == STRING_HSBDRY)
            old_string_pts += I_NumOfCurvePoints(*c) - 2;
    }
	printf("num_opt_round = %d\n\n",num_opt_round);
	
	for (i = 0; i < num_opt_round; ++i)
	{
	    status = YES;
	    if (debugging("optimize_intfc"))
		(void) printf("Optimization round %d\n",i);
	    for (c = intfc->curves; c && *c; ++c)
	    {
            /*
	    	if (hsbdry_type(*c) != MONO_COMP_HSBDRY &&
		    hsbdry_type(*c) != STRING_HSBDRY &&
		    hsbdry_type(*c) != GORE_HSBDRY) continue;
            */
            if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	    	curve = *c;
	    	nothing_done = FT_OptimizeCurveMesh(front,curve,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }

        /*
	    for (s = intfc->surfaces; s && *s; ++s)
	    {
	    	if (wave_type(*s) != ELASTIC_BOUNDARY)
		    continue;
	    	surf = *s;
	    	nothing_done = FT_OptimizeSurfMesh(front,surf,
				scaled_redist_params);
	    	status *= (int)nothing_done;
	    }
        */

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

    /*
    new_canopy_pts = 0;
	for (s = intfc->surfaces; s && *s; ++s)
	    if (wave_type(*s) == ELASTIC_BOUNDARY)
		new_canopy_pts += I_NumOfSurfPoints(*s);
    */
	
	new_string_pts = 0;
    for (c = intfc->curves; c && *c; ++c)
    {
	    if (hsbdry_type(*c) == STRING_HSBDRY)
            new_string_pts += I_NumOfCurvePoints(*c) - 2;
    }

    if (debugging("optimize_intfc"))
	{
        sprintf(gvdir,"%s/gview-after-optimize",OutName(front));
	    gview_plot_interface(gvdir,intfc);
	}
	if (debugging("trace"))
	    (void) printf("Leaving optimizeElasticCurves()\n");
}	/* end optimizeElasticStrings */

void modifyInitialization(
	Front *front)
{
	char *inname = InName(front);
	FILE *infile = fopen(inname,"r");
	char string[200],input_string[200];
	boolean bifurcate_initialization;
	
	if (!CursorAfterStringOpt(infile,"Enter yes to modify initialization:"))
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
	if (CursorAfterStringOpt(infile,"Enter yes to bifurcate initialization:"))
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

    /*
	if (CursorAfterStringOpt(infile,"Enter yes to modify RGB motion"))
    {
        //rgb_modification(front);
    }
    */

}	/* end modifyInitialization */

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
    if (CursorAfterStringOpt(infile,
            "Enter yes to modify interior interface:"))
    {
        fscanf(infile,"%s",string);
        (void) printf("%s\n",string);
        if (string[0] == 'y' || string[0] == 'Y')
        {
            CursorAfterString(infile,
                "Enter yes for translation of interior interface:");
                fscanf(infile,"%s",string);
                (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
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
            if (string[0] == 'y' || string[0] == 'Y')
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
        }
    }
	
    if (CursorAfterStringOpt(infile,
            "Enter yes to modify computational grid:"))
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
            //TODO: This never gets an array to plot,
            //      should have a third argument like
            //      FT_AddVtkScalarMovieVariable() and
            //      FT_AddVtkVectorMovieVariable() do.
            FT_AddVtkIntfcMovieVariable(front,"VMSTRESS");
            return;
        }
    }
    fclose(infile);
}	/* end initMovieStress */

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
        {
		    tri->color = log(tri->color - min_color + 1);
        }
	}
	
    /* Smoothing loop */
	intfc_surface_loop(intfc,s)
	{
	    if (Boundary(*s)) continue;
	    I_SmoothSurfColor(*s,3);
	}
}	/* end setStressColor */

static void naturalStressOfTri(
	TRI *tri,
	double ks)
{
	double tau[3];
	double sigma[3];
	double len,len0;
	double vec[3];
	double s[3],c[3];
	double b1,b2,arg,sigma1,sigma2;
	int i,j;

	for (i = 0; i < 3; ++i)
	{
	    len0 = tri->side_length0[i];
	    for (j = 0; j < 3; ++j)
	    {
            vec[j] = Coords(Point_of_tri(tri)[(i+1)%3])[j]
                      - Coords(Point_of_tri(tri)[i])[j];
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

extern void vtkPlotSurfaceStress(
    Front *front)
{
    char *outname = OutName(front);
    INTERFACE *intfc = front->interf;
    SURFACE **s;
    TRI *tri;
    POINT *p;
    char dirname[200],fname[200];
    int i,j;
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    int n,N;
    double *color;
    FILE *vfile;
    int num_tri;

    n = 0;
    sprintf(dirname,"%s/%s%s",outname,"vtk.ts",
            right_flush(front->step,7));
    if (!create_directory(dirname,NO))
        {
            printf("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
        }
    sprintf(fname,"%s/%s",dirname,"stress.vtk");
    vfile = fopen(fname,"w");
    fprintf(vfile,"# vtk DataFile Version 3.0\n");
        fprintf(vfile,"Surface stress\n");
        fprintf(vfile,"ASCII\n");
        fprintf(vfile,"DATASET UNSTRUCTURED_GRID\n");

    num_tri = 0;

    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        num_tri += (*s)->num_tri;
    }
    fprintf(vfile,"POINTS %d double\n", 3*num_tri);
    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
        for (i = 0; i < 3; ++i)
        {
            p = Point_of_tri(tri)[i];
            fprintf(vfile,"%f %f %f\n",Coords(p)[0],Coords(p)[1],
                        Coords(p)[2]);
        }
        }
    }
    fprintf(vfile,"CELLS %i %i\n",num_tri,4*num_tri);
    n = 0;
    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
        fprintf(vfile,"3 %i %i %i\n",3*n,3*n+1,3*n+2);
        n++;
        }
    }
    fprintf(vfile, "CELL_TYPES %i\n",num_tri);
        intfc_surface_loop(intfc,s)
    {
            if (Boundary(*s)) continue;
            surf_tri_loop(*s,tri)
            {
                fprintf(vfile,"5\n");
            }
        }

    fprintf(vfile, "CELL_DATA %i\n", num_tri);
        fprintf(vfile, "SCALARS von_Mises_stress double\n");
        fprintf(vfile, "LOOKUP_TABLE default\n");
    intfc_surface_loop(intfc,s)
    {
        if (Boundary(*s)) continue;
        surf_tri_loop(*s,tri)
        {
        fprintf(vfile,"%f\n",tri->color);
        }
    }
}   /* end vtkPlotSurfaceStress */

extern void gviewSurfaceStress(
    Front *front)
{
    char *outname = OutName(front);
    INTERFACE *intfc = front->interf;
    SURFACE **s;
    TRI *tri;
    POINT *p;
    double f[MAXD];
    char dirname[200];
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
    sprintf(dirname,"%s/gview/%s-ts%s",outname,"gv.stress",
            right_flush(front->step,7));
    gview_plot_color_scaled_interface(dirname,intfc);
}   /* end gviewSurfaceStress */

//TODO: This clearly does not compute the poisson_ratio -- it doesn't do anything.
//      Also should be moved to another file if implemented.
void poisson_ratio(Front *front)
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

