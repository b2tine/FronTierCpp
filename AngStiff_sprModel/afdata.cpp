#include <iFluid.h>
#include <airfoil.h>

static void naturalStressOfTri(TRI*,double);

extern void printAfExtraDada(
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
	CURVE **c;
	NODE **n;
	BOND *b;

	sprintf(filename,"%s/state.ts%s",out_name,
                        right_flush(front->step,7));
#if defined(__MPI__)
        if (pp_numnodes() > 1)
            sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
        sprintf(filename,"%s-afdata",filename);
        outfile = fopen(filename,"w");

	fprintf(outfile,"\nAirfoil extra front state data:\n");

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g %24.18g\n",sl->impulse[i],
					sr->impulse[i]);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
        }
	for (c = intfc->curves; c && *c; ++c)
	{
	    b = (*c)->first;	p = b->start;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
	    fprintf(outfile,"\n");
            fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
	    for (b = (*c)->first; b != NULL; b = b->next)
	    {
		p = b->end;
	    	sl = (STATE*)left_state(p);
	    	sr = (STATE*)right_state(p);
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",p->vel[i]);
            	fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
	    	fprintf(outfile,"\n");
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    	fprintf(outfile,"\n");
            	for (i = 0; i < dim; ++i)
                    fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    	fprintf(outfile,"\n");
	    }
	}
	for (n = intfc->nodes; n && *n; ++n)
	{
	    p = (*n)->posn;
	    sl = (STATE*)left_state(p);
	    sr = (STATE*)right_state(p);
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",p->vel[i]);
            fprintf(outfile,"%24.18g %24.18g\n",sl->pres,sr->pres);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sl->impulse[i]);
	    fprintf(outfile,"\n");
            for (i = 0; i < dim; ++i)
                fprintf(outfile,"%24.18g ",sr->impulse[i]);
	    fprintf(outfile,"\n");
	}
	fprintf(outfile,"\nCurve extra data:\n");
	for (c = intfc->curves; c && *c; ++c)
	{
	    C_PARAMS *c_params = (C_PARAMS*)(*c)->extra;
	    if (c_params == NULL)
                fprintf(outfile,"curve extra: no\n");
	    else
	    {
                fprintf(outfile,"curve extra: yes\n");
                fprintf(outfile,"point_mass = %24.18g\n",c_params->point_mass);
                fprintf(outfile,"load_mass = %24.18g\n",c_params->load_mass);
                fprintf(outfile,"load_type = %d\n",c_params->load_type);
                fprintf(outfile,"dir = %d\n",c_params->dir);
	    }
	}
	fprintf(outfile,"\nNode extra data:\n");
	for (n = intfc->nodes; n && *n; ++n)
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
	fclose(outfile);
}	/* end printAfExtraDada */

extern void readAfExtraDada(
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
	CURVE **c;
	NODE **n;
	BOND *b;
	char string[100];

        sprintf(filename,"%s-afdata",restart_name);
        infile = fopen(filename,"r");

	printf("filename = %s\n",filename);
	next_output_line_containing_string(infile,
		"Airfoil extra front state data:");

	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            if (wave_type(hs) != ELASTIC_BOUNDARY) continue;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf %lf\n",&sl->impulse[i],&sr->impulse[i]);
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
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->impulse[i]);
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
            	for (i = 0; i < dim; ++i)
                    fscanf(infile,"%lf ",&sr->impulse[i]);
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
            for (i = 0; i < dim; ++i)
                fscanf(infile,"%lf ",&sr->impulse[i]);
	}
	next_output_line_containing_string(infile,"Curve extra data:");
	for (c = intfc->curves; c && *c; ++c)
	{
	    C_PARAMS *c_params;
	    fgetstring(infile,"curve extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
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
	next_output_line_containing_string(infile,"Node extra data:");
	for (n = intfc->nodes; n && *n; ++n)
	{
	    AF_NODE_EXTRA *n_params;
	    fgetstring(infile,"node extra:");
            fscanf(infile,"%s",string);
	    if (string[0] == 'n') continue;
	    FT_ScalarMemoryAlloc((POINTER*)&n_params,sizeof(AF_NODE_EXTRA));
	    fgetstring(infile,"af_node_type =");
            fscanf(infile,"%d",(int*)&n_params->af_node_type);
	    (*n)->extra = (POINTER)n_params;
	}
}	/* end readAfExtraDada */

extern void printHyperSurfQuality(
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
		if (wave_type(*c) != ELASTIC_BOUNDARY)
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
}	/* end printHyperSurfQuality */

extern void optimizeElasticMesh(
	Front *front)
{
	INTERFACE *intfc = front->interf;
	RECT_GRID *gr = computational_grid(intfc);
	boolean nothing_done;
	int i,status;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	SCALED_REDIST_PARAMS scaled_redist_params;
	int old_string_pts,new_string_pts,old_canopy_pts,new_canopy_pts;
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	int num_opt_round;

	if (debugging("trace"))
	    (void) printf("Entering optimizeElasticMesh()\n");
	if (debugging("optimize_intfc"))
	{
	    (void) printf("Quality of mesh before optimization:\n");
	    printHyperSurfQuality(front);
	    (void) printf("Checking consistency of interface\n");
	    consistent_interface(front->interf);
	    (void) printf("Checking completed\n");
	    gview_plot_interface("gview-before-optimize",intfc);
	    if (debugging("no_optimize"))
		clean_up(0);
	}
	if (debugging("no_optimize")) return;
	if (gr->dim == 2) return;

	num_opt_round = af_params->num_opt_round;
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

	printf("num_opt_round = %d\n",num_opt_round);
	num_opt_round = 20;
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
	    gview_plot_interface("gview-after-optimize",intfc);
	    //clean_up(0);
	}
	if (debugging("trace"))
	    (void) printf("Leaving optimizeElasticMesh()\n");
}	/* end optimizeElasticMesh */

extern void modifyInitialization(
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
	
	if (CursorAfterStringOpt(infile,
            "Entering yes to modify initialization:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] != 'y' && string[0] != 'Y')
		return;
        }
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
            if (string[0] != 'y' && string[0] != 'Y')
		return;
        }
	for (i = 0; i < dim; ++i)
        {
            sprintf(input_string,"New domain limit in %d-th dimension:",i);
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
}	/* end modifyInitialization */

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
	vtkPlotSurfaceStress(front);
}	/* end gviewSurfaceStress */

static void naturalStressOfTri(
	TRI *tri,
	double ks)
{
        double tau[3];
        double sigma[3];
        int i,j;
        double len0;
        double vec[3];
        double s[3],c[3];
        double b1,b2,arg,sigma1,sigma2;
        double tcoords[3][MAXD];
        double len[3];
        double A;
        double mapcof[3][3];

        for (i = 0; i < 3; ++i)
        {
            for (j = 0; j < MAXD; ++j)
            {
                vec[j] = Coords(Point_of_tri(tri)[(i+1)%3])[j] -
                        Coords(Point_of_tri(tri)[i])[j];
                tcoords[i][j] = 0;
            }
            len[i] = Mag3d(vec);
        }
        /*transform coordinates to triangle plane with point[0] as origin*/
        tcoords[1][0] = len[0];
        tcoords[2][0] = (sqr(len[2])+sqr(len[0])-sqr(len[1]))/(2*len[0]);
        tcoords[2][1] = sqrt(sqr(len[2]) - sqr(tcoords[2][0]));

        for (i = 0; i < 3; ++i)
        {
            len0 = tri->side_length0[i];
            for (j = 0; j < 3; ++j)
            {
                vec[j] = tcoords[(i+1)%3][j] -
                        tcoords[i][j];
            }
            tau[i] = ks*(len[i] - len0);
            c[i] = vec[0]/len[i];
            s[i] = vec[1]/len[i];
        }

        /* Convert to Cartesian tensor */
        A = tri_area(tri);
        mapcof[0][0] = (tcoords[2][1]-tcoords[0][1])*
			(tcoords[1][1]-tcoords[0][1])*len[1]*len[1];
        mapcof[0][1] = (tcoords[0][1]-tcoords[1][1])*
			(tcoords[2][1]-tcoords[1][1])*len[2]*len[2];
        mapcof[0][2] = (tcoords[1][1]-tcoords[2][1])*
			(tcoords[0][1]-tcoords[2][1])*len[0]*len[0];
        mapcof[1][0] = (tcoords[2][0]-tcoords[0][0])*
			(tcoords[1][0]-tcoords[0][0])*len[1]*len[1];
        mapcof[1][1] = (tcoords[0][0]-tcoords[1][0])*
			(tcoords[2][0]-tcoords[1][0])*len[2]*len[2];
        mapcof[1][2] = (tcoords[1][0]-tcoords[2][0])*
			(tcoords[0][0]-tcoords[2][0])*len[0]*len[0];

        mapcof[2][0] = ((tcoords[2][1]-tcoords[0][1])*
			(tcoords[0][0]-tcoords[1][0])+
                        (tcoords[0][0]-tcoords[2][0])*
			(tcoords[1][1]-tcoords[0][1]))*len[1]*len[1];
        mapcof[2][1] = ((tcoords[0][1]-tcoords[1][1])*
			(tcoords[1][0]-tcoords[2][0])+
                        (tcoords[1][0]-tcoords[0][0])*
			(tcoords[2][1]-tcoords[1][1]))*len[2]*len[2];
        mapcof[2][2] = ((tcoords[1][1]-tcoords[2][1])*
			(tcoords[2][0]-tcoords[0][0])+
                        (tcoords[2][0]-tcoords[1][0])*
			(tcoords[0][1]-tcoords[2][1]))*len[0]*len[0];
        for (i = 0; i < 3; ++i)
        {
            sigma[i]=mapcof[i][0]*tau[0]+mapcof[i][1]*tau[1]
					+mapcof[i][2]*tau[2];
            sigma[i] /= 4*A*A;
        }
        sigma[2] *= 0.5;

        /* Diagonalize the stress tensor for principal directions*/
        b1 = 0.5 * (sigma[0] + sigma[1]);
        b2 = sqrt(sqr(0.5*(sigma[0] - sigma[1]))+sqr(sigma[2]));
        sigma1 = b1 + b2;
        sigma2 = b1 - b2;
        /* Use von Mises stress as a measure*/
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
}	/* end vtkSurfaceStress */

extern void PRYM(
	Front *front)
{
        INTERFACE *intfc = front->interf;
        SURFACE **s;
	CURVE **c;
	BOND *b;
	double x_min[4], y_min[4], x_max[4], y_max[4];
	double pr, dx, dy;
	int i, j, k;
	FILE *infile = fopen(InName(front),"r");
	double L[MAXD],U[MAXD];

	CursorAfterString(infile,"Enter lower bounds of the rectangle:");
	fscanf(infile,"%lf %lf",&L[0],&L[1]);
	(void) printf("%f %f\n",L[0],L[1]);
	CursorAfterString(infile,"Enter upper bounds of the rectangle:");
	fscanf(infile,"%lf %lf",&U[0],&U[1]);
	(void) printf("%f %f\n",U[0],U[1]);

	double x_org_len = U[0] - L[0];
	double y_org_len = U[1] - L[1];
	double strain = 0, stress = 0;

	/* compute Poisson ratio */
        intfc_surface_loop(intfc,s)
        {
            if (Boundary(*s)) continue;
	    i = 0;
	    surf_pos_curve_loop(*s,c)
	    {
		b = (*c)->first;
		x_min[i] = x_max[i] = b->start->_coords[0];
		y_min[i] = y_max[i] = b->start->_coords[1];
		for (b = (*c)->first; b != NULL; b = b->next)
		{
		    if (x_min[i] > b->end->_coords[0])
			x_min[i] = b->end->_coords[0];

		    if (x_max[i] < b->end->_coords[0])
			x_max[i] = b->end->_coords[0];

		    if (y_min[i] > b->end->_coords[1])
			y_min[i] = b->end->_coords[1];

		    if (y_max[i] < b->end->_coords[1])
			y_max[i] = b->end->_coords[1];
		}
		printf("curve x-d boundary: x_min: %f\t x_max: %f\n",
							x_min[i],x_max[i]);
		printf("curve y-d boundary: y_min: %f\t y_max: %f\n\n",
							y_min[i],y_max[i]);
		i++;
	    }
	    dx = x_max[2] - x_max[0] - x_org_len;
	    dy = y_org_len - (y_min[1] - y_max[3]);
	    pr = dy * x_org_len / dx / y_org_len;
	    printf("dx = %f\t dy = %f\n",dx,dy);
	    printf("x_org_len = %f\t y_org_len = %f\n",x_org_len,y_org_len);
	    printf("Poisson Ratio is:\t%f\n",pr);
	    strain = dx/x_org_len;
	    printf("Strain = %f\n",strain);
	}

	/* compute Young's modules */
	TRI *tri;
	TRI **tris;
	BOND_TRI **btris;
	POINT *pt;
	int nt;
	int side = 0;
	double force[2][3], vec[3]; 
	double len = 0;

        intfc_surface_loop(intfc,s)
        {
            if (Boundary(*s)) continue;
	    surf_pos_curve_loop(*s,c)
	    {
		double x_dis=Coords((*c)->start->posn)[0] -
				Coords((*c)->end->posn)[0];
		if (fabs(x_dis) > 1E-6)
		    continue;
		int num = 0;
		force[side][0] = force[side][1] = force[side][2] = 0.0;
		for (b = (*c)->first; b != NULL; b = b->next)
		{
		    pt = b->end;
		    btris = Btris(b);
		    tri = btris[0]->tri;
		    nt = I_FirstRingTrisAroundPoint(pt,tri,&tris);  
		    for (i = 0; i < nt; i++)
		    for (j = 0; j < 3; j++)
		    {
			if (pt == Point_of_tri(tris[i])[j])
			{
			    if (is_side_bdry(tris[i],j))
				break;
			    /* tensile stiffness part */
			    len = 0.0;
			    for (k = 0; k < 3; k++)
			    {
				vec[k] = Coords(Point_of_tri(tris[i])[j])[k] -
				     Coords(Point_of_tri(tris[i])[(j+1)%3])[k];
				len += sqr(vec[k]);
			    }		
			    len = sqrt(len);
			    for (k = 0; k < 3; k++)
			    {
				vec[k] /= len;
				force[side][k] += tris[i]->k[j] * vec[k] *
					(len - tris[i]->side_length0[j]);
			    }
			    /* angular stiffness part */
			    TRI *ltri = NULL, *rtri = NULL;
			    ltri = tris[i];
			    len = 0.0;
			    for (k = 0; k < 3; k++)
			    {
				len += sqr(Coords(Point_of_tri(ltri)[j])[k]-
					Coords(Point_of_tri(ltri)[(j+2)%3])[k]);
			    }
			    len = sqrt(len);
			    for (k = 0; k < 3; k++)
			    {
				force[side][k] += ltri->gam[j] * vec[k] *
					(len - ltri->side_length0[(j+2)%3]);
			    }
			    len = 0.0;
			    for (k = 0; k < 3; k++)
			    {
				len+=sqr(Coords(Point_of_tri(ltri)[(j+1)%3])[k]-
					Coords(Point_of_tri(ltri)[(j+2)%3])[k]);
			    }
			    len = sqrt(len);
			    for (k = 0; k < 3; k++)
			    {
				force[side][k] += ltri->gam[(j+1)%3] * vec[k] *
					(len - ltri->side_length0[(j+1)%3]);
			    }

			    rtri = Tri_on_side(ltri,j);
			    int m = 0;
			    for (m = 0; m < 3; ++m)
			    {
				if (pt == Point_of_tri(rtri)[m])
				    break;
			    }
			    len = 0.0;
			    for (k = 0; k < 3; k++)
			    {
				len += sqr(Coords(Point_of_tri(rtri)[m])[k]-
					Coords(Point_of_tri(rtri)[(m+1)%3])[k]);
			    }
			    len = sqrt(len);
			    for (k = 0; k < 3; k++)
			    {
				force[side][k] += rtri->gam[m] * vec[k] *
					(len - rtri->side_length0[m]);
			    }
			    len = 0.0;
			    for (k = 0; k < 3; k++)
			    {
				len+=sqr(Coords(Point_of_tri(rtri)[(m+1)%3])[k]-
					Coords(Point_of_tri(rtri)[(m+2)%3])[k]);
			    }
			    len = sqrt(len);
			    for (k = 0; k < 3; k++)
			    {
				force[side][k] += rtri->gam[(m+1)%3] * vec[k] *
					(len - rtri->side_length0[(m+1)%3]);
			    }
			}
		    }
		}
		side++;
	    }
	}
	for (i = 0; i < 2; i++)
	{
	    printf("\nFinal forces are : %f\t %f\t %f\n",
				force[i][0],force[i][1],force[i][2]);
	}
	stress = fabs(force[0][0]) / y_org_len;
	printf("Stress = %f\n", stress);
	printf("Strain = %f\n", strain);
	printf("Youngs E = %f\n\n\n", stress / strain);
}	/* end PRYM */
