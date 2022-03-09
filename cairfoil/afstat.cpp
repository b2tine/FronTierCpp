#include "airfoil.h"


//TODO: FUNCTION FOR DRAG COEFFICIENT OF PARACHUTE SYSTEM

static void NEW_print_airfoil_stat3d(Front*,char*);

static void print_airfoil_stat3d(Front*,char*);
static void print_airfoil_stat3d_1(Front*,char*);
static void print_airfoil_stat3d_2(Front*,char*);

static void print_strings(Front *,char *);
static void print_drag3d(Front *,char *);
static void print_rgb3d(Front *,char *);

static void print_airfoil_stat2d(Front*,char*);
static void print_airfoil_stat2d_1(Front*,char*);
static void print_airfoil_stat2d_2(Front*,char*);

static void record_stretching_length(SURFACE*,char*,double);



extern void print_airfoil_stat(
	Front *front,
	char *out_name)
{
    int dim = FT_Dimension();
    if (dim != 3) return;
    
    if (!FT_FrontContainWaveType(front,ELASTIC_BOUNDARY) &&
        !FT_FrontContainHsbdryType(front,STRING_HSBDRY)) return;

	INTERFACE* save_intfc = front->interf; 
	
	INTERFACE* elastic_intfc = nullptr;
    if (pp_numnodes() > 1)
	{
	    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	    int owner[MAXD] = {0};
	    int owner_id = af_params->node_id[0];

        int w_type[4] = {ELASTIC_BOUNDARY,ELASTIC_BAND_BOUNDARY,
            MOVABLE_BODY_BOUNDARY,NEUMANN_BOUNDARY};
        elastic_intfc = collect_hyper_surfaces(front,owner,w_type,4);
	    
        collectNodeExtra(front,elastic_intfc,owner_id);
	    front->interf = elastic_intfc;
	}	
	
	if (pp_mynode() == 0)
	{	
        print_airfoil_stat3d(front,out_name);
	}

	front->interf = save_intfc;
	
    if (elastic_intfc)
    {
        delete_interface(elastic_intfc);
    }
}	/* end print_airfoil_stat */

static void print_airfoil_stat3d(
	Front *front,
	char *out_name)
{
    //TODO: CONSOLIDATE INTO SINGLE PRINTING FUNCTION!

	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;

	switch (af_params->spring_model)
	{
	case MODEL1:
	    print_airfoil_stat3d_1(front,out_name);
	    break;
	case MODEL2:
	    print_airfoil_stat3d_2(front,out_name);//default
	    break;
	case MODEL3:
	    printf("print_airfoil_stat3d_12() not implemented!\n");
        break;
	default:
        break;
	}

	print_strings(front,out_name);
	
    print_drag3d(front,out_name);
    
    if (pp_numnodes() == 1)
    {
        print_rgb3d(front,out_name);
    }
}	/* end print_airfoil_stat3d */

/*
//TODO: WRITE THIS CONSOLIDATED VERSION
static void NEW_print_airfoil_stat3d(
	Front *front,
	char *out_name)
{
    /////////////////////////////////////////////////
	
	int j,k,nc,dim = intfc->dim;
	double payload = af_params->payload;
	double *g = af_params->gravity;
	double kg,m_g, x_sqr, side_length, vect[3];
	STATE *st;
	
    static int np,ip;
	static POINT **pts;
	static double p0[MAXD];
	
    POINT *psample;
	
    int Gmax,Posn_max;
    /////////////////////////////////////////////////
	


	INTERFACE *intfc = front->interf;
	NODE **n,*node;
	CURVE **c,*curve;
	BOND *b;
	TRI *tri;
	POINT *p;

    double esk = 0.0; //spring kinetic energy
    double esp = 0.0;
    double epi = 0.0; //spring potential energy
    double epb = 0.0;
    double egp = 0.0; //external potential energy (gravity)
    double exk = 0.0; //external kinetic energy
    double enk = 0.0; //total kinetic energy
    double canpy_area = 0.0;

    //TODO: ITERATE OVER EDGES LIKE IN STRAIN LIMITING CODE!
	SURFACE **s,*surf;
    for (s = intfc->surfaces; s && *s; ++s)
    {
        if (wave_type(*s) != ELASTIC_BOUNDARY) continue;

        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
        {
            canopy_area += tri_area(tri);
            for (j = 0; j < 3; ++j)
            {
                side_length = separation(Point_of_tri(tri)[j],Point_of_tri(tri)[(j+1)%3],3);
                x_diff = side_length - tri->side_length0[j];
                if (!is_side_bdry(tri,j))
                {
                    epi += 0.5*ks*sqr(x_diff);
                }



                double side_length = separation(Point_of_tri(tri)[j], Point_of_tri(tri)[(j+1)%3],3);

                x_sqr = 0.0;
                for (k = 0; k < 3; ++k)
                {
                    vect[k] = Coords(Point_of_tri(tri)[(j+1)%3])[k]
                        - Coords(Point_of_tri(tri)[j])[k] - tri->side_length0[j]*tri->side_dir0[j][k];

                    x_sqr += sqr(vect[k]);
                }

                if (!is_side_bdry(tri,j))
                {
                    epi += 0.5*ks*x_sqr;
                }
            }
        }

    }


    //TODO: MORE







    //When everything computed write to file
    char fname[256];
	static FILE *eskfile,*espfile,*egpfile,*efile,*exkfile,*enkfile;
	static FILE *afile,*sfile,*pfile,*vfile;
	    //static FILE *xcom_file,*vcom_file;
	    //static FILE *samplex,*sampley,*samplez;
	    //static FILE *vmaxfile;

    
    static boolean first = YES;
	if (first)
	{
	    first = NO;
        
	    sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"w");
	    sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"w");
	    sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"w");
	    sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"w");
	    sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"w");
	    sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"w");
	    
        sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"w");
	        //sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"w");
            
            //sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"w");
            //sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"w");
            //sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"w");
            //sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"w");
        

        fprintf(eskfile,"\"Spring Kinetic Energy vs. Time\"\n" "color=blue\n""thickness=1.5\n");
        fprintf(espfile,"\"Spring Potential Energy vs. Time\"\n" "color=red\n""thickness=1.5\n");
        fprintf(exkfile,"\"External Kinetic Energy vs. Time\"\n" "color=green\n""thickness=1.5\n");
        fprintf(egpfile,"\"External Potential Energy vs. Time\"\n" "color=orange\n""thickness=1.5\n");
        fprintf(enkfile,"\"Kinetic Energy vs. Time\"\n" "color=yellow\n""thickness=1.5\n");
        fprintf(efile,"\"Total Energy vs. Time\"\n" "color=navy\n""thickness=1.5\n");

        fprintf(afile,"\"Canopy Area vs. Time\"\n""color=blue\n" "thickness=1.5\n");
            //fprintf(sfile,"\"String length vs. Time\"\n""color=blue\n" "thickness=1.5\n"); //NOTE: HANDLED BY print_strings()
            
            //fprintf(pfile,"\"Payload Height vs. Time\"\n""color=blue\n" "thickness=1.5\n");
            //fprintf(vfile,"\"Payload Velocity vs. Time\"\n""color=blue\n" "thickness=1.5\n");
	        //sprintf(fname,"%s/vmax.xg",out_name); vmaxfile = fopen(fname,"w");
        
            //fprintf(xcom_file,"\"COM vs. Time\"\n""color=blue\n" "thickness=1.5\n");
            //fprintf(vcom_file,"\"V-COM vs. Time\"\n""color=blue\n" "thickness=1.5\n");
        
            //fprintf(samplex,"\"x-coords vs. Time\"\n""color=blue\n" "thickness=1.5\n");
            //fprintf(sampley,"\"y-coords vs. Time\"\n""color=red\n" "thickness=1.5\n");
            //fprintf(samplez,"\"z-coords vs. Time\"\n""color=yellow\n" "thickness=1.5\n");
    }
    else
    {
        sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"a");
        sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"a");
        sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"a");
        sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"a");
        sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"a");
        sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"a");
        
        sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"a");
            //sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"w");
            
            //sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"w");
            //sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"w");
            //sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"w");
            //sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"w");
    }


	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    double ks = af_params->ks;
    double kl = af_params->kl;
    double m_s = af_params->m_s;
    double m_l = af_params->m_l;

	double x_diff,side_length;





	fprintf(eskfile,"%16.12f  %16.12f\n", front->time, esk); fclose(eskfile);
    fprintf(espfile,"%16.12f  %16.12f\n", front->time, esp); fclose(espfile);
    fprintf(egpfile,"%16.12f  %16.12f\n", front->time, egp); fclose(egpfile);
    fprintf(exkfile,"%16.12f  %16.12f\n", front->time, exk); fclose(exkfile);
    fprintf(enkfile,"%16.12f  %16.12f\n", front->time, enk); fclose(enkfile);
    fprintf(efile,"%16.12f  %16.12f\n", front->time, esp + egp + enk); fclose(efile);

    fprintf(afile,"%16.12f  %16.12f\n", front->time, cnp_area); fclose(afile);
        //fprintf(sfile,"%16.12f  %16.12f\n", front->time, str_length); fclose(sfile);

        //fprintf(pfile,"%16.12f  %16.12f\n", front->time, pz); fclose(pfile);
        //fprintf(vfile,"%16.12f  %16.12f\n", front->time, pv); fclose(vfile);

        //fprintf(xcom_file,"%16.12f  %16.12f\n",front->time,zcom); fclose(xcom_file);
        //fprintf(vcom_file,"%16.12f  %16.12f\n",front->time,vcom); fclose(vcom_file);

}
*/

static void print_airfoil_stat3d_1(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
	NODE **n,*node;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	BOND *b;
	TRI *tri;
	POINT *p;
	static FILE *eskfile,*espfile,*egpfile,*efile,*exkfile,*enkfile;
	static FILE *vmaxfile;
	static FILE *afile,*sfile,*pfile,*vfile;
	static FILE *xcom_file,*vcom_file;
	static FILE *samplex,*sampley,*samplez;
	char fname[256];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double esk,esp,epi,epb,egp,exk,enk;
	double ks,m_s,kl,m_l,x_diff,side_length;
	int j,k,nc,dim = intfc->dim;
	double cnp_area,str_length,pz,pv;
	double zcom,vcom;
	double payload = af_params->payload;
	double *g = af_params->gravity;
	STATE *st;
	static int np,ip;
	static POINT **pts;
	POINT *psample;
	static double p0[MAXD];
	double vmax = -HUGE;
	int Gmax,Posn_max;
	static FILE *gfile;
	POINTER obj_max;

	static boolean first = YES;
	if (first)
    {
	    sprintf(fname,"%s/max_index.dat",out_name); gfile = fopen(fname,"w");
	    sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"w");
	    sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"w");
	    sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"w");
	    sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"w");
	    sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"w");
	    sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"w");
	    sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"w");
	    sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"w");
	    sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"w");
	    sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"w");
	    sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"w");
	    sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"w");
	    sprintf(fname,"%s/samplex.xg",out_name); samplex = fopen(fname,"w");
	    sprintf(fname,"%s/sampley.xg",out_name); sampley = fopen(fname,"w");
	    sprintf(fname,"%s/samplez.xg",out_name); samplez = fopen(fname,"w");
	    sprintf(fname,"%s/vmax.xg",out_name); vmaxfile = fopen(fname,"w");
            
        fprintf(eskfile,"!Spr-kinetic energy vs. time\n color=blue\n");
        fprintf(espfile,"!Spr-potentl energy vs. time\n" "color=red\n""thickness=1.5\n");
        fprintf(exkfile,"!Ext-kinetic energy vs. time\n" "color=green\n""thickness=1.5\n");
        fprintf(egpfile,"!Ext-potential energy vs. time\n" "color=orange\n""thickness=1.5\n");
        fprintf(enkfile,"!Kinetic energy vs. time\n" "color=yellow\n""thickness=1.5\n");
        fprintf(efile,"!Total energy vs. time\n" "color=navy\n""thickness=1.5\n");

        fprintf(afile,"!Canopy area vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(sfile,"!String length vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(pfile,"!Payload hight vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(vfile,"!Payload velo vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(xcom_file,"!COM vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(vcom_file,"!V-COM vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(samplex,"!x-coords vs. time\n""color=blue\n" "thickness=1.5\n");
        fprintf(sampley,"!y-coords vs. time\n""color=red\n" "thickness=1.5\n");
        fprintf(samplez,"!z-coords vs. time\n""color=yellow\n" "thickness=1.5\n");
    }
    else
    {
	    sprintf(fname,"%s/max_index.dat",out_name); gfile = fopen(fname,"a");
	    sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"a");
	    sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"a");
	    sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"a");
	    sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"a");
	    sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"a");
	    sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"a");
	    sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"a");
	    sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"a");
	    sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"a");
	    sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"a");
	    sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"a");
	    sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"a");
	    sprintf(fname,"%s/samplex.xg",out_name); samplex = fopen(fname,"a");
	    sprintf(fname,"%s/sampley.xg",out_name); sampley = fopen(fname,"a");
	    sprintf(fname,"%s/samplez.xg",out_name); samplez = fopen(fname,"a");
	    sprintf(fname,"%s/vmax.xg",out_name); vmaxfile = fopen(fname,"a");
    }

	ks = af_params->ks;
    m_s = af_params->m_s;

	esk = esp = epi = epb = egp = exk = enk = 0.0;
	cnp_area = 0.0;
	surf = NULL;
	psample = NULL;

    for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY) continue;
	    
        surf = *s;
	    zcom = center_of_mass(Hyper_surf(surf))[2];
	    vcom = center_of_mass_velo(Hyper_surf(surf))[2];
	    
        if (first)
	    {
    		np = I_NumOfSurfPoints(surf);
	    	ip = np/2;
		    FT_VectorMemoryAlloc((POINTER*)&pts,np,sizeof(POINT*));
	    }
	    else if (I_NumOfSurfPoints(surf) > np) 
	    {
            np = I_NumOfSurfPoints(surf);
            FT_FreeThese(1, pts);
            FT_VectorMemoryAlloc((POINTER*)&pts,np,sizeof(POINT*));
        }

	    I_ArrayOfSurfPoints(surf,pts);
	    psample = pts[ip];
	    for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
        {
            cnp_area += tri_area(tri);
            for (j = 0; j < 3; ++j)
            {
                side_length = separation(Point_of_tri(tri)[j],Point_of_tri(tri)[(j+1)%3],3);
                x_diff = side_length - tri->side_length0[j];
                if (!is_side_bdry(tri,j))
                {
                    epi += 0.5*ks*sqr(x_diff);
                }
            }
        }

	    unsort_surf_point(surf);
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
        {
            for (j = 0; j < 3; ++j)
            {
                p = Point_of_tri(tri)[j];
                if (sorted(p) || Boundary_point(p)) continue;
                
                for (k = 0; k < dim; ++k)
                {
                    if (fabs(p->vel[k]) > vmax)
                    {
                        vmax = fabs(p->vel[k]);
                        Gmax = Gindex(p);
                        Posn_max = 0;
                        obj_max = (POINTER)surf;
                    }
                
                    esk += 0.5*m_s*sqr(p->vel[k]);
                    egp += -g[k]*m_s*Coords(p)[k];
                
                    st = (STATE*)left_state(p);
                    exk += 0.5*m_s*sqr(st->impulse[k]);
                    enk += 0.5*m_s*sqr(p->vel[k] + st->impulse[k]);
                }
                sorted(p) = YES;
            }
        }
	}

    /*
	if (surf != NULL)
    {
        //TODO: rewrite this function
	    record_stretching_length(surf,out_name,front->time);
    }
    */

	epi *= 0.5;	//Each side is counted twice

	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == STRING_HSBDRY)
	    {
		    kl = af_params->kl;
        	m_l = af_params->m_l;
	    }
	    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    {
		    kl = af_params->ks;
        	m_l = af_params->m_s;
	    }
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
        {
            kl = af_params->kg;
            m_l = af_params->m_g;
        }
	    else
        {
		    continue;
        }
	    
        curve = *c;
	    for (b = curve->first; b != NULL; b = b->next)
        {
            x_diff = bond_length(b) - bond_length0(b);
            epb += 0.5*kl*sqr(x_diff);
            if (b != curve->last)
            {
                for (k = 0; k < dim; ++k)
                {
                    if (fabs(b->end->vel[k]) > vmax)
                    {
                        vmax = fabs(b->end->vel[k]);
                        Gmax = Gindex(b->end);
                        Posn_max = 1;
                        obj_max = (POINTER)curve;
                    }
                
                    esk += 0.5*m_l*sqr(b->end->vel[k]);
                    egp += -g[k]*m_l*Coords(b->end)[k];
                
                    st = (STATE*)left_state(b->end);
                    exk += 0.5*m_l*sqr(st->impulse[k]);
                    enk += 0.5*m_l*sqr(b->end->vel[k] + st->impulse[k]);
                }
            }
        }
	}

	esp = epi + epb;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    node = *n;
	    if (is_bdry_node(node)) continue;
	
        //TODO: Write block for RG_STRING_NODEs
        if (is_load_node(node)) 
	    {
	    	STATE *sl;
	    	pz = Coords(node->posn)[2];
	    	sl = (STATE*)left_state(node->posn);
	    	pv = sl->vel[2];
	    	for (k = 0; k < dim; ++k)
	    	{
                egp += -g[k]*payload*Coords(node->posn)[k];
                st = (STATE*)left_state(node->posn);
                exk += 0.5*payload*sqr(st->impulse[k]);
                enk += 0.5*payload*sqr(node->posn->vel[k] + st->impulse[k]);
	    	}
	    }
	    else
	    {
            if (is_gore_node(node))
            {
                m_l = af_params->m_g;
            }
            else
            {
                m_l = af_params->m_s;
            }

            for (k = 0; k < dim; ++k)
            {
                if (fabs(node->posn->vel[k]) > vmax)
                {
                    vmax = fabs(node->posn->vel[k]);
                    Gmax = Gindex(node->posn);
                    Posn_max = 2;
                    obj_max = (POINTER)node;
                }
                
                esk += 0.5*m_l*sqr(node->posn->vel[k]);
                egp += -g[k]*m_l*Coords(node->posn)[k];

                st = (STATE*)left_state(node->posn);
                exk += 0.5*m_l*sqr(st->impulse[k]);
                enk += 0.5*m_l*sqr(node->posn->vel[k] + st->impulse[k]);
            }
        }
	}

    //TODO: Monitor enk (total kinetic energy) to detect
    //      unphysical configurations of the fabric interface.
    //      Failure of the fabric solver is nearly always
    //      preceded by a spike in the spring system kinetic
    //      energy and a rapid increase in the max speed of
    //      the canopy points

	nc = 0;
    str_length = 0.0;
	for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	    str_length += curve_length(*c);
	    nc++;
	}

	if (nc != 0) str_length /= (double)nc;
	
    if (first)
	{
	    if (psample != NULL)
        {
	    	for (k = 0; k < dim; ++k)
            {
                p0[k] = Coords(psample)[k];
            }
        }
	    first = NO;
	}

	fprintf(eskfile,"%16.12f  %16.12f\n",front->time,esk); fclose(eskfile);
    fprintf(espfile,"%16.12f  %16.12f\n",front->time,esp); fclose(espfile);
    fprintf(egpfile,"%16.12f  %16.12f\n",front->time,egp); fclose(egpfile);
    fprintf(exkfile,"%16.12f  %16.12f\n",front->time,exk); fclose(exkfile);
    fprintf(enkfile,"%16.12f  %16.12f\n",front->time,enk); fclose(enkfile);
    fprintf(efile,"%16.12f  %16.12f\n",front->time,esp+egp+enk); fclose(efile);

    fprintf(afile,"%16.12f  %16.12f\n",front->time,cnp_area); fclose(afile);
    fprintf(sfile,"%16.12f  %16.12f\n",front->time,str_length); fclose(sfile);

    fprintf(pfile,"%16.12f  %16.12f\n",front->time,pz); fclose(pfile);
    fprintf(vfile,"%16.12f  %16.12f\n",front->time,pv); fclose(vfile);
    fprintf(vmaxfile,"%16.12f  %16.12f\n",front->time,vmax); fclose(vmaxfile);
	
    
    fprintf(gfile,"Max Gindex %d",Gmax);
	if (Posn_max == 0) 
	{
	    surf = (SURFACE*)obj_max;
	    fprintf(gfile," on surface type:");
	    fprintf(gfile," %s\n",f_wave_type_as_string(wave_type(surf)));
	}
	else if (Posn_max == 1) 
	{
	    curve = (CURVE*)obj_max;
	    fprintf(gfile," on curve type:");
	    fprintf(gfile," %s\n",f_hsbdry_type_as_string(hsbdry_type(curve)));
	}
	else if (Posn_max == 2) 
	{
	    node = (NODE*)obj_max;
	    fprintf(gfile," on node type:\n");
	    if (is_gore_node(node))
    		fprintf(gfile," GORE_NODE\n");
	    else if (is_load_node(node))
	    	fprintf(gfile," LOAD_NODE\n");
	    else
		    fprintf(gfile," Other NODE\n");
	}
	fflush(gfile);

    fprintf(xcom_file,"%16.12f  %16.12f\n",front->time,zcom); fclose(xcom_file);
    fprintf(vcom_file,"%16.12f  %16.12f\n",front->time,vcom); fclose(vcom_file);

    if (psample != NULL)
	{
        fprintf(samplex,"%16.12f  %16.12f\n",front->time,Coords(psample)[0] - p0[0]);
        fprintf(sampley,"%16.12f  %16.12f\n",front->time,Coords(psample)[1] - p0[1]);
        fprintf(samplez,"%16.12f  %16.12f\n",front->time,Coords(psample)[2] - p0[2]);
	}
	fclose(samplex);
	fclose(sampley);
	fclose(samplez);
}	/* end print_airfoil_stat3d_1 */

//TODO: copy in version from fabric directory when it is stable
static void print_airfoil_stat3d_2(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
	NODE **n,*node;
	CURVE **c,*curve;
	SURFACE **s,*surf;
	BOND *b;
	TRI *tri;
	POINT *p;
	
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double esk,esp,epi,epb,egp,exk,enk;
	double ks,m_s,kl,m_l,kg,m_g,x_sqr,side_length,vect[3];
	int j,k,nc,dim = intfc->dim;
	double cnp_area,str_length,pz,pv;
	double zcom,vcom;
	double payload = af_params->payload;
	double *g = af_params->gravity;


	char fname[256];

    static FILE *eskfile,*espfile,*egpfile,*efile,*exkfile,*enkfile;
	static FILE *afile,*sfile,*pfile,*vfile;
	static FILE *xcom_file,*vcom_file;
	
    static boolean first = YES;
	if (first)
	{
	    first = NO;
	    
        sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"w");
	    sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"w");
	    sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"w");
	    sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"w");
	    sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"w");
	    sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"w");
	    sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"w");
	    sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"w");
	    sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"w");
	    sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"w");
	    sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"w");
	    sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"w");
        
        fprintf(eskfile,"\"Spr-kinetic energy vs. Time\"\n");
        fprintf(espfile,"\"Spr-potentl energy vs. Time\"\n");
        fprintf(exkfile,"\"Ext-kinetic energy vs. Time\"\n");
        fprintf(egpfile,"\"Ext-potentl energy vs. Time\"\n");
        fprintf(enkfile,"\"Kinetic energy vs. Time\"\n");
        fprintf(efile,"\"Total energy vs. Time\"\n");
        fprintf(afile,"\"Canopy area vs. Time\"\n");
        fprintf(sfile,"\"String length vs. Time\"\n");
        fprintf(pfile,"\"Payload height vs. Time\"\n");
        fprintf(vfile,"\"Payload velo vs. Time\"\n");
        fprintf(xcom_file,"\"COM vs. Time\"\n");
        fprintf(vcom_file,"\"V-COM vs. Time\"\n");
    }
    else
    {
        sprintf(fname,"%s/esk.xg",out_name); eskfile = fopen(fname,"a");
	    sprintf(fname,"%s/esp.xg",out_name); espfile = fopen(fname,"a");
	    sprintf(fname,"%s/egp.xg",out_name); egpfile = fopen(fname,"a");
	    sprintf(fname,"%s/exk.xg",out_name); exkfile = fopen(fname,"a");
	    sprintf(fname,"%s/enk.xg",out_name); enkfile = fopen(fname,"a");
	    sprintf(fname,"%s/eng.xg",out_name); efile = fopen(fname,"a");
	    sprintf(fname,"%s/area.xg",out_name); afile = fopen(fname,"a");
	    sprintf(fname,"%s/str_length.xg",out_name); sfile = fopen(fname,"a");
	    sprintf(fname,"%s/payload.xg",out_name); pfile = fopen(fname,"a");
	    sprintf(fname,"%s/loadvel.xg",out_name); vfile = fopen(fname,"a");
	    sprintf(fname,"%s/xcom.xg",out_name); xcom_file = fopen(fname,"a");
	    sprintf(fname,"%s/vcom.xg",out_name); vcom_file = fopen(fname,"a");
    }
	

    ks = af_params->ks;
    m_s = af_params->m_s;

	esk = esp = epi = epb = egp = exk = enk = 0.0;
	cnp_area = 0.0;
	
	STATE *st;
    for (s = intfc->surfaces; s && *s; ++s)
	{
	    if (wave_type(*s) != ELASTIC_BOUNDARY &&
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) continue;

	    surf = *s;
	    zcom = center_of_mass(Hyper_surf(surf))[2];
	    vcom = center_of_mass_velo(Hyper_surf(surf))[2];
	    
        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
	    {
            cnp_area += tri_area(tri);
            for (j = 0; j < 3; ++j)
            {
                side_length = separation(Point_of_tri(tri)[j],
                                    Point_of_tri(tri)[(j+1)%3],3);
                x_sqr = 0.0;
                for (k = 0; k < 3; ++k)
                {
                    vect[k] = Coords(Point_of_tri(tri)[(j+1)%3])[k]
                        - Coords(Point_of_tri(tri)[j])[k] - tri->side_length0[j]*tri->side_dir0[j][k];
                    
                    x_sqr += sqr(vect[k]);
                }

                if (!is_side_bdry(tri,j))
                {
                    epi += 0.5*ks*x_sqr;
                }
            }
	    }

	    unsort_surf_point(surf);
	    for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf); tri = tri->next)
        {
            for (j = 0; j < 3; ++j)
            {
                p = Point_of_tri(tri)[j];
                if (sorted(p) || Boundary_point(p)) continue;
                for (k = 0; k < dim; ++k)
                {
                    esk += 0.5*m_s*sqr(p->vel[k]);
                    egp += -g[k]*m_s*Coords(p)[k];//TODO: domain doesn't always start at z = 0
                    st = (STATE*)left_state(p);
                    exk += 0.5*m_s*sqr(st->impulse[k]);
                    enk += 0.5*m_s*sqr(p->vel[k] + st->impulse[k]);
                }
                sorted(p) = YES;
            }
        }
	}
	
    epi *= 0.5;	//Each side is counted twice
	
    for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == STRING_HSBDRY)
	    {
	    	kl = af_params->kl;
        	m_l = af_params->m_l;
	    }
        else if (hsbdry_type(*c) == DISKGAP_STRING_HSBDRY)
	    {
	    	kl = af_params->kl_band;
        	m_l = af_params->m_l;
        }
	    else if (hsbdry_type(*c) == MONO_COMP_HSBDRY)
	    {
		    kl = af_params->ks;
        	m_l = af_params->m_s;

            b = (*c)->first;
            BOND_TRI **btris = Btris(b);
            HYPER_SURF* hs = Hyper_surf((*btris)->surface);
            if (wave_type(hs) == ELASTIC_BAND_BOUNDARY)
                kl = af_params->ks_band;

	    }
	    else if (hsbdry_type(*c) == GORE_HSBDRY)
	    {
		    kl = af_params->kg;
        	m_l = af_params->m_g;
	    }
	    else
        {
            continue;
        }
    
        curve = *c;
	    for (b = curve->first; b != NULL; b = b->next)
	    {
            x_sqr = 0.0;
            for (k = 0; k < 3; ++k)
            {
                vect[k] = Coords(b->end)[k] - Coords(b->start)[k] - bond_length0(b)*b->dir0[k];
                x_sqr += sqr(vect[k]);
            }
    
            epb += 0.5*kl*x_sqr;
            if (b != curve->last)
            {
                for (k = 0; k < dim; ++k)
                {
                    esk += 0.5*m_l*sqr(b->end->vel[k]);
                    egp += -g[k]*m_l*Coords(b->end)[k];
                    st = (STATE*)left_state(b->end);
                    exk += 0.5*m_l*sqr(st->impulse[k]);
                    enk += 0.5*m_l*sqr(b->end->vel[k] + st->impulse[k]);
                }
            }
	    }
	}

	esp = epi + epb;

	for (n = intfc->nodes; n && *n; ++n)
	{
	    node = *n;
	    if (is_bdry_node(node)) continue;
	    
        if (is_load_node(node))
	    {
	    	STATE *sl;
	    	pz = Coords(node->posn)[2];
	    	sl = (STATE*)left_state(node->posn);
	    	pv = sl->vel[2];
	    	for (k = 0; k < dim; ++k)
	    	{
		    egp += -g[k]*payload*Coords(node->posn)[k];
		    st = (STATE*)left_state(node->posn);
		    exk += 0.5*payload*sqr(st->impulse[k]);
		    enk += 0.5*payload*sqr(node->posn->vel[k] + st->impulse[k]);
	    	}
	    }
	    else
        {
            if (is_gore_node(node))
            {
                m_l = af_params->m_g;
            }
            else
            {
                m_l = af_params->m_s;
            }

            for (k = 0; k < dim; ++k)
            {
                esk += 0.5*m_l*sqr(node->posn->vel[k]);
                egp += -g[k]*m_l*Coords(node->posn)[k];
        
                st = (STATE*)left_state(node->posn);
                exk += 0.5*m_l*sqr(st->impulse[k]);
                enk += 0.5*m_l*sqr(node->posn->vel[k] + st->impulse[k]);
            }
        }
	}

	nc = 0;
    str_length = 0.0;

    for (c = intfc->curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) != STRING_HSBDRY &&
            hsbdry_type(*c) != DISKGAP_STRING_HSBDRY) continue;
	    str_length += curve_length(*c);
	    nc++;
	}

	if (nc != 0)
    {
	    str_length /= (double)nc;
    }


	fprintf(eskfile,"%16.12f  %16.12f\n",front->time,esk); fclose(eskfile);
    fprintf(espfile,"%16.12f  %16.12f\n",front->time,esp); fclose(espfile);
    fprintf(egpfile,"%16.12f  %16.12f\n",front->time,egp); fclose(egpfile);
    fprintf(exkfile,"%16.12f  %16.12f\n",front->time,exk); fclose(exkfile);
    fprintf(enkfile,"%16.12f  %16.12f\n",front->time,enk); fclose(enkfile);
    fprintf(efile,"%16.12f  %16.12f\n",front->time,esp + egp + enk); fclose(efile);

    fprintf(afile,"%16.12f  %16.12f\n",front->time,cnp_area); fclose(afile);
    fprintf(sfile,"%16.12f  %16.12f\n",front->time,str_length); fclose(sfile);
    fprintf(pfile,"%16.12f  %16.12f\n",front->time,pz); fclose(pfile);
    fprintf(vfile,"%16.12f  %16.12f\n",front->time,pv); fclose(vfile);

    fprintf(xcom_file,"%16.12f  %16.12f\n",front->time,zcom); fclose(xcom_file);
    fprintf(vcom_file,"%16.12f  %16.12f\n",front->time,vcom); fclose(vcom_file);
}	/* end print_airfoil_stat3d_2 */

static void print_strings(
	Front *front,
	char *out_name)
{
	static int dim = FT_Dimension();
	if (dim != 3) return;
	
    CURVE **curve;
	INTERFACE *intfc = front->interf;
        
    boolean status = NO;
    intfc_curve_loop(intfc, curve)
	{
        if (hsbdry_type(*curve) == STRING_HSBDRY ||
            hsbdry_type(*curve) == DISKGAP_STRING_HSBDRY)
        {
            status = YES;
        }
    }
    if (!status) return;

	char dirname[512];
	sprintf(dirname,"%s/string_info", OutName(front));
	
	static boolean first = YES;
    if (first)
	{
	    if (!create_directory(dirname, NO))
	    {
            screen("Cannot create directory %s\n",dirname);
            clean_up(ERROR);
	    }
	}

    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	static double kl = af_params->kl;
	
    char fname[512];
	FILE *length_file, *tension_file;
	
    intfc_curve_loop(intfc, curve)
	{
	    if (hsbdry_type(*curve) != STRING_HSBDRY &&
            hsbdry_type(*curve) != DISKGAP_STRING_HSBDRY) continue;
	    
        /* in the breaking strings case, one curve breaks into two. */
	    boolean new_curve = NO;
	    if (af_params->string_hash.count(Gindex(*curve)) == 0)
	    {
		    int value = af_params->string_hash.size();
		    af_params->string_hash[Gindex(*curve)] = value;
		    new_curve = YES;
	    }
	    
        sprintf(fname,"%s/str_length-%02d.xg",dirname,af_params->string_hash[Gindex(*curve)]);
	    if (first || new_curve)
        {
		    length_file = fopen(fname, "w");
        }
        else
        {
            length_file = fopen(fname, "a");
        }
        
        sprintf(fname,"%s/str_tension-%02d.xg",dirname, af_params->string_hash[Gindex(*curve)]);

	    if (first || new_curve)
        {
            tension_file = fopen(fname, "w");
        }
        else
        {
            tension_file = fopen(fname, "a");
        }
        
        BOND *bond;
	    double len0 = 0.0, len = 0.0;
	    curve_bond_loop(*curve, bond)
	    {
		    len0 += bond->length0;
		    len += distance_between_positions(Coords(bond->start),Coords(bond->end),dim);
	    }

	    double ten = (len - len0) * kl / (I_NumOfCurvePoints(*curve)-1);
	    
        fprintf(length_file,"%16.12f  %16.12f\n",front->time,len); fclose(length_file);
	    fprintf(tension_file,"%16.12f  %16.12f\n",front->time,ten); fclose(tension_file);
	}

    if (first)
    {
        first = NO;
    }
}	/* end print_strings */

//TODO: Rewrite similiar to strain limiting functions
static void record_stretching_length(
	SURFACE *surf,
	char *out_name,
	double time)
{
	int dir;
	CURVE **c,*c1,*c2;
	static FILE *lfile;
	char lname[512];
	C_PARAMS *c_params;
	POINT *p1,*p2;
	double length;

	c1 = c2 = NULL;
	for (c = surf->pos_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY &&
		(*c)->extra != NULL)
	    {
		c1 = *c;
		c_params = (C_PARAMS*)c1->extra;
	    }
	    else if (hsbdry_type(*c) == FIXED_HSBDRY)
		c2 = *c;
	}
	for (c = surf->neg_curves; c && *c; ++c)
	{
	    if (hsbdry_type(*c) == MONO_COMP_HSBDRY &&
		(*c)->extra != NULL)
	    {
		c1 = *c;
		c_params = (C_PARAMS*)c1->extra;
	    }
	    else if (hsbdry_type(*c) == FIXED_HSBDRY)
		c2 = *c;
	}
	if (c1 == NULL || c2 == NULL) return;
	if (lfile == NULL)
	{
	    sprintf(lname,"%s/length.xg",out_name);
	    lfile = fopen(lname,"w");
	    fprintf(lfile,"\"Length vs. Time\"\n");
	}
	dir = c_params->dir;
	p1 = c1->first->end;
	p2 = c2->first->end;
	length = fabs(Coords(p1)[dir] - Coords(p2)[dir]);
	fprintf(lfile,"%f  %f\n",time,length);
	fflush(lfile);
}	/* end record_stretching_length */

static void print_rgb3d(
        Front *front,
        char *out_name)
{
    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    INTERFACE *intfc = front->interf;
    SURFACE **s;
    FILE *rgxfile, *rgvfile;
    char fname[200];

    static boolean first = YES;
    intfc_surface_loop(intfc, s)
    {
        if (wave_type(*s) != MOVABLE_BODY_BOUNDARY) continue;
        
        sprintf(fname,"%s/rg_com%02d",OutName(front),body_index(*s));
        if (first)
        {
            rgxfile = fopen(fname,"w");
        }
        else
        {
            rgxfile = fopen(fname,"a");
        }

        sprintf(fname,"%s/rg_vel%02d",OutName(front),body_index(*s));
        if (first)
        {
            rgvfile = fopen(fname,"w");
        }
        else
        {
            rgvfile = fopen(fname,"a");
        }

        fprintf(rgxfile, "%16.12f  %16.12f \n", front->time,center_of_mass(*s)[2]);
        fclose(rgxfile);

        fprintf(rgvfile, "%16.12f  %16.12f \n", front->time,center_of_mass_velo(*s)[2]);
        fclose(rgvfile);
    }
        
    if (first)
    {
        first = NO;
    }
}       /* end print_rgb3d */

//TODO: REVIEW FOR ACCURACY 
static void print_drag3d(
        Front *front,
        char *out_name)
{
    INTERFACE *intfc = front->interf;
    SURFACE **s;
    TRI *tri;
    POINT *point;
    Locstate sl, sr;
    int i,j,count = 0;
    double pres_drop, area;
    char fname[200];
    FILE *dfile, *xlfile, *ylfile;
    static boolean first = YES;
    
    double (*getStateVel[3])(POINTER) = {
        getStateXvel, getStateYvel, getStateZvel
    };
    
    FILE* pafile;
    FILE *xforce, *yforce, *zforce;

    int dim = FT_Dimension();
    if (dim != 3) return;

    AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
    if (af_params->no_fluid == YES) return;
    
    /*find the freestream velocity*/
    /*freestream is the air far upstream of an aerodynamic body*/
    double free_vel[MAXD] = {0};
    double free_vel_dir[MAXD] = {0};
    double drag[MAXD]={0}, lift[MAXD]={0};
    double fvel_mag = 0.0;

    //TODO: Record the inlet/freestream velocity on intitialization instead of doing this...
    //
    //looks for the constant state dirichlet boundary (INLET)
    for (s = intfc->surfaces; s && *s; ++s)
    {
        HYPER_SURF *hs = Hyper_surf(*s);
        if (wave_type(hs) == DIRICHLET_BOUNDARY &&
            is_bdry(Surface_of_hs(hs)) && boundary_state(hs))
        {
            for (i = 0; i < dim; i++)
                free_vel[i] = getStateVel[i](boundary_state(hs));
            fvel_mag = Mag3d(free_vel);
            if (fabs(fvel_mag) > MACH_EPS) break;
        }
    }

    //printf("free_vel = %g %g %g\n",free_vel[0],free_vel[1],free_vel[2]);

    /*normalize freestream vel*/
    if (fabs(fvel_mag) > MACH_EPS)
    {
        for (i = 0; i < dim; i++)
            free_vel_dir[i] = free_vel[i]/fvel_mag;
    }
    else
    {
        /*default freestream direction is upward*/
        for (i = 0; i < dim; i++)
            free_vel_dir[i] = 0.0;
        free_vel_dir[dim-1] = 1.0;
    }

    /*compute total force on fabric canopy*/
    /*for multiparachute problem*/
    /* one file is generated for each canopy*/
    int fcount = 0;
    for (s = intfc->surfaces; s && *s; ++s)
    {
        if ((wave_type(*s) != ELASTIC_BOUNDARY && 
            wave_type(*s) != ELASTIC_BAND_BOUNDARY) || is_bdry(*s)) continue;

        /* drag force */
        sprintf(fname,"%s/drag-%d.xg",OutName(front),fcount);
        if (first)
            dfile = fopen(fname,"w");
        else
            dfile = fopen(fname,"a");
        
        /* lift force */
        sprintf(fname,"%s/xlift-%d.xg",OutName(front),fcount);
        if (first)
            xlfile = fopen(fname,"w");
        else
            xlfile = fopen(fname,"a");
        
        sprintf(fname,"%s/ylift-%d.xg",OutName(front),fcount);
        if (first)
            ylfile = fopen(fname,"w");
        else
            ylfile = fopen(fname,"a");

        /* projected area */
        sprintf(fname,"%s/parea-%d.xg",OutName(front),fcount);
        if (first)
            pafile = fopen(fname,"w");
        else
            pafile = fopen(fname,"a");

        /* total force in x direction */
        sprintf(fname,"%s/force-x-%d.xg",OutName(front),fcount);
        if (first)
            xforce = fopen(fname,"w");
        else
            xforce = fopen(fname,"a");
        
        /* total force in y direction */
        sprintf(fname,"%s/force-y-%d.xg",OutName(front),fcount);
        if (first)
            yforce = fopen(fname,"w");
        else
            yforce = fopen(fname,"a");
        
        /* total force in z direction */
        sprintf(fname,"%s/force-z-%d.xg",OutName(front),fcount);
        if (first)
            zforce = fopen(fname,"w");
        else
            zforce = fopen(fname,"a");
        
        fcount++;
        
        double parea = 0.0;
        double force[MAXD] = {0.0};
        
        for (tri = first_tri(*s); !at_end_of_tri_list(tri,*s); tri = tri->next)
        {
            pres_drop = 0.0;
            double pl = 0.0;
            double pr = 0.0;
            double vel_tri[MAXD] = {0.0};
            double centroid[MAXD] = {0.0};

            for(i = 0; i < 3; i++)
            {
                point = Point_of_tri(tri)[i];
            
                sl = (STATE*)left_state(point);
                sr = (STATE*)right_state(point);
            
                pl += getStatePres(sl);
                pr += getStatePres(sr);
            
                pres_drop += getStatePres(sl) - getStatePres(sr);
            
                for (int k = 0; k < 3; ++k)
                {
                    centroid[k] += Coords(point)[k];
                    vel_tri[k] += getStateVel[k](sl);
                }
            }

            pl /= 3.0;
            pr /= 3.0;
            pres_drop /= 3.0;
            
            for (int k = 0; k < 3; ++k)
            {
                centroid[k] /= 3.0;
                vel_tri[k] /= 3.0;
            }
            
            double unit_nor_tri[MAXD];
            auto nor_tri = Tri_normal(tri);//NOTE: Tri_normal() not unit length
            double mag_nor = Mag3d(nor_tri);
            double area_tri = tri_area(tri);
            
            double force_tri[MAXD];
            for (i = 0; i < dim; i++)
            {
                unit_nor_tri[i] = nor_tri[i]/mag_nor;
                force_tri[i] = pres_drop*unit_nor_tri[i]*area_tri;
                force[i] += force_tri[i];
            }
            parea += Dot3d(unit_nor_tri,free_vel_dir)*area_tri;//projected to xy plane

            //TODO: put inside a if (debugging(...)) block
            /*
            printf("pres_drop = %g - %g = %g  |  unit_nor_tri = %g %g %g  |"
                    "  force_tri = %g %g %g  |  tri_cen %g %g %g\n",
                    pl,pr,pres_drop,unit_nor_tri[0],unit_nor_tri[1],unit_nor_tri[2],
                    force_tri[0],force_tri[1],force_tri[2],
                    centroid[0],centroid[1],centroid[2]);*/
        }
        
        /*compute drag force and lift force*/
        double mag_drag = Dot3d(force,free_vel_dir);
        for (i = 0; i < dim; i++)
        {
            drag[i] = mag_drag*free_vel_dir[i];
            lift[i] = force[i] - drag[i];
        }

        //TODO: put inside a if (debugging(...)) block
        /*
        printf("\t force = %g %g %g  |  drag = %g %g %g  |  lift = %g %g %g\n",
                force[0],force[1],force[2],
                drag[0],drag[1],drag[2],
                lift[0],lift[1],lift[2]);*/

            //fprintf(dfile,"%16.12f  %16.12f\n",front->time,Mag3d(drag));
        fprintf(dfile,"%16.12f  %16.12f\n",front->time,drag[2]);
        fclose(dfile);

            /*fprintf(lfile,"%16.12f  %16.12f\n",front->time,Mag3d(lift));
            fclose(lfile);*/
        fprintf(xlfile,"%16.12f  %16.12f\n",front->time,lift[0]);
        fclose(xlfile);
        fprintf(ylfile,"%16.12f  %16.12f\n",front->time,lift[1]);
        fclose(ylfile);

        fprintf(xforce,"%16.12f  %16.12f\n",front->time,force[0]);
        fclose(xforce);
        fprintf(yforce,"%16.12f  %16.12f\n",front->time,force[1]);
        fclose(yforce);
        fprintf(zforce,"%16.12f  %16.12f\n",front->time,force[2]);
        fclose(zforce);

        fprintf(pafile,"%16.12f  %16.12f\n",front->time,parea);
        fclose(pafile);

    }

    first = NO;
    return;
}       /* end print_drag3d */

//SAVING THESE JUST FOR REFERENCE RIGHT NOW
static void print_airfoil_stat2d(
	Front *front,
	char *out_name)
{
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	switch (af_params->spring_model)
	{
	case MODEL1:
	    print_airfoil_stat2d_1(front,out_name);
	    break;
	case MODEL2:
	    print_airfoil_stat2d_2(front,out_name);
	    break;
	case MODEL3:
	default:
	    (void) printf("print_airfoil_stat2d_12() not implemented!\n");
	}
}	/* end print_airfoil_stat2d */

static void print_airfoil_stat2d_1(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
	CURVE **c,*curve;
	BOND *b;
	POINT *p;
	static FILE *ekfile,*epfile,*exfile,*efile,*sfile;
	char fname[256];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double ek,ep,enp;
	double kl,m_l,x_diff;
	int i,dim = intfc->dim;
	double str_length;
	STRING_NODE_TYPE start_type = af_params->start_type;
        STRING_NODE_TYPE end_type = af_params->end_type;
	double *g,payload;

	if (ekfile == NULL && pp_mynode() == 0)
        {
	    sprintf(fname,"%s/ek.xg",out_name);
            ekfile = fopen(fname,"w");
	    sprintf(fname,"%s/ep.xg",out_name);
            epfile = fopen(fname,"w");
	    sprintf(fname,"%s/ex.xg",out_name);
            exfile = fopen(fname,"w");
	    sprintf(fname,"%s/en.xg",out_name);
            efile = fopen(fname,"w");
	    sprintf(fname,"%s/str_length.xg",out_name);
            sfile = fopen(fname,"w");
            fprintf(ekfile,"\"Kinetic enegy vs. time\"\n");
            fprintf(epfile,"\"Potential enegy vs. time\"\n");
            fprintf(exfile,"\"External enegy vs. time\"\n");
            fprintf(efile,"\"Total enegy vs. time\"\n");
            fprintf(sfile,"\"String length vs. time\"\n");
        }

	kl = af_params->kl;
        m_l = af_params->m_l;
	payload = af_params->payload;
	g = af_params->gravity;
	ek = ep = enp = str_length = 0.0;
	for (c = intfc->curves; c && *c; ++c)
	{
		if (wave_type(*c) != ELASTIC_BOUNDARY &&
		    wave_type(*c) != ELASTIC_STRING)
		    continue;

		curve = *c;
		p = curve->first->start;
		if (start_type == FREE_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*m_l*sqr(p->vel[i]);
		    }
		}
		else if(start_type == LOADED_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*payload*sqr(p->vel[i]);
			enp -= payload*g[i]*Coords(p)[i];
		    }
		}
		p = curve->last->end;
		if (end_type == FREE_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*m_l*sqr(p->vel[i]);
		    }
		}
		else if(end_type == LOADED_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*payload*sqr(p->vel[i]);
			enp -= payload*g[i]*Coords(p)[i];
		    }
		}

		for (b = curve->first; b != NULL; b = b->next)
		{
		    p = b->end;
		    x_diff = bond_length(b) - bond_length0(b);
		    if (b != curve->last)
		    {
		    	for (i = 0; i < dim; ++i)
			{
		    	    ek += 0.5*m_l*sqr(p->vel[i]);
			}
		    }
		    ep += 0.5*kl*sqr(x_diff);
		    str_length += bond_length(b);
		}
	}
	if (pp_mynode() == 0)
	{
	    fprintf(ekfile,"%16.12f  %16.12f\n",front->time,ek);
            fprintf(epfile,"%16.12f  %16.12f\n",front->time,ep);
            fprintf(exfile,"%16.12f  %16.12f\n",front->time,enp);
            fprintf(efile,"%16.12f  %16.12f\n",front->time,ek+ep+enp);
            fprintf(sfile,"%16.12f  %16.12f\n",front->time,str_length);
	    fflush(ekfile);
	    fflush(epfile);
	    fflush(exfile);
	    fflush(efile);
	    fflush(sfile);
	}
}	/* end print_airfoil_stat2d_1 */

static void print_airfoil_stat2d_2(
	Front *front,
	char *out_name)
{
	INTERFACE *intfc = front->interf;
	CURVE **c,*curve;
	BOND *b;
	POINT *p;
	static FILE *ekfile,*epfile,*exfile,*efile,*sfile;
	char fname[256];
	AF_PARAMS *af_params = (AF_PARAMS*)front->extra2;
	double ek,ep,enp;
	double kl,m_l,x_diff;
	int i,dim = intfc->dim;
	double str_length;
	STRING_NODE_TYPE start_type = af_params->start_type;
        STRING_NODE_TYPE end_type = af_params->end_type;
	double *g,payload;
	double vect[MAXD],len0;

	if (ekfile == NULL && pp_mynode() == 0)
        {
	    sprintf(fname,"%s/ek.xg",out_name);
            ekfile = fopen(fname,"w");
	    sprintf(fname,"%s/ep.xg",out_name);
            epfile = fopen(fname,"w");
	    sprintf(fname,"%s/ex.xg",out_name);
            exfile = fopen(fname,"w");
	    sprintf(fname,"%s/en.xg",out_name);
            efile = fopen(fname,"w");
	    sprintf(fname,"%s/str_length.xg",out_name);
            sfile = fopen(fname,"w");
            fprintf(ekfile,"\"Kinetic enegy vs. time\"\n");
            fprintf(epfile,"\"Potential enegy vs. time\"\n");
            fprintf(exfile,"\"External enegy vs. time\"\n");
            fprintf(efile,"\"Total enegy vs. time\"\n");
            fprintf(sfile,"\"String length vs. time\"\n");
        }

	kl = af_params->kl;
        m_l = af_params->m_l;
	payload = af_params->payload;
	g = af_params->gravity;
	ek = ep = enp = str_length = 0.0;
	for (c = intfc->curves; c && *c; ++c)
	{
		if (wave_type(*c) != ELASTIC_BOUNDARY &&
		    wave_type(*c) != ELASTIC_STRING)
		    continue;

		curve = *c;
		p = curve->first->start;
		if (start_type == FREE_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*m_l*sqr(p->vel[i]);
		    }
		}
		else if(start_type == LOADED_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*payload*sqr(p->vel[i]);
			enp -= payload*g[i]*Coords(p)[i];
		    }
		}
		p = curve->last->end;
		if (end_type == FREE_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*m_l*sqr(p->vel[i]);
		    }
		}
		else if(end_type == LOADED_END)
		{
		    for (i = 0; i < dim; ++i)
		    {
		    	ek += 0.5*payload*sqr(p->vel[i]);
			enp -= payload*g[i]*Coords(p)[i];
		    }
		}

		for (b = curve->first; b != NULL; b = b->next)
		{
		    p = b->end;
		    x_diff = bond_length(b) - bond_length0(b);
		    if (b != curve->last)
		    {
		    	for (i = 0; i < dim; ++i)
			{
		    	    ek += 0.5*m_l*sqr(p->vel[i]);
			}
		    }
		    len0 = bond_length0(b);
		    x_diff = 0.0;
		    for (i = 0; i < dim; ++i)
		    {
			vect[i] = Coords(b->end)[i] - Coords(b->start)[i]
					- len0*b->dir0[i];
			x_diff += sqr(vect[i]);
		    }
		    ep += 0.5*kl*x_diff;
		    str_length += bond_length(b);
		}
	}
	if (pp_mynode() == 0)
	{
	    fprintf(ekfile,"%16.12f  %16.12f\n",front->time,ek);
            fprintf(epfile,"%16.12f  %16.12f\n",front->time,ep);
            fprintf(exfile,"%16.12f  %16.12f\n",front->time,enp);
            fprintf(efile,"%16.12f  %16.12f\n",front->time,ek+ep+enp);
            fprintf(sfile,"%16.12f  %16.12f\n",front->time,str_length);
	    fflush(ekfile);
	    fflush(epfile);
	    fflush(exfile);
	    fflush(efile);
	    fflush(sfile);
	}
}	/* end print_airfoil_stat2d_2 */
