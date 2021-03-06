#include "cgal_intfc.h"
#include "rigidbody.h"


static void prompt_for_rigid_body_params(FILE*,RG_PARAMS*); 
//static void prompt_for_rigid_body_params(int,FILE*,RG_PARAMS*); 
static void set_rgbody_params(RG_PARAMS*,HYPER_SURF*);
static void prompt_for_velocity_func(int,FILE*,RG_PARAMS*);
//static void prompt_for_velocity_func(int,char*,RG_PARAMS*);
//static void sine_vel_func(Front*,POINTER,double*,double*);

static void initSingleRigidBody(FILE*,Front*);
static void initMultiRigidBodies(FILE*,Front*,int);

static void init_rigid_sphere(FILE*,Front*);
static void init_rigid_box(FILE*,Front*);
static void init_rigid_cylinder(FILE*,Front*);
static void init_rigid_human(FILE*,Front*);

// functions for changing the human body surface
static void surf_com_translation(SURFACE*,double*);
static void surf_enlargement(SURFACE*,double);


extern void initRigidBodies(Front* front)
{
	FILE *infile = fopen(InName(front),"r");
    
    char string[100];
    if (CursorAfterStringOpt(infile,"Enter yes to add rigid body:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	}

    if (string[0] != 'y' && string[0] != 'Y')
    {
        fclose(infile);
	    return;
    }

	int num_rgb = 1;
	if (CursorAfterStringOpt(infile,"Enter the number of rigid bodies:"))
	{
	    fscanf(infile,"%d",&num_rgb);
	    (void) printf("%d\n",num_rgb);
	}

	if (num_rgb == 1)
	    initSingleRigidBody(infile,front);
	else if (num_rgb > 1)
	    initMultiRigidBodies(infile,front,num_rgb);


	RG_PARAMS* rgb_params = (RG_PARAMS*)front->extra3;

    CURVE **c;
    SURFACE **s;

    for (s = front->interf->surfaces; s && *s; ++s)
    {
        if (wave_type(*s) == MOVABLE_BODY_BOUNDARY ||
            wave_type(*s) == NEUMANN_BOUNDARY)
        {
            if (wave_type(*s) == NEUMANN_BOUNDARY)
                rgb_params->is_fixed = true;
            prompt_for_rigid_body_params(infile,rgb_params);
            set_rgbody_params(rgb_params,Hyper_surf(*s));
        }
    }
    
    fclose(infile);
}

/*
extern void rgb_init(Front* front, RG_PARAMS* rgb_params)
{
        int dim = FT_Dimension();
        if (dim == 1) return;

        CURVE **c;
        SURFACE **s;
        char* inname = InName(front);

        if (dim == 2)
        {
            for (c = front->interf->curves; c && *c; ++c)
            {
                if (wave_type(*c) == MOVABLE_BODY_BOUNDARY)
                {
                    prompt_for_rigid_body_params(dim,inname,rgb_params);
                    set_rgbody_params(rgb_params,Hyper_surf(*c));
                }
            }
        }
        else
        {
            for (s = front->interf->surfaces; s && *s; ++s)
            {
                if (wave_type(*s) == MOVABLE_BODY_BOUNDARY)
                {
                    prompt_for_rigid_body_params(dim,inname,rgb_params);
                    set_rgbody_params(rgb_params,Hyper_surf(*s));
                }
            }
        }

}*/       /* end rgb_init */


static void initSingleRigidBody(
	FILE *infile,
	Front *front)
{
        char string[100];
	(void) printf("Available type of rigid body include:\n");
	(void) printf("\tSphere (S)\n");
	(void) printf("\tBox (B)\n");
	(void) printf("\tHuman (H)\n");
	(void) printf("\tCylinder (C)\n");
	CursorAfterString(infile,"Enter type of rigid body:");
	fscanf(infile,"%s",string);
	(void) printf("%s\n",string);
	switch (string[0])
	{
	case 's':
	case 'S':
	    init_rigid_sphere(infile, front);
            break;
    case 'b':
    case 'B':
        init_rigid_box(infile, front);
        break;
	case 'c':
	case 'C':
	    init_rigid_cylinder(infile, front);
	    break;
	case 'h':
	case 'H':
	    init_rigid_human(infile, front);
	    break;
    default:
        (void) printf("Unknow type of rigid body!\n");
        clean_up(ERROR);
        }
}	/* end initSingleRigidBody */

static void initMultiRigidBodies(
	FILE *infile,
	Front *front,
	int num_rgb)
{
	int i;
	char string[100];
	INTERFACE *cur_intfc = current_interface();

	for (i = 0; i < num_rgb; ++i)
	{
	    sprintf(string,"For rigid body %d", i+1);
	    CursorAfterString(infile, string);
	    printf("\n");
	    initSingleRigidBody(infile,front);
	}
	set_current_interface(cur_intfc);
}	/* initMultiRigidBodies */

static void init_rigid_sphere(
	FILE *infile,
	Front *front)
{
        char string[100];
        double cen[MAXD];
        double radius,radii[MAXD];
        int w_type;
        int i,dim = FT_Dimension();
        int neg_comp,pos_comp;
        SURFACE *surf;

        CursorAfterString(infile,"Enter center of the sphere:");
        fscanf(infile,"%lf %lf %lf",cen,cen+1,cen+2);
        (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
        CursorAfterString(infile,"Enter radius of the sphere:");
        fscanf(infile,"%lf",&radius);
        (void) printf("%f\n",radius);
        for (i = 0; i < dim; ++i) radii[i] = radius;

        (void) printf("Rigid body can be fixed (F) or Movable (M)\n");
        (void) printf("The default is Movable (M)\n");
        w_type = MOVABLE_BODY_BOUNDARY;
        neg_comp = SOLID_COMP;
        pos_comp = GAS_COMP2;
        if (CursorAfterStringOpt(infile,"Type yes if the rigid body is fixed:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                w_type = NEUMANN_BOUNDARY;
        }
        bool cgal_mesh = false;
        if (CursorAfterStringOpt(infile,"Type yes to use CGAL for rigid body:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                cgal_mesh = true;
        }

        if (cgal_mesh)
        {
            CGAL_MakeEllipsoidalSurf(front,cen,radii,neg_comp,pos_comp,w_type,
                                        1,&surf);
        }
        else
        {
            FT_MakeEllipticSurf(front,cen,radii,neg_comp,pos_comp,w_type,2,
                                        &surf);
        }
}	/* end init_rigid_sphere */

static void init_rigid_box(
	FILE *infile,
	Front *front)
{
	char string[100];
	double cen[MAXD], edge[MAXD];
	int w_type;
	int neg_comp,pos_comp;
	SURFACE *surf;

        CursorAfterString(infile,"Enter center of the box:");
        fscanf(infile,"%lf %lf %lf",cen,cen+1,cen+2);
        (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
        CursorAfterString(infile,"Enter edges of the box:");
        fscanf(infile,"%lf %lf %lf",edge,edge+1,edge+2);
        (void) printf("%f %f %f\n",edge[0],edge[1],edge[2]);
        (void) printf("Rigid body can be fixed (F) or Movable (M)\n");
        (void) printf("The default is Movable (M)\n");
        w_type = MOVABLE_BODY_BOUNDARY;
        neg_comp = SOLID_COMP;
        pos_comp = LIQUID_COMP2;
        if (CursorAfterStringOpt(infile,"Type yes if the rigid body is fixed:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                w_type = NEUMANN_BOUNDARY;
        }


        bool cgal_mesh = false;
        if (CursorAfterStringOpt(infile,"Type yes to use CGAL for rigid body:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                cgal_mesh = true;
        }

        if (cgal_mesh)
        {
            CGAL_MakeCuboidSurf(front,cen,edge,neg_comp,pos_comp,w_type,1,
                                    &surf);
        }
        else
        {
            FT_MakeCuboidSurf(front,cen,edge,neg_comp,pos_comp,w_type,2,&surf);
        }
}	/* end init_rigid_box */

static void init_rigid_cylinder(
	FILE *infile,
	Front *front)
{
	char string[100];
	double cen[MAXD], radius, height;
    int idir;
	SURFACE* surf;
	int w_type;
	int neg_comp,pos_comp;

        CursorAfterString(infile,"Enter center of the cylinder:");
            fscanf(infile,"%lf %lf %lf",cen,cen+1,cen+2);
            (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
        CursorAfterString(infile,"Enter radius of the cylinder:");
            fscanf(infile,"%lf",&radius);
            (void) printf("%f\n",radius);
        CursorAfterString(infile,"Enter height of the cylinder:");
            fscanf(infile,"%lf",&height);
            (void) printf("%f\n",height);
        CursorAfterString(infile,"Enter central axis of the cylinder:");
            fscanf(infile,"%ld",&idir);
            (void) printf("%d\n",idir);
        (void) printf("Rigid body can be fixed (F) or Movable (M)\n");
        (void) printf("The default is Movable (M)\n");
        w_type = MOVABLE_BODY_BOUNDARY;
            neg_comp = SOLID_COMP;
            pos_comp = LIQUID_COMP2;
        if (CursorAfterStringOpt(infile,
                    "Type yes if the rigid body is fixed:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                w_type = NEUMANN_BOUNDARY;
        }

        bool cgal_mesh = false;
        if (CursorAfterStringOpt(infile,"Type yes to use CGAL for rigid body:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                cgal_mesh = true;
        }

        if (cgal_mesh)
        {
            CGAL_MakeCylindricalSurf(front,cen,radius,height,
                    idir,neg_comp,pos_comp,w_type,1,&surf);
            
        }
        else
        {
	        FT_MakeCylinderSurf(front,cen,radius,height/2,
                    idir,neg_comp,pos_comp,w_type,&surf);
        }
}	/* end init_rigid_cylinder */

static void init_rigid_human(
	FILE *infile,
	Front *front)
{
	char vtk_name[100], string[100];
	double cen[MAXD], coeff;
	int w_type;
	int neg_comp,pos_comp;
	SURFACE *surf;

	sprintf(string,"Enter the vtk file name for human body:");
	CursorAfterString(infile,string);
	fscanf(infile,"%s",vtk_name);
	(void) printf("%s\n",vtk_name);
	neg_comp = SOLID_COMP;
	pos_comp = LIQUID_COMP2;
	read_vtk_surface(front->interf,neg_comp,pos_comp,vtk_name,&surf);
	if(CursorAfterStringOpt(infile,"Enter center of the human body:"))
	{
	    fscanf(infile,"%lf %lf %lf",cen,cen+1,cen+2);
	    (void) printf("%f %f %f\n",cen[0],cen[1],cen[2]);
	    surf_com_translation(surf, cen);
	}
	if(CursorAfterStringOpt(infile,"Enter enlargement coefficient:"))
	{
	    fscanf(infile,"%lf",&coeff);
	    (void) printf("%f\n",coeff);
	    surf_enlargement(surf, coeff);
	}
	(void) printf("Rigid body can be fixed (F) or Movable (M)\n");
	(void) printf("The default is Movable (M)\n");
	wave_type(surf) = MOVABLE_BODY_BOUNDARY;
	if (CursorAfterStringOpt(infile,
	    		"Type yes if the rigid body is fixed:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] == 'y' || string[0] == 'Y')
	        wave_type(surf) = NEUMANN_BOUNDARY;
	}
}	/* end init_rigid_human */

static void surf_com_translation(
	SURFACE *surf,
	double *cen)
{
	TRI *tri;
	double com[3] = {0.0};
	int count = 0;
	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (int i = 0; i < 3; ++i)
	    {
	        POINT *p = Point_of_tri(tri)[i];
		if (sorted(p)) continue;
	    	for (int j = 0; j < 3; ++j)
	    	    com[j] += Coords(p)[j];
	    	count++;
	        sorted(p) = YES;
	    }
	}
	for (int j = 0; j < 3; ++j)
	    com[j] /= count;
	printf("com = %f %f %f \n", com[0], com[1], com[2]);
	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                POINT *p = Point_of_tri(tri)[i];
                if (sorted(p)) continue;
                for (int j = 0; j < 3; ++j)
                    Coords(p)[j] = Coords(p)[j] - com[j] +  cen[j];
                sorted(p) = YES;
            }
        }
}

static void surf_enlargement(
	SURFACE *surf,
	double coeff)
{
	TRI *tri;
	double com[3] = {0.0};
	int count = 0;
	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
	{
	    for (int i = 0; i < 3; ++i)
	    {
	        POINT *p = Point_of_tri(tri)[i];
		if (sorted(p)) continue;
	    	for (int j = 0; j < 3; ++j)
	    	    com[j] += Coords(p)[j];
	    	count++;
	        sorted(p) = YES;
	    }
	}
	for (int j = 0; j < 3; ++j)
            com[j] /= count;
	unsort_surf_point(surf);
	surf_tri_loop(surf,tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                POINT *p = Point_of_tri(tri)[i];
                if (sorted(p)) continue;
                for (int j = 0; j < 3; ++j)
                    Coords(p)[j] = (Coords(p)[j] - com[j]) * coeff + com[j];
                sorted(p) = YES;
            }
        }
}

static void prompt_for_rigid_body_params(
        FILE* infile,
        RG_PARAMS* rgb_params)
{
        char msg[100],s[100],ss[100];
        static int count = 1;

        int dim = rgb_params->dim;
        boolean is_preset_motion = NO;
        double mag_dir = 0.0;

        if (debugging("rgbody"))
            (void) printf("Enter prompt_for_rigid_body_params()\n");

        if (count == 1)
        {
            rgb_params->no_fluid = NO;
            if (CursorAfterStringOpt(infile,
                "Entering yes to turn off fluid solver: "))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                    rgb_params->no_fluid = YES;
            }
        }

        rgb_params->body_index = count++;
        if (rgb_params->is_fixed) return;

        sprintf(s, "For rigid body %d", rgb_params->body_index);
        CursorAfterString(infile, s); printf("\n");
        long idpos = ftell(infile);

        CursorAfterString(infile,"Type yes if motion is preset: ");
        fscanf(infile,"%s",s);
        (void) printf("%s\n",s);
        if (s[0] == 'y' || s[0] == 'Y')
        {
            (void) printf("Available preset motion types are:\n");
            (void) printf("\tPRESET_TRANSLATION\n");
            (void) printf("\tPRESET_ROTATION\n");
            (void) printf("\tPRESET_MOTION (general)\n");
            CursorAfterString(infile,"Enter type of preset motion: ");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[7])
            {
            case 'M':
                rgb_params->motion_type = PRESET_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = PRESET_COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = PRESET_TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = PRESET_ROTATION;
                break;
            default:
                (void) printf("Unknown type of preset motion!\n");
                clean_up(ERROR);
            }
            (void) fseek(infile,idpos,SEEK_SET);
        }
        else
        {
            (void) printf("Available dynamic motion types are:\n");
            (void) printf("\tFREE_MOTION\n");
            (void) printf("\tCOM_MOTION (center of mass)\n");
            (void) printf("\tTRANSLATION\n");
            (void) printf("\tROTATION\n");
            CursorAfterString(infile,"Enter type of dynamic motion:");
            fscanf(infile,"%s",s);
            (void) printf("%s\n",s);
            switch(s[0])
            {
            case 'F':
                rgb_params->motion_type = FREE_MOTION;
                break;
            case 'C':
                rgb_params->motion_type = COM_MOTION;
                break;
            case 'T':
                rgb_params->motion_type = TRANSLATION;
                break;
            case 'R':
                rgb_params->motion_type = ROTATION;
                break;
            default:
                (void) printf("Unknown type of motion!\n");
                clean_up(ERROR);
            }
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            mag_dir = 0.0;
            CursorAfterString(infile,"Enter the direction of motion:");
            for (int i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->translation_dir[i]);
                (void) printf("%f ",rgb_params->translation_dir[i]);
                mag_dir += sqr(rgb_params->translation_dir[i]);
            }
            (void) printf("\n");
            mag_dir = sqrt(mag_dir);
            for (int i = 0; i < dim; ++i)
                rgb_params->translation_dir[i] /= mag_dir;
            (void) fseek(infile,idpos,SEEK_SET);
        }
    
        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION)
        {
            sprintf(msg,"Enter the total mass for rigid body:");
            CursorAfterString(infile,msg);
            fscanf(infile,"%lf",&rgb_params->total_mass);
            (void) printf("%f\n",rgb_params->total_mass);
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == COM_MOTION ||
            rgb_params->motion_type == TRANSLATION ||
            rgb_params->motion_type == PRESET_MOTION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            sprintf(msg,"Enter the initial center of mass for rigid body:");
            CursorAfterString(infile,msg);
            for (int i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
                (void) printf("%f ",rgb_params->center_of_mass[i]);
            }
            (void) fseek(infile,idpos,SEEK_SET);

            (void) printf("\n");
            sprintf(msg,"Enter the initial center of mass velocity:");
            CursorAfterString(infile,msg);
            for (int i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->cen_of_mass_velo[i]);
                (void) printf("%f ",rgb_params->cen_of_mass_velo[i]);
            }
            (void) printf("\n");
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (rgb_params->motion_type == PRESET_ROTATION)
        {
            /* 2D preset rotation is always about the z-axis */
            /* 3D preset rotation axis along rotation_dir */
            if (dim == 3)
            {
                mag_dir = 0.0;
                CursorAfterString(infile,"Enter the direction of rotation:");
                for (int i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                    (void) printf("%f ",rgb_params->rotation_dir[i]);
                    mag_dir += sqr(rgb_params->rotation_dir[i]);
                }
                (void) printf("\n");
                mag_dir = sqrt(mag_dir);
                for (int i = 0; i < dim; ++i)
                    rgb_params->rotation_dir[i] /= mag_dir;
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (int i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
            (void) fseek(infile,idpos,SEEK_SET);

            /* Center of axis is the coordinate of a point on the axis */
            CursorAfterString(infile,"Enter rotation center:");
            for (int i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                (void) printf("%f ",rgb_params->rotation_cen[i]);
            }
            (void) fseek(infile,idpos,SEEK_SET);

            (void) printf("\n");
            CursorAfterString(infile,"Enter preset angular velocity:");
            fscanf(infile,"%lf",&rgb_params->angular_velo);
            (void) printf("%f\n",rgb_params->angular_velo);
            (void) fseek(infile,idpos,SEEK_SET);
        
            if (dim == 3)
            {
                /* used to update the maximum speed in 3D cases */
                for (int i = 0; i < dim; ++i)
                    rgb_params->p_angular_velo[i] = rgb_params->angular_velo
                                    * rgb_params->rotation_dir[i];
            }
        }

        if (rgb_params->motion_type == ROTATION)
        {
            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an point:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    sprintf(msg,"Enter rotation center:");
                    CursorAfterString(infile,msg);
                    for (int i = 0; i < dim; ++i)
                    {
                        fscanf(infile,"%lf",&rgb_params->rotation_cen[i]);
                        (void) printf("%f ",rgb_params->rotation_cen[i]);
                    }
                    (void) printf("\n");
                }
                (void) fseek(infile,idpos,SEEK_SET);
            }

            if (CursorAfterStringOpt(infile,
                "Type yes if rigid body will rotate about an axis:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    /* For 2D, it is always about the z-axis */
                    if (dim == 3)
                    {
                        sprintf(msg,"Enter direction of the axis:");
                        CursorAfterString(infile,msg);
                        for (int i = 0; i < dim; ++i)
                        {
                            fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                            (void) printf("%f ",rgb_params->rotation_dir[i]);
                            mag_dir += sqr(rgb_params->rotation_dir[i]);
                        }
                        mag_dir = sqrt(mag_dir);
                        for (int i = 0; i < dim; ++i)
                            rgb_params->rotation_dir[i] /= mag_dir;
                        (void) printf("\n");
                    }
                }
            }
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (rgb_params->motion_type == FREE_MOTION ||
            rgb_params->motion_type == ROTATION)
        {
            CursorAfterString(infile,"Enter the moment of inertial: ");
            /*
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->moment_of_inertial);
                (void) printf("%f\n",rgb_params->moment_of_inertial);
            }
            else if (dim == 3)*/
            if (dim == 3)
            {
                for (int i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_moment_of_inertial[i]);
                    (void) printf("%f ",rgb_params->p_moment_of_inertial[i]);
                }
                (void) printf("\n");
            }
            (void) fseek(infile,idpos,SEEK_SET);

            CursorAfterString(infile,"Enter initial angular velocity: ");
            /*
            if (dim == 2)
            {
                fscanf(infile,"%lf",&rgb_params->angular_velo);
                (void) printf("%f\n",rgb_params->angular_velo);
            }
            else if (dim == 3)*/
            if (dim == 3)
            {
                for (int i = 0; i < dim; ++i)
                {
                    fscanf(infile,"%lf",&rgb_params->p_angular_velo[i]);
                    (void) printf("%f ",rgb_params->p_angular_velo[i]);
                }
                (void) printf("\n");
                /* initialize the euler parameters */
                rgb_params->euler_params[0] = 1.0;
                for (int i = 1; i < 4; ++i)
                    rgb_params->euler_params[i] = 0.0;
            }
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (rgb_params->motion_type == PRESET_ROTATION ||
            rgb_params->motion_type == PRESET_MOTION ||
            rgb_params->motion_type == PRESET_TRANSLATION)
        {
            rgb_params->vparams = NULL;
            rgb_params->vel_func = NULL;
            if (CursorAfterStringOpt(infile,
                "Type yes to use prescribed velocity function:"))
            {
                fscanf(infile,"%s",s);
                (void) printf("%s\n",s);
                if (s[0] == 'y' || s[0] == 'Y')
                {
                    prompt_for_velocity_func(dim,infile,rgb_params);
                }
            }
            (void) fseek(infile,idpos,SEEK_SET);
        }

        if (debugging("rgbody"))
            (void) printf("Leaving prompt_for_rigid_body_params()\n");
}       /* end prompt_for_rigid_body_params */

static void set_rgbody_params(
        RG_PARAMS* rg_params,
        HYPER_SURF* hs)
{
        int i,dim = rg_params->dim;
	body_index(hs) = rg_params->body_index;
        total_mass(hs) = rg_params->total_mass;
        mom_inertial(hs) = rg_params->moment_of_inertial;
        angular_velo(hs) = rg_params->angular_velo;
        motion_type(hs) = rg_params->motion_type;
        vparams(hs) = rg_params->vparams;
        vel_func(hs) = rg_params->vel_func;
        surface_tension(hs) = 0.0;
        for (i = 0; i < dim; ++i)
        {
            center_of_mass(hs)[i] = rg_params->center_of_mass[i];
            center_of_mass_velo(hs)[i] = rg_params->cen_of_mass_velo[i];
            rotation_center(hs)[i] = rg_params->rotation_cen[i];
            translation_dir(hs)[i] = rg_params->translation_dir[i];
            if (dim == 3)
            {
                rotation_direction(hs)[i] = rg_params->rotation_dir[i];
                p_mom_inertial(hs)[i] = rg_params->p_moment_of_inertial[i];
                p_angular_velo(hs)[i] = rg_params->p_angular_velo[i];
            }
        }
        if (dim == 3)
        {
            for (i = 0; i < 4; i++)
                euler_params(hs)[i] = rg_params->euler_params[i];
        }
}       /* end set_rgbody_params */

static void prompt_for_velocity_func(
        int dim,
        FILE *infile,
        RG_PARAMS *rgb_params)
{
        clean_up(1);

        //TODO: No need to bring in other iFluid data structures at this time..
        //      eventually will attempt to unify the iFluid and cFluid interfaces.

    /*
    FILE *infile = fopen(inname,"r");
    char s[100];
    int i;
    
    static TIME_DEPENDENT_PARAMS *td_params;

    FT_ScalarMemoryAlloc((POINTER*)&td_params,
                    sizeof(TIME_DEPENDENT_PARAMS));
    (void) printf("Available prescribed velocity function types are:\n");
    (void) printf("\tSINE:\n");
    CursorAfterString(infile,
            "Enter the prescribed velocity function type:");
    fscanf(infile,"%s",s);
    (void)printf("%s\n",s);
    if (s[0] == 's' || s[0] == 'S')
    {
        td_params->td_type = SINE_FUNC;
        CursorAfterString(infile,"Enter velocity amplitude:");
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf ",&td_params->v_amp[i]);
            (void) printf("%f ",td_params->v_amp[i]);
        }
        (void) printf("\n");
        CursorAfterString(infile,"Enter oscillation frequency:");
        fscanf(infile,"%lf ",&td_params->omega);
        (void) printf("%f\n",td_params->omega);
        td_params->omega *= 2.0*PI;
        CursorAfterString(infile,"Enter initial phase:");
        fscanf(infile,"%lf ",&td_params->phase);
        (void) printf("%f\n",td_params->phase);
        td_params->phase *= PI/180.0;
        rgb_params->vparams = (POINTER)td_params;
        rgb_params->vel_func = sine_vel_func;
    }
    else
    {
        (void) printf("Unknown type of time-dependent function!\n");
        clean_up(ERROR);
    }
    fclose(infile);
    */

}       /* end prompt_for_velocity_func */

/*
static void sine_vel_func(
        Front* front,
        POINTER vparams,
        double *coords,
        double *velo)
{
        int i;
        int dim = front->rect_grid->dim;
        TIME_DEPENDENT_PARAMS *td_params;
        double time = front->time;

        td_params = (TIME_DEPENDENT_PARAMS*)vparams;
        for (i = 0; i < dim; ++i)
            velo[i] = td_params->v_amp[i] * sin(td_params->omega*time
                        + td_params->phase);
}*/       /* end sine_vel_func */

