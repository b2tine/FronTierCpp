#include <FronTier.h>

static void surf_com_translation(SURFACE*,double*);
static void surf_enlargement(SURFACE*,double);

static void initSingleRigidBody(FILE*,Front*);
static void initMultiRigidBodies(FILE*,Front*,int);
static void init_rigid_sphere(FILE*,Front*);
static void init_rigid_box(FILE*,Front*);
static void init_rigid_human(FILE*,Front*);
static void init_rigid_cylinder(FILE*,Front*);


void initRigidBody(
	Front *front)
{
	FILE *infile = fopen(InName(front),"r");
        char string[100];
	int num_rgb;

	if (CursorAfterStringOpt(infile,"Enter yes to add rigid body:"))
	{
	    fscanf(infile,"%s",string);
	    (void) printf("%s\n",string);
	    if (string[0] != 'y' && string[0] != 'Y')
		return;
	}
	else
	    return;

	num_rgb = 1;
	if (CursorAfterStringOpt(infile,"Enter the number of rigid bodies:"))
	{
	    fscanf(infile,"%d",&num_rgb);
	    (void) printf("%d\n",num_rgb);
	}

	if (num_rgb == 1)
	    initSingleRigidBody(infile,front);
	else
	    initMultiRigidBodies(infile,front,num_rgb);
	fclose(infile);
}

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
	case 'h':
	case 'H':
	    init_rigid_human(infile, front);
	    break;
	case 'c':
	case 'C':
	    init_rigid_cylinder(infile, front);
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
        pos_comp = LIQUID_COMP2;
        if (CursorAfterStringOpt(infile,
                            "Type yes if the rigid body is fixed:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                w_type = NEUMANN_BOUNDARY;
        }
        FT_MakeEllipticSurf(front,cen,radii,neg_comp,pos_comp,
                            w_type,2.0,&surf);
	return;
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
        if (CursorAfterStringOpt(infile,
                            "Type yes if the rigid body is fixed:"))
        {
            fscanf(infile,"%s",string);
            (void) printf("%s\n",string);
            if (string[0] == 'y' || string[0] == 'Y')
                w_type = NEUMANN_BOUNDARY;
        }
        FT_MakeCuboidSurf(front,cen,edge,neg_comp,pos_comp,w_type,2,&surf);
	return;
}	/* end init_rigid_box */

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
	return;
}	/* end init_rigid_human */

static void init_rigid_cylinder(
	FILE *infile,
	Front *front)
{
	char string[100];
	double cen[MAXD], radius, height;
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
	(void) printf("Rigid body can be fixed (F) or Movable (M)\n");
	(void) printf("The default is Movable (M)\n");
	wave_type(surf) = MOVABLE_BODY_BOUNDARY;
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
	FT_MakeCylinderSurf(front,cen,radius,height/2,2,neg_comp,pos_comp,
					w_type,&surf);
	return;
}	/* end init_rigid_cylinder */

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
