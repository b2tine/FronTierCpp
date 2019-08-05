#include "rigidbody.h"


static void prompt_for_rigid_body_params(int,char*,RG_PARAMS*); 
static void set_rgbody_params(RG_PARAMS*,HYPER_SURF*);
static void prompt_for_velocity_func(int,char*,RG_PARAMS*);
//static void sine_vel_func(Front*,POINTER,double*,double*);



EXPORT void rgb_init(Front*front, RG_PARAMS* rgb_params)
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

}       /* end rgb_init */

static void prompt_for_rigid_body_params(
        int dim,
        char* inname,
        RG_PARAMS* rgb_params)
{
    int i;
    char msg[100],s[100],ss[100];
    FILE *infile = fopen(inname,"r");
    boolean is_preset_motion = NO;
    double mag_dir;
    static int count = 1;

    if (debugging("rgbody"))
        (void) printf("Enter prompt_for_rigid_body_params()\n");

    if( count == 1 )
    {
        rgb_params->dim = dim;
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
        (void) printf("\tFREE_MOTION:\n");
        (void) printf("\tCOM_MOTION (center of mass):\n");
        (void) printf("\tTRANSLATION:\n");
        (void) printf("\tROTATION:\n");
        CursorAfterString(infile,"Enter type of dynamic motion: ");
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
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&rgb_params->translation_dir[i]);
            (void) printf("%f ",rgb_params->translation_dir[i]);
            mag_dir += sqr(rgb_params->translation_dir[i]);
        }
        (void) printf("\n");
        mag_dir = sqrt(mag_dir);
        for (i = 0; i < dim; ++i)
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
        for (i = 0; i < dim; ++i)
        {
            fscanf(infile,"%lf",&rgb_params->center_of_mass[i]);
            (void) printf("%f ",rgb_params->center_of_mass[i]);
        }
        (void) fseek(infile,idpos,SEEK_SET);

        (void) printf("\n");
        sprintf(msg,"Enter the initial center of mass velocity:");
        CursorAfterString(infile,msg);
        for (i = 0; i < dim; ++i)
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
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                (void) printf("%f ",rgb_params->rotation_dir[i]);
                mag_dir += sqr(rgb_params->rotation_dir[i]);
            }
            (void) printf("\n");
            mag_dir = sqrt(mag_dir);
            for (i = 0; i < dim; ++i)
                rgb_params->rotation_dir[i] /= mag_dir;
            /* initialize the euler parameters */
            rgb_params->euler_params[0] = 1.0;
            for (i = 1; i < 4; ++i)
                rgb_params->euler_params[i] = 0.0;
        }
        (void) fseek(infile,idpos,SEEK_SET);

        /* Center of axis is the coordinate of a point on the axis */
        CursorAfterString(infile,"Enter rotation center:");
        for (i = 0; i < dim; ++i)
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
            for (i = 0; i < dim; ++i)
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
                for (i = 0; i < dim; ++i)
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
                    for (i = 0; i < dim; ++i)
                    {
                        fscanf(infile,"%lf",&rgb_params->rotation_dir[i]);
                        (void) printf("%f ",rgb_params->rotation_dir[i]);
                        mag_dir += sqr(rgb_params->rotation_dir[i]);
                    }
                    mag_dir = sqrt(mag_dir);
                    for (i = 0; i < dim; ++i)
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
            for (i = 0; i < dim; ++i)
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
            for (i = 0; i < dim; ++i)
            {
                fscanf(infile,"%lf",&rgb_params->p_angular_velo[i]);
                (void) printf("%f ",rgb_params->p_angular_velo[i]);
            }
            (void) printf("\n");
            /* initialize the euler parameters */
            rgb_params->euler_params[0] = 1.0;
            for (i = 1; i < 4; ++i)
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
                prompt_for_velocity_func(dim,inname,rgb_params);
            }
        }
        (void) fseek(infile,idpos,SEEK_SET);
    }
    fclose(infile);

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
        center_of_mass_velo(hs)[i] =
                            rg_params->cen_of_mass_velo[i];
        rotation_center(hs)[i] =
                            rg_params->rotation_cen[i];
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
        char *inname,
        RG_PARAMS *rgb_params)
{
    checkpoint("Not Implemented, Exiting");
    clean_up(1);

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

