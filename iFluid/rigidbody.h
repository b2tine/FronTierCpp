#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include <FronTier.h>

struct RG_PARAMS
{
    int dim;
	boolean no_fluid;		        /* For benchmark tests */
	int	body_index;		            /* Body index */
    double  total_mass;             /* Total mass */
    double  moment_of_inertial;     /* Moment of inertial about the axis */
    double  center_of_mass[MAXD];   /* Center of mass */
    double  rotation_dir[MAXD];     /* Direction of rotation */
	double	translation_dir[MAXD];	    /* Restricted direction of motion */
    double  rotation_cen[MAXD];     /* Center of rotation */
    double  cen_of_mass_velo[MAXD]; /* Center of mass velocity */
    double  angular_velo;           /* Angular velocity of rotation */
	double  p_moment_of_inertial[MAXD];
	double  p_angular_velo[MAXD];
	double  euler_params[4];
    double  old_euler_params[4];
	void	(*vel_func)(Front*,POINTER,double*,double*);
	POINTER vparams;
    MOTION_TYPE motion_type;
};


// rigidbody.cpp
extern void initRigidBody(Front *front);
extern void rgb_init(Front*,RG_PARAMS*);
extern void resetRigidBodyVelocity(Front *front);
//extern void rgb_modification(Front*,RG_PARAMS*);
extern void printRigidBodyMeshQuality(Front *front);


#endif
