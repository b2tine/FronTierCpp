#include "collid.h"



double CollisionSolver3d::compute_acceleration(
         SPRING_VERTEX *sv)
{
    for (int k = 0; k < dim; ++k)
        accel[k] = 0.0;

    //acceleration due to elastic stretching force
    for (int j = 0; j < sv->num_nb; ++j)
    {
        double len = 0.0;
        double vec[MAXD];
        for (int k = 0; k < m_dim; ++k)
        {
            vec[k] = sv->x_nb[j][k] - sv->x[k];
            len += vec[k]*vec[k];
        }
        len = sqrt(len);

        double dL = len - sv->len0[j];

        //zero compressive stress
        if (dL <= 0.0) continue;

        for (int k = 0; k < dim; ++k)
        {
            accel[k] += sv->k[j]*dL*vec[k]/len/sv->m;
        }
    }

    //acceleration due to elastic bending and damping forces
    for (int k = 0; k < dim; ++k)
    {
        accel[k] += sv->bendforce[k]/sv->m;
        accel[k] -= sv->lambda*(sv->v[k] - sv->ext_impul[k])/sv->m;
    }

    for (int k = 0; k < dim; ++k)
    {
        sv->f[k] = accel[k]*sv->m;
    }

    for (int k = 0; k < dim; ++k)
    {
        accel[k] += sv->ext_accel[k] + sv->fluid_accel[k] + sv->other_accel[k];
    }

}   /* end compute_acceleration */



