#include "implicit.h"

void WenoFlux(
	int mesh_size, 
	double *u_old, 
	double *cflux,
	double dx, 
	double dt) 
{
	//double *u1, *u2, *u3;
	double *flux;
	int i;
	int nrad = 3; /* 5th-order weno */
	int extend_size = mesh_size + 2 * nrad;
	double *u_extend;
	    //double lambda = -dt/dx;

	//FT_VectorMemoryAlloc((POINTER*)&u1,mesh_size,sizeof(double));
	//FT_VectorMemoryAlloc((POINTER*)&u2,mesh_size,sizeof(double));
	//FT_VectorMemoryAlloc((POINTER*)&u3,mesh_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&u_extend,extend_size,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&flux,mesh_size,sizeof(double));

	/* Set the value on extended mesh */
	for (i = 0; i < mesh_size; i++) {
	    u_extend[i + nrad] = u_old[i];
	}
	for (i = 1; i <= nrad; i++) {
	    u_extend[nrad - i] = u_old[mesh_size - 1 - i];
	    u_extend[nrad + mesh_size - 1 + i] = u_old[i];
	}

	Weno5_Get_Flux(u_extend,flux,dx,mesh_size);
	for (i = 0; i < mesh_size; i++)
    {
		cflux[i] = flux[i];
    }

    /*
    //u^1
	Weno5_Get_Flux(u_extend,flux,dx,mesh_size);
	for (i = 0; i < mesh_size; i++)
    {
		u1[i] = u_old[i] + 0.5*dt*flux[i];
    }
    */
	
    //FT_FreeThese(5,u1,u2,u3,u_extend,flux);

    FT_FreeThese(2,u_extend,flux);
}

void implicitSolver(
	int mesh_size, 
	double *u_old, 
	double *u_new,
	double *source,//explicitly computed convective flux
	double dx, 
	double dt) 
{
    //temp hard coded variables for prototyping
    double u_left = 0.0;
    double u_right = 0.0;
    double mu = 0.1;

    int ilower = 1;
    int iupper = mesh_size;

    PETSc solver;
    solver.Create(ilower,iupper,3,3);
    solver.Reset_A();
    solver.Reset_b();
    solver.Reset_x();

    for (int i = 1; i < mesh_size; ++i)
    {
        int I = i;
    
        double alpha = mu*dt/dx/dx;
        double beta = 1.0 + 2.0*mu*dt/dx/dx;
        double rhs = u_old[i] - source[i]*dt;

        int I_nb;
        double coeff_nb = -alpha;
        double coeff = beta;

        for (int nb = 0; nb < 2; ++nb)
        {
            I_nb = (nb == 0) ? I-1 : I+1;

            if (I_nb == 0)
            {
                //left side dirichlet condition
                rhs -= coeff_nb*u_left;
            }
            else if (I_nb == mesh_size)
            {
                //right side dirichlet condition
                rhs -= coeff_nb*u_right;
            }
            else
            { 
                //no boundary
                solver.Add_A(I,I_nb,coeff_nb);
            }
        }
            
        solver.Add_A(I,I,beta);
        solver.Add_b(I,rhs);
    }


    solver.SetMaxIter(500);
    solver.SetTol(1.0e-08);
    solver.Solve();

    double* sol = new double[mesh_size];
    solver.Get_x(sol);

    //u_new[0] = u_left;
    for (int i = 0; i < mesh_size; ++i)
    {
        u_new[i] = sol[i];
    }
    //u_new[mesh_size-1] = u_right;

    delete[] sol;
}

