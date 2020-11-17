#include "keps.h"

//TODO: eventually move into 2d and subclasses

void KE_CARTESIAN::implicitComputeK2d(COMPONENT sub_comp)
{
    int icoords[MAXD], icoords_nb[MAXD];
    double crx_coords[MAXD];
	double nor[MAXD];
    COMPONENT comp;
    
    const GRID_DIRECTION dir[3][2] = {
        {WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}
    };

    //double *gamma = field->gamma;
	double *Pk = field->Pk;
    double *K = field->k;
	double *mu_t = field->mu_t;

	double rho = eqn_params->rho;
	double nu = eqn_params->mu/eqn_params->rho;
    double Cmu = eqn_params->Cmu;
    double C1 = eqn_params->C1;
    double C2 = eqn_params->C2;
    double Cbc = eqn_params->Cbc;
	double delta_k = eqn_params->delta_k;
	double delta_eps = eqn_params->delta_eps;
    double l0 = eqn_params->l0;
	
	INTERFACE* grid_intfc = front->grid_intfc;
	HYPER_SURF_ELEMENT *hse;
	HYPER_SURF *hs;
	STATE *intfc_state;


    start_clock("computeK");
    if (debugging("trace")) printf("Entering implicitComputeK2d()\n");

    setIndexMap(sub_comp);

    if (debugging("trace"))
    {
        //TODO: What was the purpose of this check?
        //
        //      Appears to be an unfinished feature.
        //      should domain_size be getting multiplied by top_g[i] or something similar??
        int domain_size = 1;
        printf("ilower = %d  iupper = %d\n",ilower,iupper);
        for (int i = 0; i < dim; ++i)
            printf("domain_size = %d\n",domain_size);
    }

    PETSc solver;
    solver.Create(ilower, iupper-1, 5, 5);

    start_clock("set_coefficients");

    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        int ic = d_index2d(i,j,top_gmax);
        comp = top_comp[ic];

        if (comp != sub_comp) continue;
        
        icoords[0] = i;
        icoords[1] = j;
        int I = ij_to_I[i][j];

        double aII = 0.0;
        double RHS = 0.0;

        //gammaK updated using previous values (time n)
        double nu_t = mu_t[ic]/rho;
        double gammaK = std::max(Cmu*K[ic]/nu_t,0.0);

        Karray[ic] = (1.0 - gammaK*m_dt)*K[ic] + Pk[ic]*m_dt;

        for (int l = 0; l < dim; ++l)
        {
            double alpha = 0.5*m_dt*field->vel[l][ic]/top_h[l];
            double beta = m_dt*nu_t/delta_k/sqr(top_h[l]);
                //double beta = m_dt*(nu + nu_t/delta_k)/sqr(top_h[l]);
        
            Karray[ic] -= 2.0*beta*K[ic];
        
            for (int nb = 0; nb < 2; ++nb)
            {
                for (int n = 0; n < dim; ++n)
                    icoords_nb[n] = icoords[n];

                int sign = pow(-1,nb);

                icoords_nb[l] =
                    (nb == 0) ? icoords[l] - 1 : icoords[l] + 1;

                int icnb = d_index(icoords_nb,top_gmax,dim);

                crx_status = FT_StateStructAtGridCrossing(front,
                        grid_intfc,icoords,dir[l][nb],comp,
                        (POINTER*)&intfc_state,&hs,crx_coords);
                
                if (!crx_status) 
                {
                    //no cross
                    Karray[ic] += sign*alpha*K[icnb] + beta*K[icnb];
                    solver.Add_A(I,I_nb,coeff_nb);
                }
                else if (wave_type(hs) == NEUMANN_BOUNDARY ||
                         wave_type(hs) == MOVABLE_BODY_BOUNDARY)
                {
                    double v_tan[MAXD];
                    setSlipBoundary(icoords,l,nb,comp,hs,intfc_state,field->vel,v_tan);

                    double mag_vtan = Magd(v_tan,dim);

                    //friction velocity
                    double u_t = std::max(
                            pow(eqn_params->Cmu,0.25)*sqrt(std::max(field->k[ic],0.0)),
                                            mag_vtan/eqn_params->y_p);
                    
                    //double u_t = 0.41*mag_vtan/log(eqn_params->E*eqn_params->y_p);
                    double K_nb = u_t*u_t/sqrt(eqn_params->Cmu);
                    Karray[ic] += sign*alpha*K_nb + beta*K_nb;

                    //RHS += ...
                }
                else if (wave_type(hs) == DIRICHLET_BOUNDARY)
                {
                    // Inlet/Outlet on rectangular boundaries only
                    if (boundary_state_function(hs) &&
                    strcmp(boundary_state_function_name(hs),
                                    "flowThroughBoundaryState") == 0)
                    {
                        //OUTLET: nor dot grad(k) = 0
                        double K_nb = K[ic];
                        Karray[ic] += sign*alpha*K_nb + beta*K_nb;
                    }
                    else
                    {
                        //INLET: k = Cbc*|u|^2 where Cbc \in [0.003,0.01]
                        double sqrmag_vel = Dot2d(intfc_state->vel,intfc_state->vel);
                        double K_nb = eqn_params->Cbc*sqrmag_vel;
                        Karray[ic] += sign*alpha*K_nb + beta*K_nb;
                    }

                    //RHS += ...
                }
                else
                {
                    printf("UNKNOWN BOUNDARY TYPE!\n"); LOC();
                    clean_up(EXIT_FAILURE);
                }
            
            }//nb

        }//l
    
        solver.Add_A(I,I,aII);
        solver.Add_b(I,RHS);

    }//i,j
            

    stop_clock("set_coefficients");

    
    double* x;
    int num_iter = 0;
    double rel_residual = 0;

    solver.SetMaxIter(500);
    solver.SetTol(1e-8);
    
    start_clock("petsc_solve");
    solver.Solve();
    stop_clock("petsc_solve");
    
    solver.GetNumIterations(&num_iter);
    solver.GetFinalRelativeResidualNorm(&rel_residual);

    if (debugging("PETSc"))
    {
        printf("KE_CARTESIAN::implicitComputeK2d(): "
                "num_iter = %d, rel_residual = %g \n",num_iter,rel_residual);
    }

    
    //TODO: not parallelized intentionally?
        //FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));

    int size = iupper - ilower;
    FT_VectorMemoryAlloc((POINTER*)&x,size,sizeof(double));

    solver.Get_x(x);

    start_clock("scatter_data");
    for (int j = jmin; j <= jmax; ++j)
    for (int i = imin; i <= imax; ++i)
    {
        I = ij_to_I[i][j];
        ic = d_index2d(i,j,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
            Karray[ic] = x[I-ilower];
        else
            Karray[ic] = 0.0;
    }

	FT_ParallelExchGridArrayBuffer(Karray,front,NULL);
    
    for (int j = 0; j <= top_gmax[1]; ++j)
    for (int i = 0; i <= top_gmax[0]; ++i)
    {
        ic = d_index2d(i,j,top_gmax);
        comp = cell_center[ic].comp;
        if (comp == sub_comp)
            K[ic] = Karray[ic];
    }
    
    stop_clock("scatter_data");
    
    FT_FreeThese(1,x);

    if (debugging("trace")) printf("Leaving implicitcomputeK2d()\n");
    stop_clock("computeK");
}































