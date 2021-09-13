#include "iFluid.h"

#include <algorithm>


struct STAT_ENTRY
{
    double coords[MAXD];
    double val[MAXD];
};


struct BWFSTEP_LESSTHAN
{
    bool operator()(const STAT_ENTRY& A, const STAT_ENTRY& B) const
    {
        if (A.coords[0] == B.coords[0])
        {
            return A.coords[1] > B.coords[1];
        }
        return A.coords[0] < B.coords[0];
    }
};


static void printBackwardFacingStepStats(Front* front, POINTER params);





void Incompress_Solver_Smooth_Basis::printProblemSpecificStats()
{
    IF_PROB_TYPE prob_type = iFparams->prob_type;
    switch (prob_type)
    {
        case BACKWARD_FACING_STEP:
            printBackwardFacingStepStats(front,iFparams->prob_params);
            break;
        default:
            break;
    }
}

//TODO:
static void printBackwardFacingStepStats(Front* front, POINTER params)
{
    ////////////////////////////////////////////////////////////////
    printf("\n\tERROR printBackwardFacingStepStats() : Function not complete!\n");
    return;
    ////////////////////////////////////////////////////////////////

    int dim = FT_Dimension();
    if (dim != 2) return;

    static bool first = true;
    static char dirname[512];
    if (first)
    {
        sprintf(dirname,"%s/bwfstep-data",OutName(front));
        if (!create_directory(dirname,NO))
        {
            screen("Cannot create directory %s\n",dirname);
            LOC(); clean_up(EXIT_FAILURE);
        }
        first = false;
    }


    BWFSTEP_PARAMS* bfstep_params = (BWFSTEP_PARAMS*)params;
    double H = bfstep_params->H;

    std::vector<STAT_ENTRY> stat_vec;

    INTERFACE* intfc = front->interf;
    CURVE** c;
    BOND* b;

    intfc_curve_loop(intfc,c)
    {
        if (wave_type(*c) == NEUMANN_BOUNDARY) //&& is_bdry(*c))
        {
            HYPER_SURF* hs = Hyper_surf(*c);
            curve_bond_loop(*c,b)
            {
                STATE* state;
                POINT* pt = b->end;
                
                if (ifluid_comp(negative_component(hs)))
                {
                    state = (STATE*)left_state(pt);
                }
                else if (ifluid_comp(positive_component(hs)))
                {
                    state = (STATE*)right_state(pt);
                }
                else
                {
                    printf("printBackwardFacingStepStats() ERROR: "
                            "no fluid component on hypersurface\n");
                    LOC(); clean_up(EXIT_FAILURE);
                }

                STAT_ENTRY stat;
                for (int i = 0; i < dim; ++i)
                {
                    stat.coords[i] = Coords(pt)[i];
                    stat.val[i] = state->shear_force[i];
                    //TODO: using the state not working ....
                }
                stat_vec.push_back(stat);
            }
        }
    }

    std::sort(stat_vec.begin(),stat_vec.end(),BWFSTEP_LESSTHAN());
    
    std::string fname = dirname;
    fname += "/skinfriction-coef-" + std::to_string(front->step) + ".txt";
    FILE* outfile = fopen(fname.c_str(),"a");

    for (auto it = stat_vec.begin(); it != stat_vec.end(); ++it)
    {
        fprintf(outfile,"%g %g\n",it->coords[0],it->val[0]);
    }
    fclose(outfile);
}


