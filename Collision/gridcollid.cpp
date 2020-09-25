#include "gridcollid.h"

//Default values: bdry = false 
void getCandidates(INTERFACE* intfc, bool bdry)
{
    //Table table;
    SURFACE***** surfs;
    BOND***** bonds;
    TRI***** tris;

    int*** num_bonds;
    int*** num_tris;


    boolean status;
    if (intfc->modified || intfc->table->new_grid)
    {
        status = make_tri_lists(intfc);
        if (status != YES)
        {
            printf("make_tri_lists() ERROR"); LOC();
            clean_up(EXIT_FAILURE);
            //return NO;
        }
        
        status = make_bond_lists3d(intfc);
        if (status != YES)
        {
            printf("make_bond_lists3d() ERROR"); LOC();
            clean_up(EXIT_FAILURE);
            //return NO;
        }
        intfc->modified = NO;
        intfc->table->new_grid = NO;
    }

    int xmax = topological_grid(intfc).gmax[0];
    int ymax = topological_grid(intfc).gmax[1];
    int zmax = topological_grid(intfc).gmax[2];

/*
*	The information that follows is set up in the interface table defined
*	in int.h, and is updated by the function make_tri_lists(), which "bins"
*	the surface elements (TRI's) into grid cubes determined by the
*	topological grid.
*/

	surfs = intfc->table->surfaces;
    bonds = intfc->table->bonds;
    tris = intfc->table->tris;
	num_tris = intfc->table->num_of_tris;
	comps = intfc->table->compon3d;

/*
*	The following code sets up a loop over all of the "bins" (i.e. grid
*	cubes, or blocks of the topological grid.)
*/

	for (iz = 0; iz < zmax; ++iz)
	for (iy = 0; iy < ymax; ++iy)
	for (ix = 0; ix < xmax; ++ix)
	{
	    /* This grid cube must be passed through by the front...  */

	    if (comps[iz][iy][ix] != ONFRONT)
		continue;

	    /* and it must contain 2 or more TRI's.     */
	    
	    int nt = num_tris[iz][iy][ix];
	    if (nt < 2)
		continue;

	    /*
	    *    Set t,s equal to, respectively, the list of TRI's and the
	    *    list of SURFACE's for this grid cube.
	    */
	    
	    TRI** t = tris[iz][iy][ix];
	    SURFACE** s = surfs[iz][iy][ix];
	    
	        /* Loop over all pairs of Tris: */
	    for (int i = 0, t0 = t, s0 = s; i < nt - 1; ++i, ++t0, ++s0)
	    {

		if (is_bdry(*s0) && !bdry)
		    continue;

		for (int j = i+1, t1 = t+j, s1 = s+j; j < nt; ++j, ++t1, ++s1)
		{

		    if (*t1 == Tri_on_side01(*t0) ||
		        *t1 == Tri_on_side12(*t0) ||
                *t1 == Tri_on_side20(*t0)) continue;
		    ct0 = *t0;  ct1 = *t1;
		    cs0 = *s0;  cs1 = *s1;
		    
		    if (is_bdry(cs1) && !bdry)
			continue;
		}

}


//sort_tris_on_blocks(intfc);
//
//
//
