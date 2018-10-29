#if !defined(_FCELL_H)
#define _FCELL_H

#include <front/fdecs.h>
struct _CELL_INFO_3D{
	// Input variables 
	boolean is_corner;
	boolean is_side;
	boolean is_face;
	boolean is_soln_solute;
	boolean is_center_solute;
	boolean unstable_cell;
	double solute_volume;
	int ix[2][2][2];
	int iy[2][2][2];
	int iz[2][2][2];
	int comps[2];
	int nv[2];
	int num_comps;
	int comp[2][2][2];
	int icrds[2][2][2][MAXD];
	double crds[2][2][2][MAXD];
	int soln_comp;
	int cell_comp;
	double crx_coords[3][2][2][MAXD];
	
	// output variables 
	boolean nb_flag[3][3][3];
	double	nb_frac[3][3][3];
	double area[3][3][3];
	double orphan;
	double full_cell_vol;
};
typedef struct _CELL_INFO_3D CELL_INFO_3D;



#endif
