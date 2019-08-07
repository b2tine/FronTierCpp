/***************************************************************
FronTier is a set of libraries that implements different types of 
Front Traking algorithms. Front Tracking is a numerical method for 
the solution of partial differential equations whose solutions have 
discontinuities.  

Copyright (C) 1999 by The University at Stony Brook. 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
****************************************************************/


/*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include "fabric.h"

extern void Fabric_InitFronTier(
        Front *front,
        F_BASIC_DATA *f_basic)
{
        static LEVEL_FUNC_PACK level_func_pack;
        static AF_PARAMS af_params;

        front->extra2 = (POINTER)&af_params;
        FT_ReadSpaceDomain(f_basic->in_name,f_basic);
        FT_StartUp(front,f_basic);
        FT_InitIntfc(front,&level_func_pack);
}       /* end Fabric_InitFronTier */

extern void Fabric_InitModules(
        Front *front)
{
        FILE *infile = fopen(InName(front),"r");

        CursorAfterString(infile,"Start parameters for modules");
        printf("\n");

        initParachuteModules(front);
        optimizeElasticMesh(front);
        set_equilibrium_mesh(front);

        fclose(infile);
}       /* end Fabric_InitModules */
