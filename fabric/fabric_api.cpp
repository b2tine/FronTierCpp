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

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <string>


#ifdef __cplusplus
extern "C" {
#endif


//void Fabric_Init(int argc, char* argv[])
void Fabric_Init(char* inname)
{
    checkpoint("Entering Fabric_Init()");

    static Front front;
	static F_BASIC_DATA f_basic;

    std::ifstream infile(inname);
    std::vector<std::string> argstrings;

    //TODO: need to make sure file is open
    //      getting stuck in read loop
    checkpoint(1);
    int nargs = 0;
    while (!infile.eof())
    {
        std::string curr;
        infile >> curr;
        argstrings.push_back(curr);
        nargs++;
    }
    checkpoint(2);

    argstrings.insert(argstrings.begin(),"dummyarg");
    
    char* args[nargs];
    for (int i = 0; i < nargs; ++i)
    {
        args[i] = new char[argstrings[i].length()+1];
        std::strcpy(args[i],argstrings[i].c_str());
    }
    
    for (int i = 0; i < nargs; ++i)
    {
        std::cout << args[i] << " ";
    }
    std::cout << "\n";
    
    checkpoint("Leaving Fabric_Init()");
    return;

	//FT_Init(argc,argv,&f_basic);

	/* Initialize basic computational data */

        //Fabric_InitFronTier(&front,&f_basic);
        //Fabric_InitModules(&front);
	    
        //FT_Draw(&front);

	//clean_up(0);
}


#ifdef __cplusplus
}
#endif


