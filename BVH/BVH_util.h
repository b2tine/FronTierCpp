#ifndef BVH_UTIL_H
#define BVH_UTIL_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
//NOTE: CGAL/Surface_mesh.h has macro conflict with the 
//      fdecs.h macro "node_type"

#include "BVH.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include <limits>
#include <utility>
#include <map>


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using cgalPoint3 = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<cgalPoint3>;
using Vertex_index = Mesh::Vertex_index;
using Face_index = Mesh::Face_index;


void TriMeshOFF2MonoCompSurf(Front*, Mesh*);
void TriMeshOFF2Surf(INTERFACE*,COMPONENT, COMPONENT, Mesh*, SURFACE**);

std::pair<std::vector<double>,std::vector<double> >
getInputMeshDimensionsWithPad(Mesh* mesh, double pad);



#endif
