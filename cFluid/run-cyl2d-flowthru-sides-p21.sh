#!/bin/bash



MPIEXEC="/opt/petsc-3.13.4-dbg/bin/mpiexec"
#MPIEXEC="/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec"


$MPIEXEC -np 2 ./cFluid -d 2 -p 2 1 -i in-cyl2d-heatflux-vremansgs-flowthru-sides-coarse -o out-cyl2d-flowthru-sides-p21 &
