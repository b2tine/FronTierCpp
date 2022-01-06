#!/bin/bash



MPIEXEC="/opt/petsc-3.13.4-dbg/bin/mpiexec"
#MPIEXEC="/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec"


$MPIEXEC -np 2 ./cFluid -d 3 -p 2 1 1 -i in-ball3d-heatflux-vreman-ftsides-coarse -o out-ball3d-flowthru-sides-p211 &

