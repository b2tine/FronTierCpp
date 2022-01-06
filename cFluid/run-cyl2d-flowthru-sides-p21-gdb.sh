#!/bin/bash


# When the command below is executed, two new xterm terminal windows will appear.

# In each xterm terminal type 'run' and hit enter.
# The parallel run will not start until the 'run' command has been executed in both terminals.



MPIEXEC="/opt/petsc-3.13.4-dbg/bin/mpiexec"
#MPIEXEC="/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec"


$MPIEXEC -np 2 xterm -e gdb --args ./cFluid -d 2 -p 2 1 -i in-cyl2d-heatflux-vremansgs-flowthru-sides-coarse -o out-cyl2d-flowthru-sides-p21-gdb
