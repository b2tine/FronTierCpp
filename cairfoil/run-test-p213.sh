#!/bin/bash


MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec
#MPIEXEC=/opt/petsc-3.13.4-dbg/bin/mpiexec


$MPIEXEC -np 6 ./cairfoil -d 3 -p 2 1 3 -i $1 -o out-test-p213 &
