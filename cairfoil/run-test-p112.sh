#!/bin/bash


#MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec
MPIEXEC=/opt/petsc-3.13.4-dbg/bin/mpiexec


$MPIEXEC -np 2 ./cairfoil -d 3 -p 1 1 2 -i in-C9-Vent-Ball-vremansgs-ftsides-no-collsn-test -o out-test-p112 &
