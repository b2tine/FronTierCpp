#!/bin/bash


MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec
#MPIEXEC=/opt/petsc-3.13.4-dbg/bin/mpiexec


$MPIEXEC -np 2 ./parachute -d 3 -p 1 1 2 -i in-C9-nVnG-ball-VREMAN-v10-coarse-nocollsn -o out-test-p112 &
