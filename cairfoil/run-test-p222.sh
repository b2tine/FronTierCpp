#!/bin/bash


MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec


$MPIEXEC -np 8 ./cairfoil -d 3 -p 2 2 2 -i $1 -o out-test-p222 &
