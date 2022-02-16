#!/bin/bash


MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec


INFILE=$1
OUTDIR=$2

$MPIEXEC -np 8 ./cairfoil -d 3 -p 2 2 2 -i $INFILE -o $OUTDIR &
