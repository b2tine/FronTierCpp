#!/bin/bash


MPIEXEC=/usr/local/pkg/petsc-3.13.4-dbg/bin/mpiexec


INFILE=in-DGB-Ball-vremansgs-ftsides-ssoutflow-nocollsn
OUTDIR=out-dgb-p222

$MPIEXEC -np 8 ./cairfoil -d 3 -p 2 2 2 -i $INFILE -o $OUTDIR &
