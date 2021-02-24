#!/bin/bash

/opt/petsc-3.13.4-dbg/bin/mpirun -np 2 ./fabric_Cplus -d 3 -p 2 1 1 -i in-fabric-on-ball-pptest -o out-pptest &
