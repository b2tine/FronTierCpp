#!/bin/bash

#SBATCH --job-name=mpi_job_test      # Job name
#SBATCH --ntasks=2                   # Number of MPI tasks (i.e. processes
#SBATCH --cpus-per-task=1            # Number of cores per MPI task
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --output=mpi_test_%j.log     # Path to the standard output and error files relative to the working directory

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module load mvapich2/2.3.3-gcc-8.3.1
export LD_LIBRARY_PATH="/opt/petsc/petsc-3.13.4-dbg/lib:/usr/lib64:$LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH

#Need to use mpiexec when using mvapich
mpiexec -np 2 ./iFluid -d 2 -p 2 1 -i in-cyl2d-VREMAN-advterm -o out-cyl2d-test


#The command below does not work with mvapich, but should with openmpi or mpich.

    #srun --mpi=pmi2 ./iFluid -d 2 -p 2 1 -i in-cyl2d-VREMAN-advterm -o out-cyl2d-test
