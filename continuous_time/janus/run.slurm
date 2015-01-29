#!/bin/bash
#SBATCH --time=00:01:00
#SBATCH --job-name=test
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

module load R/r-3.1.1

mpirun -np 1 R --no-save < mpi.R 
