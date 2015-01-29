#!/bin/bash
#PBS -N mpi_test
#PBS -l walltime=00:01:00
#PBS -l nodes=1:ppn=1

#The following will also load dependent modules, incl openmpi
module load R/r-3.0.1

#Change directory to the working directory
cd $PBS_O_WORKDIR

orterun -np 1 -hostfile $PBS_NODEFILE R --vanilla -f EEID_analysis.R
