#!/bin/bash
#PBS -N mpi_test
#PBS -l walltime=10:00:00
#PBS -l nodes=2:ppn=12

#The following will also load dependent modules, incl openmpi
module load R

#Change directory to the working directory
cd $PBS_O_WORKDIR

orterun -np 1 -hostfile $PBS_NODEFILE R --vanilla -f EEID_analysis.R
