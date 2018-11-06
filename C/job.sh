#!/bin/bash
#SBATCH --job-name=omp25
#SBATCH --output=res_omp25.txt
#SBATCH --error=err_omp25.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --time=24:00:00
#SBATCH --account=pr005014_taskp
#SBATCH --partition=taskp

srun ./LIF_2D_Classic_OpenMP initial2.csv 2
