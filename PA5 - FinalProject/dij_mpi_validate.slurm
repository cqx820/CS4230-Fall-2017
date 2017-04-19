#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o mpi_validate_file
#SBATCH --export=ALL
#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

module load gcc openmpi
mpicc -O2 -g -Wall -o dij_mpi Dijkstra_MPI.c
mpiexec -n 1 ./dij_mpi -v

