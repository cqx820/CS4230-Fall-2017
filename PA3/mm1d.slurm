#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o mm1d.file
#SBATCH --export=ALL
#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

module load gcc openmpi
mpicc -o MM1d MM1d.cpp
python randomMatrix.py 4 4
python randomMatrix2.py 4 4
mpiexec -n 2 ./MM1d

