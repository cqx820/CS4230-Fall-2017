#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o cannon.file
#SBATCH --export=ALL
#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

module load gcc openmpi
mpicc -o cannon cannon.c -std=c99
python randomMatrix.py 1024 1024
python randomMatrix2.py 1024 1024
mpiexec -n 5 ./cannon -t 1024 1024

