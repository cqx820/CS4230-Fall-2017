#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o seq_file
#SBATCH --export=ALL
#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

gcc -o dij_seq Dijkstra_Seq.c
gcc -o generator adj_mat_gen.c
./generator 1024 1024
./dij_seq

