#!/bin/csh
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o seq_validate_file
#SBATCH --export=ALL
#SBATCH --account=soc-kp
#SBATCH --partition=soc-kp

gcc -o dij_seq Dijkstra_Seq.c
./dij_seq -v

