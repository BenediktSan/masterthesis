#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=10G
#SBATCH --job-name=2D_5P_fast
#SBATCH --output=/work/smbtsand/output/2D_5P_fast.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main2D_5P_fast.out "2D"

echo "sbatch: STOP"


