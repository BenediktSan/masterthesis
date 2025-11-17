#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=11:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=9P_norm
#SBATCH --output=/work/smbtsand/9P_norm.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main9P_norm.out

echo "sbatch: STOP"


