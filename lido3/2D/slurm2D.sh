#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --time=06:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=2D_5P_no_norm.out
#SBATCH --output=/work/smbtsand/output/B2D_5P_no_norm.out
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main2D_no_norm.out "2D" "5"

echo "sbatch: STOP"


