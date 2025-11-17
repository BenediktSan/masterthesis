#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=2D_5P_no_norm_xxz
#SBATCH --output=/work/smbtsand/output/2D_5P_no_norm_xxz.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main2D_no_norm_xxz.out "2D" "5"

echo "sbatch: STOP"


