#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=B111_6P_no_norm_long
#SBATCH --output=/work/smbtsand/output/B111_6P_no_norm_long.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun mainB111_no_norm_long.out "3D" "6"

echo "sbatch: STOP"


