#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=B111_7P_no_norm
#SBATCH --output=/work/smbtsand/output/B111_7P_no_norm.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun mainB111_no_norm.out "3D" "7"

echo "sbatch: STOP"


