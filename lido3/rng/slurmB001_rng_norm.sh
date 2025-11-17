#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=16:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=B001_6P_rng_vec_norm
#SBATCH --output=/work/smbtsand/output/B001_6P_rng_vec_norm.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun mainB001_rng_vec_norm.out "3D" "6"

echo "sbatch: STOP"


