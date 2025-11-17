#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=1-12:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=5G
#SBATCH --job-name=B111_6P_rng_vec_no_norm.out
#SBATCH --output=/work/smbtsand/output/B111_6P_rng_vec_no_norm.out
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun mainB111_6P_rng_vec_no_norm.out

echo "sbatch: STOP"
