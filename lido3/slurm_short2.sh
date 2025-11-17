#!/bin/bash -l
#SBATCH --partition=med
#SBATCH --time=04:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=10G
#SBATCH --job-name=6P
#SBATCH --output=/work/smbtsand/output/6P.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main6P.out


echo "sbatch: STOP"


