#!/bin/bash -l
#SBATCH --partition=short
#SBATCH --time=00:03:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=myjob
#SBATCH --output=/work/smbtsand/Masterarbeit/output/output.txt
#SBATCH --mail-type=FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
mkdir /work/smbtsand/output
cd /work/smbtsand/Masterarbeit/main

mpicxx -std=c++11 main.cpp -o main.out

echo "sbatch: STOP"