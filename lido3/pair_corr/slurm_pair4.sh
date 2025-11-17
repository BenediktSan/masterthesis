#!/bin/bash -l
#SBATCH --partition=ultralong
#SBATCH --time=5-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=8P_pair
#SBATCH --output=/work/smbtsand/output/8P_pair.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main_pair.out "3D" "8"

echo "sbatch: STOP"


