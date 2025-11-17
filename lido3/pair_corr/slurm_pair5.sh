#!/bin/bash -l
#SBATCH --partition=ultralong
#SBATCH --time=10-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=9P_pair
#SBATCH --output=/work/smbtsand/output/9P_pair.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main_pair.out "3D" "9"

echo "sbatch: STOP"


