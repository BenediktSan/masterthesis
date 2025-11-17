#!/bin/bash -l
#SBATCH --partition=ultralong
#SBATCH --time=3-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=B001_9P_no_norm
#SBATCH --output=/work/smbtsand/output/B001_9P_no_norm.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun mainB001_no_norm.out "3D" "9"

echo "sbatch: STOP"


