#!/bin/bash -l
#SBATCH --partition=ultralong
#SBATCH --time=14-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 
#SBATCH --mem=0
#SBATCH --job-name=2D_10P_no_norm_long
#SBATCH --output=/work/smbtsand/output/2D_10P_no_norm_long.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main2D_no_norm_long.out "2D" "10"

echo "sbatch: STOP"


