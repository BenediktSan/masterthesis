#!/bin/bash -l
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=30 --ntasks=600 --cpus-per-task=1 --constraint=cstd01|cstd02
#SBATCH --mem=0
#SBATCH --job-name=6P_B011_pair
#SBATCH --output=/work/smbtsand/output/6P_B011_pair.txt
#SBATCH --mail-type=BEGIN,END,FAIL

echo "sbatch: START SLURM_JOB_ID $SLURM_JOB_ID on $SLURMD_NODENAME"
cd /work/smbtsand/Masterarbeit/main

mpirun main_B011_pair.out "3D" "6"

echo "sbatch: STOP"


