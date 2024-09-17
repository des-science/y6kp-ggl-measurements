#!/bin/bash 
#SBATCH -A des 
#SBATCH -C cpu 
#SBATCH -q regular 
#SBATCH -t 3:00:00 
#SBATCH --nodes=1

module load python
source activate sompz

srun python run_ggl.py
