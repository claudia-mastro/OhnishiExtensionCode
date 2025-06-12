#!/bin/bash
  
#SBATCH --job-name=LM2500
#SBATCH --partition=scavenge
#SBATCH --time=1-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-500
#SBATCH --error=error.%j.err
#SBATCH --error=/home/cim24/palmer_scratch/Logs/error.%j.err
#SBATCH --output=/home/cim24/palmer_scratch/Logs/output.%j.out
module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'LM2_array.R' $SLURM_ARRAY_TASK_ID 400 25 11 
