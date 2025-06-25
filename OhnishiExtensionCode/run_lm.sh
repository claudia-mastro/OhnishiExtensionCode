#!/bin/bash
  
#SBATCH --job-name=LM2_50k
#SBATCH --partition=scavenge
#SBATCH --time=12:00:00
#SBATCH --mem=20G
#SBATCH --array=1-500
#SBATCH --error=/home/cim24/palmer_scratch/Logs/error.%A_%a.err
#SBATCH --output=/home/cim24/palmer_scratch/Logs/output.%A_%a.out
module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'ls_eff_LM2.R' $SLURM_ARRAY_TASK_ID 100 25 11 
