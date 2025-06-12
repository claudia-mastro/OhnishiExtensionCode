#!/bin/bash
  
#SBATCH --job-name=ceffs
#SBATCH --partition=scavenge
#SBATCH --time=1-00:00
#SBATCH --mem=20G
#SBATCH --array=1-200

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'conditional_effects.R' $SLURM_ARRAY_TASK_ID
