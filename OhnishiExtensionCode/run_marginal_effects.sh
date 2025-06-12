#!/bin/bash
  
#SBATCH --job-name=meffs
#SBATCH --partition=long
#SBATCH --time=14-00:00
#SBATCH --mem=20G
#SBATCH --array=1-200

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'marginal_effects.R' $SLURM_ARRAY_TASK_ID
