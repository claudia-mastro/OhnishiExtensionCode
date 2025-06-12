#!/bin/bash
  
#SBATCH --job-name=GP1
#SBATCH --partition=week
#SBATCH --time=3-12:00:00
#SBATCH --mem=80G
#SBATCH --array=1-100

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'GP1_array.R' $SLURM_ARRAY_TASK_ID 

