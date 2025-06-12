#!/bin/bash
  
#SBATCH --job-name=CADEA_bias
#SBATCH --partition=scavenge
#SBATCH --time=01:00:00
#SBATCH --mem=5G
#SBATCH --array=1-500

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'CASES_bias_LM.R' $SLURM_ARRAY_TASK_ID 20 25 11
