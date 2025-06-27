#!/bin/bash
  
#SBATCH --job-name=GP2_bias
#SBATCH --partition=scavenge
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH --array=1-500
#SBATCH --error=/home/cim24/palmer_scratch/Logs/error.%A_%a.err
#SBATCH --output=/home/cim24/palmer_scratch/Logs/output.%A_%a.out

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'effect_bias.R' $SLURM_ARRAY_TASK_ID 

