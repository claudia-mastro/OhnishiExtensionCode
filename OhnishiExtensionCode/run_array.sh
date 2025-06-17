#!/bin/bash
  
#SBATCH --job-name=GP1
#SBATCH --partition=week
#SBATCH --time=7-00:00:00
#SBATCH --mem=500G
#SBATCH --array=1
#SBATCH --error=/home/cim24/palmer_scratch/Logs/error.%A_%a.err
#SBATCH --output=/home/cim24/palmer_scratch/Logs/output.%A_%a.out

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'GP1_array.R' $SLURM_ARRAY_TASK_ID 

