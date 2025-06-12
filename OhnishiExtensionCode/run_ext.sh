#!/bin/bash
  
#SBATCH --job-name=Y_bias
#SBATCH --partition=day
#SBATCH --time=1-00:00
#SBATCH --mem=20G

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'calc_Y_pred_GP1.R' 
