#!/bin/bash
  
#SBATCH --job-name=ceff_bias
#SBATCH --partition=scavenge
#SBATCH --time=1:00:00
#SBATCH --mem=20G

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'effect_bias_LM.R' 
