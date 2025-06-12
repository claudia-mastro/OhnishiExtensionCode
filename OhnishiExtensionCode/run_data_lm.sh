#!/bin/bash
  
#SBATCH --job-name=LM_emploidur
#SBATCH --partition=scavenge
#SBATCH --time=1-00:00
#SBATCH --mem=50G

module load R-bundle-Bioconductor/3.18-foss-2022b-R-4.3.2

srun Rscript ./'LM_array_data.R'
