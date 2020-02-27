#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

#Rscript calculate_r2.R  huygue_gwas_140_chr1_TMP X1.55246035
Rscript calculate_r2.R ${chr} ${bp}

