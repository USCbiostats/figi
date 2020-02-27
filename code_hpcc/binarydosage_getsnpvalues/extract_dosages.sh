#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

#Rscript extract_dosages.R ${chr} /staging/dvc/andreeki/GetSNPValues/aspirin_bin2/${file}

Rscript extract_dosages.R ${chr} /staging/dvc/andreeki/gwas_ld_annotation/${file}
