#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=2
#SBATCH --output=./logs/BinaryDosage_ukbiobank_chr%a.log

time Rscript -e 'library(BinaryDosage); args <- commandArgs(trailingOnly=T); VCFtoBD(paste0(args[1], "'".vcf"'"), paste0(args[1], "'".bdose"'"))' ukbiobank_chr${SLURM_ARRAY_TASK_ID}
