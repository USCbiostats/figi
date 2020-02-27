#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --array=1-3
#SBATCH --output=/auto/rcf-40/andreeki/FIGI_code/Code_hpcc/VCF_to_BinaryDosage/logs/WriteBD_ukb2_chr%a.log

cd /staging/dvc/andreeki/BD
time Rscript -e 'library(BinaryDosage); args <- commandArgs(trailingOnly=T); VCFtoBD(paste0(args[1], "'".vcf"'"), paste0(args[1], "'".bdose"'"))' UKB2_chr${SLURM_ARRAY_TASK_ID}
