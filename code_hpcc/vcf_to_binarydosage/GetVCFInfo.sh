#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=GetVCFInfo_UKB_chr5.log

Rscript GetVCFInfo.R ukbiobank_chr1.vcf