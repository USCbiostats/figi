#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=16-22
#SBATCH --output=./logs/vcftools_ukbiobank_chr%a.log

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/bdose
OUTF=/auto/pmd-02/figi/HRC_BDose

cd ${OUT}

# subset and write VCF file
vcftools --gzvcf ${REF}/ukbiobank_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --snps ./ukbiobank/ukbiobank_subset_chr${SLURM_ARRAY_TASK_ID}.txt --recode --out ${OUT}/ukbiobank_chr${SLURM_ARRAY_TASK_ID}
mv ukbiobank_chr${SLURM_ARRAY_TASK_ID}.recode.vcf ukbiobank_chr${SLURM_ARRAY_TASK_ID}.vcf


#Rscript -e 'library(BinaryDosage); args <- commandArgs(trailingOnly=T); VCFtoBD(paste0(args[1], "'".vcf"'"), paste0(args[1], "'"_newest.bdose"'"))' ukbiobank_chr5

#Rscript -e 'library(BinaryDosage); args <- commandArgs(trailingOnly=T); VCFtoBD(paste0(args[1], "'".vcf"'"), paste0(args[1], "'".bdose"'"))' $1_chr$2

#rm ${OUT}/$1_chr$2.vcf
#mv ${OUT}/$1_chr$2.bdose ${OUTF}