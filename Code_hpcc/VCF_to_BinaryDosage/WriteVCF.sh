#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-9
#SBATCH --output=./logs/WriteVCF_ukb1_chr%a.log

#-----------------------------------------------------------------------------#
# extract gzipped vcf file to vcf file for BinaryDosage package
#-----------------------------------------------------------------------------#
REF=/auto/pmd-02/figi/HRC_extract/ukbiobank/UKB1/vcf
OUT=/staging/dvc/andreeki/BD
cd ${OUT}

# subset and write VCF file
vcftools --gzvcf ${REF}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz --recode --out ${OUT}/ukb1_chr${SLURM_ARRAY_TASK_ID}
mv ukb1_chr${SLURM_ARRAY_TASK_ID}.recode.vcf ukb1_chr${SLURM_ARRAY_TASK_ID}.vcf


