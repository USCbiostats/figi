#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-2

#-----------------------------------------------------------------------------#
# extract gzipped vcf file to vcf file for BinaryDosage package
#-----------------------------------------------------------------------------#
REF=/auto/pmd-02/figi/HRC_extract/ukbiobank
OUT=/staging/dvc/andreeki/BD
cd ${OUT}

# subset and write VCF file
vcftools --gzvcf ${REF}/chr${SLURM_ARRAY_TASK_ID}.vcf.gz --recode --out ${OUT}/ukbiobank_chr${SLURM_ARRAY_TASK_ID}
mv ukbiobank_chr${SLURM_ARRAY_TASK_ID}.recode.vcf ukbiobank_chr${SLURM_ARRAY_TASK_ID}.vcf


