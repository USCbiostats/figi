#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-17
#SBATCH --output=./logs/WriteVCF_reach_chr%a.log

#-----------------------------------------------------------------------------#
# extract gzipped vcf file to vcf file for BinaryDosage package
#-----------------------------------------------------------------------------#
REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/BD
#cd ${OUT}

# subset and write VCF file
vcftools --gzvcf ${REF}/reach_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --recode --out ${OUT}/reach_chr${SLURM_ARRAY_TASK_ID}
mv ${OUT}/reach_chr${SLURM_ARRAY_TASK_ID}.recode.vcf ${OUT}/reach_chr${SLURM_ARRAY_TASK_ID}.vcf


