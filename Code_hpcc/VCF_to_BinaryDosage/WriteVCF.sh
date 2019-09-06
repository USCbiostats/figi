#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-3
#SBATCH --constraint=IB
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --output=./logs/WriteVCF_UKB2_chr%a.log


#-----------------------------------------------------------------------------#
# extract gzipped vcf file to vcf file for BinaryDosage package
#-----------------------------------------------------------------------------#
REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/BD
#cd ${OUT}

# subset and write VCF file
vcftools --gzvcf ${REF}/UKB2_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --recode --out ${OUT}/UKB2_chr${SLURM_ARRAY_TASK_ID}
mv ${OUT}/UKB2_chr${SLURM_ARRAY_TASK_ID}.recode.vcf ${OUT}/UKB2_chr${SLURM_ARRAY_TASK_ID}.vcf


