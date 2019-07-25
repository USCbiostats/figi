#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22

# Replace original IDs with VCFIDs (ID_platform)
REF=/auto/pmd-02/figi/HRC_extract
OUT=/auto/pmd-02/figi/HRC_VCF_SampleRename

#bcftools reheader -s /auto/pmd-02/figi/andreeki/Data/FIGI_Sample_Rename/FIGI_Sample_Rename_VCF_08162018_${batch}.txt ${REF}/${batch}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -o ${OUT}/${batch}_chr${SLURM_ARRAY_TASK_ID}.vcf.gz

bcftools reheader -s /auto/pmd-02/figi/andreeki/Data/FIGI_Sample_Rename/FIGI_Sample_Rename_VCF_20190723_ukbiobank.txt ${REF}/${batch}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -o ${OUT}/${batch}_chr${SLURM_ARRAY_TASK_ID}.vcf.gz

