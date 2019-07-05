#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=12,15

#Rscript job01.R ${SLURM_ARRAY_TASK_ID} GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_vcfid GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr${SLURM_ARRAY_TASK_ID}

Rscript job01.R ${SLURM_ARRAY_TASK_ID} GetSNPValues_asp_ref_age_ref_imp_sex_study_gxe_PC1-3_N_72820_vcfid GetSNPValues_asp_ref_age_ref_imp_sex_study_gxe_PC1-3_N_72820_chr${SLURM_ARRAY_TASK_ID}
