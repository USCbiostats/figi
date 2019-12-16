#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22
#SBATCH --output=../logs/Plink_aspirin_ldclump_chiSqGxE_chr%a.log

INN=/auto/pmd-01/andreeki/clump
INN=/staging/dvc/andreeki/clump
OUT=/staging/dvc/andreeki/clump_results/aspirin_72269

plink \
    --bed ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.bed \
    --bim ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}_fix.bim \
    --fam ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.fam \
    --memory 8000 \
    --clump ${OUT}/Plink_ldclump_aspirin_age_sex_pc3_studygxe_72269_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}.txt \
    --clump-p1 1 --clump-p2 1 \
    --out   ${OUT}/Plink_ldclump_aspirin_age_sex_pc3_studygxe_72269_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}

