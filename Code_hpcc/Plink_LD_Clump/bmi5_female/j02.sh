#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22
#SBATCH --output=../logs/Plink_bmi5_female_ldclump_chiSqGxE_chr%a.log

INN=/auto/pmd-01/andreeki/clump
INN=/staging/dvc/andreeki/clump
OUT=/staging/dvc/andreeki/clump_results/bmi5_female

plink \
    --bed ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.bed \
    --bim ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}_fix.bim \
    --fam ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.fam \
    --memory 8000 \
    --clump ${OUT}/Plink_bmi5_female_ldclump_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}.txt \
    --clump-p1 1 --clump-p2 1 \
    --out   ${OUT}/Plink_bmi5_female_ldclump_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}

