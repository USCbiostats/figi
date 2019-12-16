#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22
#SBATCH --output=../logs/Plink_asp_ref_ldclump_chiSqGxE_controls_chr%a.log

INN=/auto/pmd-01/andreeki/clump
INN=/staging/dvc/andreeki/clump
OUT=/staging/dvc/andreeki/clump_results/asp_ref

plink \
    --bfile ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID} \
    --memory 8000 \
    --clump ${OUT}/Plink_asp_ref_ldclump_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}.txt \
    --clump-p1 1 --clump-p2 1 \
    --out   ${OUT}/Plink_asp_ref_ldclump_chiSqGxE_controls_chr${SLURM_ARRAY_TASK_ID}

