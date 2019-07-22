#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-2
#SBATCH --output=../logs/Plink_folate_dietqc2_ldclump_chiSqGxE_chr%a.log

INN=/staging/dvc/andreeki/clump
OUT=/staging/dvc/andreeki/clump_results

plink --bed ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.bed --bim ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}_fix.bim --fam ${INN}/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.fam --memory 12000 --clump ${OUT}/Plink_clump_chiSqGxE_folate_dietqc2_chr${SLURM_ARRAY_TASK_ID}.txt --clump-p1 1 --clump-p2 1 --out ${OUT}/Plink_clump_chiSqGxE_folate_dietqc2_chr${SLURM_ARRAY_TASK_ID}

