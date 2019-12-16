#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=9,11-15
#SBATCH --output=Plink_aspirin_ldclump_chiSqGxE_chr%a.log

# run genome wide clumping (just to see)

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/clump
cd ${OUT}

#plink --bed corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.bed --bim corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}_fix.bim --fam corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.fam --clump Plink_aspirin_ldclump_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}.txt --clump-p1 1 --clump-p2 1 --out clump_aspirin_chiSqGxE_chr${SLURM_ARRAY_TASK_ID}

plink --bed corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.bed --bim corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}_fix.bim --fam corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.fam --memory 12000 --clump Plink_aspirin_ldclump_chiSqG_chr${SLURM_ARRAY_TASK_ID}.txt --clump-p1 0.001 --clump-p2 0.1   --out clump_aspirin_chiSqG_chr${SLURM_ARRAY_TASK_ID}

#sbatch --output=corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray' j01.sh