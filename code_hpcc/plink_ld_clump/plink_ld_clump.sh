#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22

#INN=/auto/pmd-01/andreeki/clump
INN=/staging/dvc/andreeki/clump
OUT=/staging/dvc/andreeki/clump_results/${exposure}

plink \
    --bed ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.bed \
    --bim ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.bim \
    --fam ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.fam \
    --memory 8000 \
    --clump ${OUT}/FIGI_v2.3_gxeset_${exposure}_basic_covars_gxescan_chr${SLURM_ARRAY_TASK_ID}_chiSqGxE_ldclump.txt \
    --clump-p1 1 --clump-p2 1 \
    --out   ${OUT}/FIGI_v2.3_gxeset_${exposure}_basic_covars_gxescan_chr${SLURM_ARRAY_TASK_ID}_chiSqGxE_ldclump

plink \
    --bed ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.bed \
    --bim ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.bim \
    --fam ${INN}/corect_oncoarray_controls_chr${SLURM_ARRAY_TASK_ID}.fam \
    --memory 8000 \
    --clump ${OUT}/FIGI_v2.3_gxeset_${exposure}_basic_covars_gxescan_chr${SLURM_ARRAY_TASK_ID}_chiSqG_ldclump.txt \
    --clump-p1 1 --clump-p2 1 \
    --out   ${OUT}/FIGI_v2.3_gxeset_${exposure}_basic_covars_gxescan_chr${SLURM_ARRAY_TASK_ID}_chiSqG_ldclump
