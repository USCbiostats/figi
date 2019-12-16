#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-22

# extract same markers from the KGP files (by chromosome)
# uses plink

REF=/auto/pmd-02/figi/andreeki/Refs/1KGP
OUT=/staging/dvc/andreeki/pca_ibd
cd ${OUT}

plink --bfile ${REF}/kgp.chr${SLURM_ARRAY_TASK_ID}.biallelic --extract ${OUT}/FIGI_PC_Backbone_Sample_30K_rsID.txt --memory 16000 --make-bed --out ${OUT}/tmp/kgp_backbone_chr${SLURM_ARRAY_TASK_ID}

