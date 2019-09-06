#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=./logs/j04.log

OUT=/auto/pmd-02/figi/PCA
cd ${OUT}

# create 2 sets:
    # gwas + KGP
    # gxe + KGP

# gwas set + kgp
plink --bfile FIGI_GwasSet_190729 --memory 16000 --bmerge ${OUT}/files/kgp_backbone --make-bed --out FIGI_GwasSet_KGP_190729
# gxe set + kgp
plink --bfile FIGI_GxESet_190729 --memory 16000 --bmerge ${OUT}/files/kgp_backbone --make-bed --out FIGI_GxESet_KGP_190729

