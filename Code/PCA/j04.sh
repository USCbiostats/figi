#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=./logs/j04.log

OUT=/staging/dvc/andreeki/pca_ibd
cd ${OUT}

# create 2 sets:
    # gwas + KGP
    # gxe + KGP

# gwas set + kgp
plink --bfile FIGI_GwasSet --memory 16000 --bmerge kgp_backbone --make-bed --out FIGI_GwasSet_KGP
# gxe set + kgp
plink --bfile FIGI_GxESet --memory 16000 --bmerge kgp_backbone --make-bed --out FIGI_GxESet_KGP

