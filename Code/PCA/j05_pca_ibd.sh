#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=./logs/j05_PCA_IBD.log

plink2 --bfile FIGI_GwasSet     --pca 20 approx --out FIGI_GwasSet
plink2 --bfile FIGI_GwasSet_KGP --pca 20 approx --out FIGI_GwasSet_KGP
plink2 --bfile FIGI_GxESet      --pca 20 approx --out FIGI_GxESet
plink2 --bfile FIGI_GxESet_KGP  --pca 20 approx --out FIGI_GxESet_KGP

king -b FIGI_GwasSet.bed --cpus 8 --related --degree 2 --prefix FIGI_GwasSet
king -b FIGI_GxESet.bed --cpus 8 --related --degree 2 --prefix FIGI_GxESet