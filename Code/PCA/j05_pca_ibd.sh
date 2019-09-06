#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=./logs/j05_PCA_IBD.log

OUT=/auto/pmd-02/figi/PCA

plink2 --bfile ${OUT}/FIGI_GwasSet_190729     --pca 20 approx --out ${OUT}/FIGI_GwasSet_190729
plink2 --bfile ${OUT}/FIGI_GwasSet_KGP_190729 --pca 20 approx --out ${OUT}/FIGI_GwasSet_KGP_190729
plink2 --bfile ${OUT}/FIGI_GxESet_190729      --pca 20 approx --out ${OUT}/FIGI_GxESet_190729
plink2 --bfile ${OUT}/FIGI_GxESet_KGP_190729  --pca 20 approx --out ${OUT}/FIGI_GxESet_KGP_190729

king -b ${OUT}/FIGI_GwasSet_190729.bed --cpus 8 --related --degree 2 --prefix ${OUT}/FIGI_GwasSet_190729
king -b ${OUT}/FIGI_GxESet_190729.bed --cpus 8 --related --degree 2 --prefix ${OUT}/FIGI_GxESet_190729
