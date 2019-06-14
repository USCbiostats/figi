#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --job-name=j05
#SBATCH --output=j05.log

REF=/auto/pmd-02/figi/andreeki/Refs/1KGP
OUT=/staging/dvc/andreeki/PCA

cd ${OUT}

# merge kgp chromosomes 1-22
for j in {1..22}; do echo ${OUT}/tmp/kgp_backbone_chr$j; done > ${OUT}/mergelist_kgp.txt
plink --merge-list ${OUT}/mergelist_kgp.txt --make-bed --out ${OUT}/kgp_backbone

# (Note - Grace instructions had pruning step by 1KGP super-population)
# (not doing that here)

# for PCA calculation - convert HRC SNP IDs into rsIDs to merge with KGP SNPs
awk '{print $1, $14}' hrc_kgp_innerjoin_rsid.txt > update_ids.txt # recall this file is created in the SNP_Sample_30K.R file

plink --bfile ALL_Merged_PCA --memory 16000 --update-name update_ids.txt --make-bed --out TMP2
plink --bfile TMP2 --memory 16000 --bmerge kgp_backbone --make-bed --out ALL_Merged_PCA_KGP


plink --bfile ALL_Merged_PCA --memory 16000 --keep ~/FIGI_Sample_Rename_VCF_08162018_ALL_drops.txt --make-bed --out ALL_Merged_PCA_drops