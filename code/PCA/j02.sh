#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB

# Merge chromosome plink files from j01.sh (by imputation batch)

#OUT=/staging/dvc/andreeki/pca_ibd
OUT=/auto/pmd-02/figi/PCA

# write out a list of files to merge to use with plink
for chr in {1..22}
do 
    echo ${OUT}/files/${batch}_backbone_chr${chr}
done > ${OUT}/files/mergelist_${batch}.txt

# call plink
plink --merge-list ${OUT}/files/mergelist_${batch}.txt --memory 8000 --make-bed --out ${OUT}/files/${batch}_backbone

#--------------------
# example submission
#sbatch --output=./logs/j02_ccfr_1m_1mduo_reimpute.log --export=batch='ccfr_1m_1mduo_reimpute' j02.sh
#sbatch --output=./logs/j02_corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray_nonEUR_reimpute' j02.sh
