#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

# Merge chromosome plink files from j01.sh (by imputation batch)

OUT=/staging/dvc/andreeki/pca_ibd
cd ${OUT}

# write out a list of files to merge to use with plink
for chr in {1..22}
do 
    echo ${OUT}/tmp/${batch}_backbone_chr${chr}
done > ${OUT}/tmp/mergelist_${batch}.txt

# call plink
plink --merge-list ${OUT}/tmp/mergelist_${batch}.txt --make-bed --out ${OUT}/${batch}_backbone

#--------------------
# example submission
#sbatch --output=./logs/j02_ccfr_1m_1mduo_reimpute.log --export=batch='ccfr_1m_1mduo_reimpute' j02.sh
#sbatch --output=./logs/j02_corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray_nonEUR_reimpute' j02.sh