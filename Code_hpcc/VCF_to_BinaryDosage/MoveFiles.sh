#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

wdir="/auto/pmd-02/figi/HRC_BDose"

#cp ukbiobank_chr16.bdose ${wdir}
#cp ukbiobank_chr17.bdose ${wdir}
#cp ukbiobank_chr18.bdose ${wdir}
#cp ukbiobank_chr19.bdose ${wdir}
#cp ukbiobank_chr20.bdose ${wdir}
#cp ukbiobank_chr21.bdose ${wdir}
#cp ukbiobank_chr22.bdose ${wdir}

cp corect_oncoarray_nonEUR_reimpute_chr*.bdose ${wdir}