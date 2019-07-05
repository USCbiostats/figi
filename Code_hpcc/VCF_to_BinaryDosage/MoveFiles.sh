#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

wdir="/auto/pmd-02/figi/HRC_BDose"

cp ukbiobank_chr8.bdose ${wdir}
cp ukbiobank_chr9.bdose ${wdir}
cp ukbiobank_chr10.bdose ${wdir}
cp ukbiobank_chr11.bdose ${wdir}
cp ukbiobank_chr12.bdose ${wdir}
cp ukbiobank_chr13.bdose ${wdir}
cp ukbiobank_chr14.bdose ${wdir}
cp ukbiobank_chr15.bdose ${wdir}
cp ukbiobank_chr16.bdose ${wdir}
cp ukbiobank_chr17.bdose ${wdir}
cp ukbiobank_chr18.bdose ${wdir}
cp ukbiobank_chr19.bdose ${wdir}
cp ukbiobank_chr20.bdose ${wdir}
cp ukbiobank_chr21.bdose ${wdir}
cp ukbiobank_chr22.bdose ${wdir}

#cp FIGI_GECCO*.bdose ${wdir}
