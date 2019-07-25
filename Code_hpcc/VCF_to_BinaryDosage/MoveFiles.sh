#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --output=/auto/rcf-40/andreeki/FIGI_code/Code_hpcc/VCF_to_BinaryDosage/logs/MoveFiles_ukb1_chr%a.log

wdir="/auto/pmd-02/figi/HRC_BDose"

cp /staging/dvc/andreeki/BD/ukb2_chr1.bdose ${wdir}
cp /staging/dvc/andreeki/BD/ukb2_chr2.bdose ${wdir}
cp /staging/dvc/andreeki/BD/ukb2_chr3.bdose ${wdir}
cp /staging/dvc/andreeki/BD/ukb2_chr4.bdose ${wdir}
cp /staging/dvc/andreeki/BD/ukb2_chr5.bdose ${wdir}
cp /staging/dvc/andreeki/BD/ukb2_chr6.bdose ${wdir}
