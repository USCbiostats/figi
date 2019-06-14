#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

OUT=/staging/dvc/andreeki/posthoc

cd ${OUT}

#tabix -p vcf ${REF}/cytosnp_comb_chr22.vcf.gz

echo "Start time is ---- $(date)"

epacts single --vcf ${OUT}/FIGI_ALL_chr22_topSNPs.vcf.gz \
              --ped ~/FIGI_Covariates_08272018_EPACTS.ped \
              --min-maf 0.01 --chr 22 --pheno DISEASE --cov age_ref --cov SEX --cov PC1 --cov PC2 --cov PC3 --cov PC4 --test b.wald --field DS --out FIGI_ALL_chr22_topSNPs.out --run 8

echo "End time is ---- $(date)"
