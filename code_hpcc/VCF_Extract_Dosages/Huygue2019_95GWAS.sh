#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-22

vcftools --gzvcf /auto/pmd-02/figi/HRC_ImpSrv_Output/ccfr_oncoarray_postqc/extract/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz --positions GWAS95_vcftools.txt --recode --out ccfr_oncoarray_GWAS95_chr${SLURM_ARRAY_TASK_ID}