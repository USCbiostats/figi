#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-22
#SBATCH --output=BinaryDosage_corect_oncoarray_nonEUR_reimpute_chr%a.log

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/bdose
OUTF=/auto/pmd-02/figi/HRC_BDose

cd ${OUT}

time zcat ${REF}/corect_oncoarray_nonEUR_reimpute_chr${SLURM_ARRAY_TASK_ID}.vcf.gz > ${OUT}/corect_oncoarray_nonEUR_reimpute_chr${SLURM_ARRAY_TASK_ID}.vcf

time Rscript -e 'library(BinaryDosage); args <- commandArgs(trailingOnly=T); VCFtoBD(paste0(args[1], "'".vcf"'"), paste0(args[1], "'".bdose"'"))' corect_oncoarray_nonEUR_reimpute_chr${SLURM_ARRAY_TASK_ID}

#rm ${OUT}/$1_chr$2.vcf
#mv ${OUT}/$1_chr$2.bdose ${OUTF}