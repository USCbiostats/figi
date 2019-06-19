#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-9

# plink has a function to clump GWAS hits (similar to priority pruner)
# will once again use corect_oncoarray to calculate LD --- check with the PIs on this
# alternatively, can merge plink files across imputation batches, will just take time and disk space.

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/clump
cd ${OUT}

plink --vcf ${REF}/${batch}_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --double-id --snps-only --biallelic-only --keep-allele-order --make-bed --out ${OUT}/${batch}_chr${SLURM_ARRAY_TASK_ID}

#sbatch --output=corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray' j01.sh