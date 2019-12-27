#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1,2,9,10,13
#SBATCH --output=/auto/rcf-40/andreeki/FIGI_code/Code_hpcc/Plink_LD_Clump/logs/job01_chr%a.log

# plink has a function to clump GWAS hits (similar to priority pruner)
# will once again use corect_oncoarray to calculate LD --- check with the PIs on this
# alternatively, can merge plink files across imputation batches, will just take time and disk space.

# for this script, create correct oncoarray plink files. 
# it's outdated, since I did this for controls only, files are stored in the pmd drive

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
OUT=/staging/dvc/andreeki/clump
cd ${OUT}

plink --vcf ${REF}/${batch}_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --double-id --snps-only --biallelic-only --keep-allele-order --memory 8000 --make-bed --out ${OUT}/${batch}_chr${SLURM_ARRAY_TASK_ID}

#sbatch --output=corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray' j01.sh
