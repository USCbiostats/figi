#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-21
#SBATCH --output=oncoarray_clinicaltrials_bgzip_chr%a.log

cd $SLURM_SUBMIT_DIR
#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/corect_oncoarray_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --keep IDs_FIRE3_Mavericc_Tribe_vcftools.txt --recode --out oncoarray_fire3_mavericc_tribe_chr${SLURM_ARRAY_TASK_ID}
#oncoarray_fire3_mavericc_tribe_chr10.recode.vcf



bgzip -c oncoarray_fire3_mavericc_tribe_chr${SLURM_ARRAY_TASK_ID}.recode.vcf > oncoarray_fire3_mavericc_tribe_chr${SLURM_ARRAY_TASK_ID}.recode.vcf.gz
~