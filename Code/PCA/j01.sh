#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-22

# Extract snp backbone (random sample) from VCFs (all batches + chromosomes)
# ${batch} should be imputation batch name, export during sbatch command
# also export log file (slurm doesn't support exported variables from being used in #SBATCH lines..)

REF=/auto/pmd-02/figi/HRC_VCF_SampleRename
#OUT=/staging/dvc/andreeki/pca_ibd
OUT=/auto/pmd-02/figi/PCA
cd ${OUT}

vcftools --gzvcf ${REF}/${batch}_chr${SLURM_ARRAY_TASK_ID}.vcf.gz --positions ${OUT}/FIGI_PC_Backbone_Sample_30K.txt --recode --out ${OUT}/files/${batch}_backbone_chr${SLURM_ARRAY_TASK_ID}

plink --vcf ${OUT}/files/${batch}_backbone_chr${SLURM_ARRAY_TASK_ID}.recode.vcf --double-id --snps-only --biallelic-only --keep-allele-order --memory 8000 --make-bed --out ${OUT}/files/${batch}_backbone_chr${SLURM_ARRAY_TASK_ID}


# example submission
#sbatch --output=ccfr_1m_1mduo_reimpute_chr${SLURM_ARRAY_TASK_ID}.log --export=batch='ccfr_1m_1mduo_reimpute' j01.sh

#sbatch --output=corect_oncoarray_nonEUR_reimpute.log --export=batch='corect_oncoarray_nonEUR_reimpute' j01.sh
