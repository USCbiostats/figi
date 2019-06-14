#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

plink --vcf ${batch}_chr${chr}.recode.vcf --double-id --keep-allele-order --make-bed --out ./plink/${batch}_chr${chr}

cp ./plink/${batch}_chr${chr}.bim ./plink/${batch}_chr${chr}_backup.bim

Rscript fix_bim_file.R ./plink/${batch}_chr${chr}

#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.vcf.gz --positions PriorityPruner_noUKB_chr${chr}_snplist.txt --keep PriorityPruner_noUKB_10000_samplelist.txt --recode --out ${batch}_chr${chr}