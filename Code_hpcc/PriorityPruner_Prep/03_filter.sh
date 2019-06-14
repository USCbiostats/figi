#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

plink --bfile ./plink/${batch}_chr${chr} --extract ./plink/PriorityPruner_noUKB_ALL_snplist.txt  --make-bed --out ./plink/${batch}_chr${chr}_filter

#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.vcf.gz --positions PriorityPruner_noUKB_chr${chr}_snplist.txt --keep PriorityPruner_noUKB_10000_samplelist.txt --recode --out ${batch}_chr${chr}