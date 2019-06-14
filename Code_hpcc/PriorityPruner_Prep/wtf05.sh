#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti


for chr in {1..22}
do
    plink --bfile ./plink/final_chr${chr} --recode transpose --out ./plink/final_transpose_chr${chr}
    Rscript tfamFix.R ./plink/final_transpose_chr${chr}
done


#vcftools --gzvcf /auto/pmd-02/figi/HRC_VCF_SampleRename/${batch}_chr${chr}.vcf.gz --positions PriorityPruner_noUKB_chr${chr}_snplist.txt --keep PriorityPruner_noUKB_10000_samplelist.txt --recode --out ${batch}_chr${chr}