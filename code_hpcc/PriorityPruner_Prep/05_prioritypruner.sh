#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

# just merge the files, the transpose to tped (per chromosome)

java -jar PriorityPruner.jar --tped ./plink/final_transpose_chr${chr}.tped --tfam ./plink/final_transpose_chr${chr}_sexfix.tfam --snp_table priorityPruner_pvals_chr${chr}.txt --r2 0.5 --chr ${chr} --out final_prioritypruner_pvalGxE_chr${chr}


#plink --merge-list ./plink/merge_list_chr${chr}.txt --keep-allele-order --make-bed --out ./plink/final_chr${chr}