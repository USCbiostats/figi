#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti

# just merge the files, the transpose to tped (per chromosome)

plink --merge-list ./plink/merge_list_chr${chr}.txt --keep-allele-order --make-bed --out ./plink/final_chr${chr}

# transpose using plink (tped format)
# (insert example command here)