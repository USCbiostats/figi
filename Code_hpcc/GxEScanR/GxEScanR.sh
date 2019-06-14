#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-22
#SBATCH --output=./logs/GxEScanR_GxE_asp_ref_sex_age_pc10_studygxe_noUKB_58109_chr%a.log

Rscript GxEScanR_20190117.R ${SLURM_ARRAY_TASK_ID}