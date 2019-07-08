#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=5
#SBATCH --job-name=GxEScanR_sex
#SBATCH --output=GxEScanR_sex_ukb_original_test_chr%a.log

Rscript ukb_step01.R ${SLURM_ARRAY_TASK_ID}