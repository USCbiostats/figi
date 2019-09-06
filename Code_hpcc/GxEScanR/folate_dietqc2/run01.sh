#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-22
#SBATCH --output=./logs/GxEScanR_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447_binCovF_chr%a.log

Rscript job01.R ${SLURM_ARRAY_TASK_ID}
