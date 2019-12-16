#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --array=1-22
#SBATCH --output=/staging/dvc/andreeki/GxE/logs/GxEScanR_nsaids_sex_age_pc3_studygxe_70335_binCovF_chr%a.log

Rscript job01.R ${SLURM_ARRAY_TASK_ID} /staging/dvc/andreeki/GxE/nsaids/FIGI_GxESet_nsaids_sex_age_pc3_studygxe_70335
