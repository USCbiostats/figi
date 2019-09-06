#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --array=4-5
#SBATCH --output=/staging/dvc/andreeki/GxE/logs/GxEScanR_asp_ref_sex_age_pc3_studygxe_72820_binCovF_chr%a.log

Rscript job01.R ${SLURM_ARRAY_TASK_ID} /staging/dvc/andreeki/GxE/asp_ref/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820
