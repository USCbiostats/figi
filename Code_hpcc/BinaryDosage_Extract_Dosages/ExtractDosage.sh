#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=2,3,6,7,10,16,17,18,21,22
#SBATCH --output=GxEScanR_GxE_aspref_FIGI_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript ExtractDose.R ${SLURM_ARRAY_TASK_ID}


