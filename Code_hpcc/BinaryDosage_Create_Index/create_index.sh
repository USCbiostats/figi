#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-9
#SBATCH --output=GetBDoseInfo_FIGI_ALL_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript GetBDoseInfo.R FIGI_ALL_chr${SLURM_ARRAY_TASK_ID}