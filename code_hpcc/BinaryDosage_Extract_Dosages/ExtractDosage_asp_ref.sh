#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-8,10-12,14,15,17-20
#SBATCH --output=./logs/ExtractDose_asp_ref_FIGI_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript ExtractDose_asp_ref.R ${SLURM_ARRAY_TASK_ID}