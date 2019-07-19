#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=22
#SBATCH --output=./logs/MergeBDose_FIGI_CORECT_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript MergeBDose_CORECT.R ${SLURM_ARRAY_TASK_ID}
