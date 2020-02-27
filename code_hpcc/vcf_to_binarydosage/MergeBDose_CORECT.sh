#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=1-5
#SBATCH --output=./logs/MergeBDose_CORECT_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript MergeBDose_CORECT.R ${SLURM_ARRAY_TASK_ID}
