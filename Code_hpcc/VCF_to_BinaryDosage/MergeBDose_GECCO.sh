#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=10,14
#SBATCH --output=MergeBDose_FIGI_GECCO_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript MergeBDose_GECCO.R ${SLURM_ARRAY_TASK_ID}