#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --array=4-5
#SBATCH --output=./logs/MergeBDose_FIGI_chr%a.log

cd $SLURM_SUBMIT_DIR
Rscript MergeBDose.R ${SLURM_ARRAY_TASK_ID}

