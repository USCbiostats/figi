#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=2

Rscript test.R ${SLURM_ARRAY_TASK_ID} 
