#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=18,20,21
#SBATCH --output=/auto/rcf-40/andreeki/FIGI_code/Code_hpcc/Get_Rsq_Estimate/logs/GetRsqEstimateMinimac_chr%a.log

Rscript job01.R ${SLURM_ARRAY_TASK_ID}
