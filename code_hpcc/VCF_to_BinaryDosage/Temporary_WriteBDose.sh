#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --array=22
#SBATCH --output=/auto/rcf-40/andreeki/FIGI_code/Code_hpcc/VCF_to_BinaryDosage/logs/Temporary_WriteBD_UKB1_chr%a.log

time Rscript Temporary_WriteBDose.R ${SLURM_ARRAY_TASK_ID}
