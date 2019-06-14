#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=5
#SBATCH --output=GetBDoseInfo_FIGI_ALL_chr%a.log

#Rscript GetBDoseInfo.R /auto/pmd-02/figi/HRC_BDose/FIGI_ALL_chr${SLURM_ARRAY_TASK_ID}.bdose
cd $SLURM_SUBMIT_DIR
#Rscript GetBDoseInfo.R /auto/pmd-02/figi/HRC_BDose/corsa_axiom_chr${SLURM_ARRAY_TASK_ID}.bdose
Rscript GetBDoseInfo.R /staging/dvc/andreeki/bdose/ukbiobank_chr${SLURM_ARRAY_TASK_ID}