#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --array=1-3
#SBATCH --output=./logs/GetBDoseInfo_FIGI_chr%a.log

#cd $SLURM_SUBMIT_DIR

Rscript GetBDoseInfo.R /auto/pmd-02/figi/HRC_BDose/FIGI_chr${SLURM_ARRAY_TASK_ID}
#Rscript GetBDoseInfo.R /auto/pmd-02/figi/HRC_BDose/corsa_axiom_chr${SLURM_ARRAY_TASK_ID}.bdose
#Rscript GetBDoseInfo.R /staging/dvc/andreeki/bdose/ukbiobank_chr${SLURM_ARRAY_TASK_ID}

