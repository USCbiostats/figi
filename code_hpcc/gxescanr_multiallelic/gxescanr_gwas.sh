#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --account=lc_dvc
#SBATCH --partition=conti
#SBATCH --constraint=IB
#SBATCH --mail-user=andreeki@usc.edu
#SBATCH --mail-type=END
#SBATCH --array=1-22

Rscript gxescanr.R ${SLURM_ARRAY_TASK_ID} /staging/dvc/andreeki/gxescanr/${exposure}/FIGI_v2.3_gwasset_basic_covars_gxescan
