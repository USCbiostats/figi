#!/bin/bash

#sbatch --job-name=gxescanr_asp_ref --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan.log --export=exposure=asp_ref gxescanr.sh 

#sbatch --job-name=gxe_alcohol_mod --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan.log --export=exposure=alcoholc_moderate gxescanr.sh
#sbatch --job-name=gxe_alcohol_heavy --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan.log --export=exposure=alcoholc_heavy gxescanr.sh

sbatch --job-name=gxe_hrt_ref_pm --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_gxescan.log --export=exposure=hrt_ref_pm gxescanr.sh

#sbatch --job-name=gxe_asp_ref --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan.log --export=exposure=asp_ref gxescanr.sh

#sbatch --job-name=gwas --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gwasset_basic_covars_gxescan.log --export=exposure=gwas gxescanr_gwas.sh


