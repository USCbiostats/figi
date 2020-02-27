#!/bin/bash

#sbatch --job-name=gxescanr_asp_ref --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan.log --export=exposure=asp_ref gxescanr.sh 


#sbatch --job-name=gxe_hrt_ref_pm --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_gxescan.log --export=exposure=hrt_ref_pm gxescanr.sh

#sbatch --job-name=gxe_asp_ref --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan.log --export=exposure=asp_ref gxescanr.sh

#sbatch --job-name=gxe_aspirin --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan.log --export=exposure=aspirin gxescanr.sh
#sbatch --job-name=gxe_nsaids --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_nsaids_basic_covars_gxescan.log --export=exposure=nsaids gxescanr.sh

#sbatch --job-name=gwas --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gwasset_basic_covars_gxescan.log --export=exposure=gwas gxescanr_gwas.sh

#sbatch --job-name=gxe_hrt_ref_pm2 --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_gxescan.log --export=exposure=hrt_ref_pm2 gxescanr.sh

#sbatch --job-name=gxe_eo_ref_pm_gxe --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_gxescan.log --export=exposure=eo_ref_pm_gxe gxescanr.sh

#sbatch --job-name=gxe_ep_ref_pm_gxe --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_gxescan.log --export=exposure=ep_ref_pm_gxe gxescanr.sh



# Alcohol
#sbatch --job-name=gxe_alcohol_mod --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan.log --export=exposure=alcoholc_moderate gxescanr.sh
#sbatch --job-name=gxe_alcohol_heavy --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan.log --export=exposure=alcoholc_heavy gxescanr.sh
#sbatch --job-name=gxe_alcohol_heavy_vs_moderate --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_gxescan.log --export=exposure=alcoholc_heavy_vs_moderate gxescanr.sh

# BMI
#sbatch --job-name=gxe_bmi5 --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_bmi5_basic_covars_gxescan.log --export=exposure=bmi5 gxescanr.sh

# height10
sbatch --job-name=gxe_height10 --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_height10_basic_covars_gxescan.log --export=exposure=height10 gxescanr.sh

# T2D
#sbatch --job-name=gxe_diab --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_diab_basic_covars_gxescan.log --export=exposure=diab gxescanr.sh

# Meat Intake
#sbatch --job-name=gxe_redmeatqc2 --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_gxescan.log --export=exposure=redmeatqc2 gxescanr.sh
#sbatch --job-name=gxe_procmeatqc2 --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_procmeatqc2_basic_covars_gxescan.log --export=exposure=procmeatqc2 gxescanr.sh

# Folate
#sbatch --job-name=gxe_foltot --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_folate_totqc2_basic_covars_gxescan.log --export=exposure=folate_totqc2 gxescanr.sh
#sbatch --job-name=gxe_foldiet --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_folate_dietqc2_basic_covars_gxescan.log --export=exposure=folate_dietqc2 gxescanr.sh

# Calcium 
#sbatch --job-name=gxe_catot  --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_calcium_totqc2_basic_covars_gxescan.log --export=exposure=calcium_totqc2 gxescanr.sh
#sbatch --job-name=gxe_cadiet --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_calcium_dietqc2_basic_covars_gxescan.log --export=exposure=calcium_dietqc2 gxescanr.sh

# Smoking
#sbatch --job-name=gxe_smk_ever  --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_smk_ever_basic_covars_gxescan.log --export=exposure=smk_ever gxescanr.sh
#sbatch --job-name=gxe_smk_aveday --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_smk_aveday_basic_covars_gxescan.log --export=exposure=smk_aveday gxescanr.sh
#sbatch --job-name=gxe_smk_pkyr  --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_smk_pkyr_basic_covars_gxescan.log --export=exposure=smk_pkyr gxescanr.sh


# fruit vegetable fiber
#sbatch --job-name=gxe_fiber  --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_fiberqc2_basic_covars_gxescan.log --export=exposure=fiberqc2 gxescanr.sh
#sbatch --job-name=gxe_vegetable --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_vegetableqc2_basic_covars_gxescan.log --export=exposure=vegetableqc2 gxescanr.sh
#sbatch --job-name=gxe_fruit  --output=/staging/dvc/andreeki/gxescanr/logs/FIGI_v2.3_gxeset_fruitqc2_basic_covars_gxescan.log --export=exposure=fruitqc2 gxescanr.sh
