#!/bin/bash

#sbatch --job-name=clump --output=/staging/dvc/andreeki/clump_results/logs/plink_asp_ref_ldclump_chiSqGxE_chr%a.log --export=exposure=asp_ref plink_ld_clump.sh
#sbatch --job-name=clump_alc_mod --output=/staging/dvc/andreeki/clump_results/logs/plink_alcoholc_moderate_ldclump_chiSqGxE_chr%a.log --export=exposure=alcoholc_moderate plink_ld_clump.sh
#sbatch --job-name=clump_alc_heavy --output=/staging/dvc/andreeki/clump_results/logs/plink_alcoholc_heavy_ldclump_chiSqGxE_chr%a.log --export=exposure=alcoholc_heavy plink_ld_clump.sh


# ---- NSAIDs ---- #
#sbatch --job-name=ld_asp_ref --output=/staging/dvc/andreeki/clump_results/logs/plink_asp_ref_ldclump_chiSqGxE_chr%a.log --export=exposure=asp_ref plink_ld_clump.sh
#sbatch --job-name=ld_nsaids --output=/staging/dvc/andreeki/clump_results/logs/plink_nsaids_ldclump_chr%a.log --export=exposure=nsaids plink_ld_clump.sh
#sbatch --job-name=ld_aspirin --output=/staging/dvc/andreeki/clump_results/logs/plink_aspirin_ldclump_chr%a.log --export=exposure=aspirin plink_ld_clump.sh

# ---- alcoholc ---- #
#sbatch --job-name=clump_alc_h_vs_m --output=/staging/dvc/andreeki/clump_results/logs/plink_alcoholc_heavy_vs_moderate_ldclump_chr%a.log --export=exposure=alcoholc_heavy_vs_moderate plink_ld_clump.sh
#sbatch --job-name=clump_alc_m --output=/staging/dvc/andreeki/clump_results/logs/plink_alcoholc_moderate_ldclump_chr%a.log --export=exposure=alcoholc_moderate plink_ld_clump.sh
#sbatch --job-name=clump_alc_h --output=/staging/dvc/andreeki/clump_results/logs/plink_alcoholc_heavy_ldclump_chr%a.log --export=exposure=alcoholc_heavy plink_ld_clump.sh

# ---- red/proc meat ----#
#sbatch --job-name=clump_red --output=/staging/dvc/andreeki/clump_results/logs/plink_redmeatqc2_ldclump_chr%a.log --export=exposure=redmeatqc2 plink_ld_clump.sh
#sbatch --job-name=clump_proc --output=/staging/dvc/andreeki/clump_results/logs/plink_procmeatqc2_ldclump_chr%a.log --export=exposure=procmeatqc2 plink_ld_clump.sh
#sbatch --job-name=clump_fiber --output=/staging/dvc/andreeki/clump_results/logs/plink_fiberqc2_ldclump_chr%a.log --export=exposure=fiberqc2 plink_ld_clump.sh
#sbatch --job-name=clump_fruit --output=/staging/dvc/andreeki/clump_results/logs/plink_fruitqc2_ldclump_chr%a.log --export=exposure=fruitqc2 plink_ld_clump.sh
#sbatch --job-name=clump_veg --output=/staging/dvc/andreeki/clump_results/logs/plink_vegetableqc2_ldclump_chr%a.log --export=exposure=vegetableqc2 plink_ld_clump.sh


# ---- bmi5 ---- #
#sbatch --job-name=clump_bmi5 --output=/staging/dvc/andreeki/clump_results/logs/plink_bmi5_ldclump_chr%a.log --export=exposure=bmi5 plink_ld_clump.sh

# ---- diab ---- #
#sbatch --job-name=clump_diab --output=/staging/dvc/andreeki/clump_results/logs/plink_diab_ldclump_chr%a.log --export=exposure=diab plink_ld_clump.sh

# ---- hrt ---- #
sbatch --job-name=clump_hrt --output=/staging/dvc/andreeki/clump_results/logs/plink_hrt_ref_pm2_ldclump_chr%a.log --export=exposure=hrt_ref_pm2 plink_ld_clump.sh
sbatch --job-name=clump_eo --output=/staging/dvc/andreeki/clump_results/logs/plink_eo_ref_pm_gxe_ldclump_chr%a.log --export=exposure=eo_ref_pm_gxe  plink_ld_clump.sh
sbatch --job-name=clump_ep --output=/staging/dvc/andreeki/clump_results/logs/plink_ep_ref_pm_gxe_ldclump_chr%a.log --export=exposure=ep_ref_pm_gxe plink_ld_clump.sh
