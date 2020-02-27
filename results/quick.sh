#!/bin/bash

# ---- asp_ref ---- #
#Rscript preprocess_gxescanr_results.R asp_ref FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan
Rscript create_qq_manhattan_plots.R asp_ref FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3
#Rscript create_two_step_plots.R asp_ref FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan age_ref_imp sex study_gxe PC1 PC2 PC3

# ---- alcoholc_moderate ---- #
#Rscript preprocess_gxescanr_results_multiallelic.R alcoholc_moderate FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan 
#Rscript create_qq_manhattan.R alcoholc_moderate
#Rscript create_two_step_plots.R alcoholc_moderate FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe PC1 PC2 PC3

# ---- alcoholc_heavy ---- #
#Rscript preprocess_gxescanr_results_multiallelic.R alcoholc_heavy FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan
#Rscript create_two_step_plots.R alcoholc_heavy FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe PC1 PC2 PC3

#Rscript preprocess_gxescanr_results_multiallelic.R hrt_ref_pm FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_gxescan 
#Rscript preprocess_gxescanr_results.R hrt_ref_pm2 FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_gxescan

#Rscript preprocess_gxescanr_results_multiallelic.R eo_ref_pm_gxe FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_gxescan
#Rscript preprocess_gxescanr_results_multiallelic.R ep_ref_pm_gxe FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_gxescan
