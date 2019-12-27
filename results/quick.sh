#!/bin/bash

# ---- asp_ref ---- #
#Rscript preprocess_gxescanr_results.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820
#Rscript create_qq_manhattan_plots.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820 age_ref_imp sex study_gxe PC1 PC2 PC3
#Rscript create_qq_manhattan_plots_ld_clumped.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820 age_ref_imp sex study_gxe PC1 PC2 PC3

# ---- alcoholc_moderate ---- #
#Rscript preprocess_gxescanr_results.R alcoholc_moderate FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan 
#Rscript create_qq_manhattan_plots.R alcoholc_moderate FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe PC1 PC2 PC3

# ---- alcoholc_heavy ---- #
#Rscript preprocess_gxescanr_results.R alcoholc_heavy FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan
#Rscript create_qq_manhattan_plots.R alcoholc_heavy FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe PC1 PC2 PC3

Rscript preprocess_gxescanr_results.R asp_ref FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan 
