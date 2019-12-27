#!/bin/bash

# ---- asp_ref ---- #
Rscript preprocess_gxescanr_results.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820
Rscript create_qq_manhattan_plots.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820 age_ref_imp sex study_gxe PC1 PC2 PC3
Rscript create_qq_manhattan_plots_ld_clumped.R asp_ref FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820 age_ref_imp sex study_gxe PC1 PC2 PC3

# ---- aspirin ---- #
Rscript preprocess_gxescanr_results.R aspirin FIGI_GxESet_aspirin_age_sex_pc3_studygxe_72269 
Rscript create_qq_manhattan_plots.R aspirin FIGI_GxESet_aspirin_age_sex_pc3_studygxe_72269 age_ref_imp sex study_gxe PC1 PC2 PC3
Rscript create_qq_manhattan_plots_ld_clumped.R aspirin FIGI_GxESet_aspirin_age_sex_pc3_studygxe_72269  age_ref_imp sex study_gxe PC1 PC2 PC3

# ---- nsaids ---- #
Rscript preprocess_gxescanr_results.R nsaids FIGI_GxESet_nsaids_age_sex_pc3_studygxe_70335
Rscript create_qq_manhattan_plots.R nsaids FIGI_GxESet_nsaids_age_sex_pc3_studygxe_70335 age_ref_imp sex study_gxe PC1 PC2 PC3
Rscript create_qq_manhattan_plots_ld_clumped.R nsaids FIGI_GxESet_nsaids_age_sex_pc3_studygxe_70335  age_ref_imp sex study_gxe PC1 PC2 PC3

# ---- gwas ---- #


# ---- locuszoom plots ---- #
