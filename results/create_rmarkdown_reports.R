#=============================================================================#
# create post-harmonization and results
#=============================================================================#


# ---- BMI ---- 
# bmi5
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'bmi5',
  is_exposure_categorical = F,
  energy_adj = F,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/bmi5_results.html")




# ---- Alcohol ---- 
# alcoholc_heavy
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'alcoholc_heavy',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/alcoholc_heavy_results.html")

# alcoholc_moderate
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'alcoholc_moderate',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/alcoholc_moderate_results.html")

# alcoholc_heavy_vs_moderate
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'alcoholc_heavy_vs_moderate',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/alcoholc_heavy_vs_moderate_results.html")



# ---- NSAIDS ---- 
# asp_ref
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'asp_ref',
  is_exposure_categorical = T,
  energy_adj = F,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/asp_ref_results.html")

# aspirin
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'aspirin',
  is_exposure_categorical = T,
  energy_adj = F,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/aspirin_results.html")

# nsaids
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'nsaids',
  is_exposure_categorical = T,
  energy_adj = F,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/nsaids_results.html")





# ---- Smoking ---- 
# smk_ever
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'smk_ever',
  is_exposure_categorical = T,
  energy_adj = F,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/smk_ever_results.html")


# smk_aveday
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'smk_aveday',
  is_exposure_categorical = F,
  energy_adj = F,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/smk_aveday_results.html")


# smk_pkyr
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'smk_pkyr',
  is_exposure_categorical = F,
  energy_adj = F,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/smk_pkyr_results.html")



# ---- Red and Processed meat ----
# redmeatqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'redmeatqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/redmeatqc2_results.html")


# procmeatqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'procmeatqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/procmeatqc2_results.html")




# ---- Fruits, vegetables, fiber ----

# fruitqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'fruitqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/fruitqc2_results.html")


# fiberqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'fiberqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/fiberqc2_results.html")


# vegetableqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'vegetableqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/vegetableqc2_results.html")



# ---- HRT ---- 
# NO SEX STRATIFIED
# hrt_ref_pm
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'hrt_ref_pm',
  is_exposure_categorical = T,
  energy_adj = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/hrt_ref_pm_results.html")


# hrt_ref_pm2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'hrt_ref_pm2',
  is_exposure_categorical = T,
  energy_adj = F,
  table1_by_e = T, 
  additional_analyses = F,
  covariates_suffix = 'age_ref_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/hrt_ref_pm2_results.html")


# eo_ref_pm_gxe
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'eo_ref_pm_gxe',
  is_exposure_categorical = T,
  energy_adj = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/eo_ref_pm_gxe_results.html")

# ep_ref_pm_gxe
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'ep_ref_pm_gxe',
  is_exposure_categorical = T,
  energy_adj = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/ep_ref_pm_gxe_results.html")


# ---- T2D ---- #
# diab
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'diab',
  is_exposure_categorical = T,
  energy_adj = F,
  covariates_suffix = 'age_ref_imp_sex_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/diab_results.html")


# ---- Folate ----
# folate_totqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'folate_totqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  additional_analyses = T,
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/folate_totqc2_results.html")


# folate_dietqc2
rm(list = ls())
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'folate_dietqc2',
  is_exposure_categorical = T,
  energy_adj = T,
  table1_by_e = F, 
  covariates_suffix = 'age_ref_imp_sex_energytot_imp_study_gxe_pc1_pc2_pc3'
), output_file = "~/Dropbox/FIGI/Results/folate_dietqc2_results.html")
