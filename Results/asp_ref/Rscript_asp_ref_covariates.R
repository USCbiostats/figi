#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
# FIGI Analysis 08/08/2019 - reimputed UKB + reach
# FIGI Analysis 12/12/2019
#
# Create covariate tables for GxEScanR
# 
# ------ Notes ------ #
# use gxe set as determined in `gxe` variable (in addition to drop == 0)
# 	- gxe == 1 removes nonEUR studies
#		- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
# study_gxe = adj. for study + platform
# remove case-only or control-only studies (causes gxescanr to crash)
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
rm(list = ls())
figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")


cov <- format_data_glm(figi_gwas, 'asp_ref', T, 0, c("energytot_imp"))
cov_gxescan <- format_data_gxescan(cov, 'asp_ref')



#-----------------------------------------------------------------------------#
# OLD CODE
#-----------------------------------------------------------------------------#
# cov_gxe <- figi %>%
#   filter(drop == 0 & gxe == 1) %>% 
#   inner_join(pc, by = c('vcfid' = 'IID'))
# 
# # checks
# xtabs(~ cov_gxe$asp_ref, addNA = T)
# xtabs(~ cov_gxe$sex, addNA = T)
# xtabs(~ cov_gxe$study_gxe, addNA = T)
# xtabs(~ cov_gxe$outc + cov_gxe$asp_ref, addNA = T)
# 
# 
# # remove missing asp_ref, recode covariates
# cov_full <- cov_gxe %>% 
#   mutate(outcome = ifelse(outc == "Case", 1, 0),
#          age_ref_imp = as.numeric(age_ref_imp),
#          sex = ifelse(sex == "Female", 0, 1),
#          study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
#          asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) 
# 
# cov <- cov_full %>% 
#   filter(asp_ref != "") %>% 
#   dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3)), asp_ref) %>% 
#   filter(complete.cases(.))
# 
# table(cov$study_gxe, cov$outcome)
# sort(unique(cov$study_gxe))
# drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
#   filter(Freq == 0)
# 
# cov <- filter(cov, !study_gxe %in% unique(drops$Var1))
# 
# cov_full_subset <- cov_full %>% 
#   filter(vcfid %in% cov$vcfid)
# 
# 
# 
# 
# 
# # create indicator variables for GxEScanR
# cov_gxescan <- cov
# for(t in unique(cov_gxescan$study_gxe)) {
#   cov_gxescan[paste0(t)] <- ifelse(cov_gxescan$study_gxe==t,1,0)
# }
# 
# # checks (what studies are included and how many)
# check <- dplyr::select(cov_gxescan, unique(cov_gxescan$study_gxe)) %>% 
#   summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
#   rownames_to_column()
# 
# # clean up , make Kentucky reference study_gxe (arbitrary)
# cov_gxescan <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -asp_ref, asp_ref)
# 
# # save files
# # saveRDS(cov,         file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72269_GLM.rds', version = 2)
# # saveRDS(cov_gxescan, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72269.rds",     version = 2)
# 
# # save files
# # remember that you only recoded outcome, asp_ref, sex, etc. So yes it contains all covariates but it's stll specific to the exposure in question. I'll create individual files now but eventually can create a single global, cleaned E dataset
# saveRDS(cov_full_subset, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset_asp_ref_age_ref_imp_sex_study_gxe_pc1_pc2_pc3_full.rds', version = 2)
# saveRDS(cov,             file = '~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset_asp_ref_age_ref_imp_sex_study_gxe_pc1_pc2_pc3_glm.rds', version = 2)
# saveRDS(cov_gxescan,     file = "~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset_asp_ref_age_ref_imp_sex_study_gxe_pc1_pc2_pc3.rds",      version = 2)

