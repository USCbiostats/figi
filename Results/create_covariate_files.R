#=============================================================================#
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


#-----------------------------------------------------------------------------#
# create exposure phenotype files for gxescan
#
# ------ Notes ------ #
# I check counts against old dataset. Gwas and Gxe match in numbers
# (minus the single exclusion from CCFR_4)
#
# remember that you do NOT need to have a reference study indicator variable
# when running models. But I am using a reference study just in case
# GxEScanR complains. The reference study will always be the first
# study in a e main effect dataset (alphabetical, typically ATBC)
#
# looks like additional epi data from lccs and hawaiiccs, confirm with Yi
# (this is specific for NSAIDs since I know from previous runs)
#   - HawaiiCCS_AD N = 1264
#   - LCCS N = 2085
#
#
# IMPORANT
# - this step requires some thought, make sure it conforms to analysis plan 
#   decisions
# - be mindful of covariates
#-----------------------------------------------------------------------------#

wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp")) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude)
  saveRDS(cov, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"), version = 2)
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan.rds"), version = 2)
}


# NSAIDs
wrap(figi_gwas, 'asp_ref')
wrap(figi_gwas, 'nsaids')
wrap(figi_gwas, 'aspirin')

# smoking
wrap(figi_gwas, 'smoke')
wrap(figi_gwas, 'smk_ever')
wrap(figi_gwas, 'smk_pkyrqc2')


# Alcohol
# split by moderate vs heavy drinking
# note that nondrinker counts vary slightly between subsets of moderate and heavy drinking
# because of zero cell count drops (alcoholc_moderate keeps MECC_1 whereas heavy drinking analysis drops it)
wrap(figi_gwas, 'alcoholc_heavy', vars_to_exclude = c())
wrap(figi_gwas, 'alcoholc_moderate', vars_to_exclude = c())

# tmp <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds")
# table(tmp$alcoholc_moderate)
# tmp <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_glm.rds")
# table(tmp$alcoholc_heavy)

# diabetes
wrap(figi_gwas, 'diab', vars_to_exclude = c())

# dietary variables
wrap(figi_gwas, 'calcium_totqc2', vars_to_exclude = c())
wrap(figi_gwas, 'calcium_dietqc2',  vars_to_exclude = c())
wrap(figi_gwas, 'calcium_supp',  vars_to_exclude = c())
wrap(figi_gwas, 'folate_totqc2',  vars_to_exclude = c())
wrap(figi_gwas, 'folate_dietqc2', vars_to_exclude = c())
wrap(figi_gwas, 'folate_sup_yn', vars_to_exclude = c())
wrap(figi_gwas, 'fiberqc2', vars_to_exclude = c())
wrap(figi_gwas, 'vegetableqc2', vars_to_exclude = c())
wrap(figi_gwas, 'fruitqc2',  vars_to_exclude = c())
wrap(figi_gwas, 'redmeatqc2',  vars_to_exclude = c())
wrap(figi_gwas, 'procmeatqc2',  vars_to_exclude = c())


# HRT
wrap(figi_gwas, 'hrt_ref_pm', vars_to_exclude = c("energytot_imp", "sex"))
wrap(figi_gwas, 'eo_ref_pm', vars_to_exclude = c("energytot_imp", "sex"))
wrap(figi_gwas, 'ep_ref_pm', vars_to_exclude = c("energytot_imp", "sex"))

# numeric variables
wrap(figi_gwas, 'bmi', is_e_categorical = F, vars_to_exclude = c())
wrap(figi_gwas, 'bmi5', is_e_categorical = F,  vars_to_exclude = c())
wrap(figi_gwas, 'heightcm', is_e_categorical = F,  vars_to_exclude = c())


# tmp <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_bmi_basic_covars_glm.rds")
# drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>% filter(Freq == 0)
# 
# wtf <- filter(figi_gwas, gxe == 1, !is.na(bmi)) %>% 
#   mutate(outcome = as.numeric(outcome) - 1)
# drops <- data.frame(table(wtf$outcome, wtf$study_gxe)) %>% filter(Freq == 0)
# wtf <- filter(wtf, !study_gxe %in% drops$Var2)


