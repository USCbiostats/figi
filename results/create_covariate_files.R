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
# GWAS set
#
# similar set up as GxE
# use studyname rather than study_gxe (?)
#
#-----------------------------------------------------------------------------#
cov_gwas <- figi_gwas %>% 
  filter(outc != "Other", 
         EUR_subset == 1) %>% 
  mutate(outcome = ifelse(outcome == "Control", 0, 1),
         sex = ifelse(sex == "Female", 0, 1)) %>% 
  dplyr::select(vcfid, outcome, sex, studyname, paste0(rep("PC", 3), seq(1,3)), age_ref_imp) %>% 
  filter(complete.cases(.))

drops <- data.frame(table(cov_gwas$outcome, cov_gwas$studyname)) %>%
  filter(Freq <= 0)
cov_gwas <- filter(cov_gwas, !studyname %in% unique(drops$Var2)) %>%
  dplyr::mutate(studyname = fct_drop(studyname)) 

# indicator variables for GxEScanR
cov_gwas_gxescan <- cov_gwas
ref_study <- as.character(unique(cov_gwas_gxescan[, 'studyname'])[1])
for(t in unique(cov_gwas_gxescan$studyname)) {
  cov_gwas_gxescan[paste0(t)] <- ifelse(cov_gwas_gxescan$studyname==t,1,0)
}

# clean up , make Kentucky reference study_gxe (arbitrary)
cov_gwas_gxescan <- dplyr::select(cov_gwas_gxescan, -ref_study, -studyname, -age_ref_imp, age_ref_imp)
# saveRDS(cov_gwas_gxescan, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gwasset_basic_covars_gxescan.rds", version = 2)





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
# IMPORANT
# - this step requires some thought, make sure it conforms to analysis plan 
#   decisions
# - be mindful of covariates
#-----------------------------------------------------------------------------#

wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp"), studies_to_exclude) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, eur_only = T) %>% 
    filter(!study_gxe %in% exclude) %>% 
    mutate(study_gxe = fct_drop(study_gxe))
  saveRDS(cov, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan.rds"), version = 2)
}




# ------ Alcohol ------
# https://docs.google.com/document/d/19IPCQsD6soQ4E2PnQ3MFuxO0n3gpTOW0VjxkI4ZRRQM/edit#
#
# split by moderate vs heavy drinking
# note that nondrinker counts vary slightly between subsets of moderate and heavy drinking
# because of zero cell count drops (alcoholc_moderate keeps MECC_1 whereas heavy drinking analysis drops it)


# exclusions
exclude <- c("ASTERISK", "ESTHER_VERDI", "PPS3", "PPS4", "PHS")

wrap(figi_gwas, 'alcoholc_moderate', vars_to_exclude = c(), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$alcoholc)

wrap(figi_gwas, 'alcoholc_heavy', vars_to_exclude = c(), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$alcoholc)

wrap(figi_gwas, 'alcoholc', vars_to_exclude = c(), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$alcoholc)

# explore studies dropped due to small cell sizes
# y <- filter(figi_gwas, gxe == 1, EUR_subset == 1, !is.na(alcoholc_moderate)) %>%
#   mutate(outcome = fct_drop(outcome))
# drops <- data.frame(table(y$outcome, y[, 'alcoholc_moderate'], y$study_gxe)) %>%
#   filter(Freq == 0)



# ------ HRT ------
# https://drive.google.com/file/d/1WoUnA_mh6oKb6UAC6-X_-MLEx47yAZ6r/view?usp=sharing
#
# HawaiiCCS has updated data but still high missingness for HRT. exclude for now

# HRT
exclude = c("HawaiiCCS_AD")
wrap(figi_gwas, 'hrt_ref_pm', vars_to_exclude = c("energytot_imp", "sex"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm_full_glm.rds")
table(x$meno, useNA = 'ifany')
summary(as.numeric(x$age_ref))
summary(hist(x$age_ref))


wrap(figi_gwas, 'eo_ref_pm', vars_to_exclude = c("energytot_imp", "sex"))
wrap(figi_gwas, 'ep_ref_pm', vars_to_exclude = c("energytot_imp", "sex"))





# NSAIDs
wrap(figi_gwas, 'asp_ref')
wrap(figi_gwas, 'nsaids')
wrap(figi_gwas, 'aspirin')

# smoking
wrap(figi_gwas, 'smoke')
wrap(figi_gwas, 'smk_ever')
wrap(figi_gwas, 'smk_pkyrqc2')


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



# numeric variables
wrap(figi_gwas, 'bmi', is_e_categorical = F, vars_to_exclude = c())
wrap(figi_gwas, 'bmi5', is_e_categorical = F,  vars_to_exclude = c())
wrap(figi_gwas, 'heightcm', is_e_categorical = F,  vars_to_exclude = c())



#-----------------------------------------------------------------------------#
# same as above, but keep full covariate set
#-----------------------------------------------------------------------------#

wrap_fullcovar <- function(d, exposure) {
  x <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"))
  out <- d %>% 
    filter(vcfid %in% x$vcfid)
  saveRDS(out, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_full_glm.rds"), version = 2)
  
}

# NSAIDs
wrap_fullcovar(figi_gwas, 'asp_ref')
wrap_fullcovar(figi_gwas, 'nsaids')
wrap_fullcovar(figi_gwas, 'aspirin')

# smoking
wrap_fullcovar(figi_gwas, 'smoke')
wrap_fullcovar(figi_gwas, 'smk_ever')
wrap_fullcovar(figi_gwas, 'smk_pkyrqc2')

# Alcohol
# split by moderate vs heavy drinking
# note that nondrinker counts vary slightly between subsets of moderate and heavy drinking
# because of zero cell count drops (alcoholc_moderate keeps MECC_1 whereas heavy drinking analysis drops it)
wrap_fullcovar(figi_gwas, 'alcoholc_heavy')
wrap_fullcovar(figi_gwas, 'alcoholc_moderate')

# diabetes
wrap_fullcovar(figi_gwas, 'diab')

# dietary variables
wrap_fullcovar(figi_gwas, 'calcium_totqc2')
wrap_fullcovar(figi_gwas, 'calcium_dietqc2')
wrap_fullcovar(figi_gwas, 'calcium_supp')
wrap_fullcovar(figi_gwas, 'folate_totqc2')
wrap_fullcovar(figi_gwas, 'folate_dietqc2')
wrap_fullcovar(figi_gwas, 'folate_sup_yn')
wrap_fullcovar(figi_gwas, 'fiberqc2')
wrap_fullcovar(figi_gwas, 'vegetableqc2')
wrap_fullcovar(figi_gwas, 'fruitqc2')
wrap_fullcovar(figi_gwas, 'redmeatqc2')
wrap_fullcovar(figi_gwas, 'procmeatqc2')


# HRT
wrap_fullcovar(figi_gwas, 'hrt_ref_pm')
wrap_fullcovar(figi_gwas, 'eo_ref_pm')
wrap_fullcovar(figi_gwas, 'ep_ref_pm')

# numeric variables
wrap_fullcovar(figi_gwas, 'bmi')
wrap_fullcovar(figi_gwas, 'bmi5')
wrap_fullcovar(figi_gwas, 'heightcm')


