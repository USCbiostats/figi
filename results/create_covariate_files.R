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
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)




#-----------------------------------------------------------------------------#
# GWAS set
#
# similar set up as GxE
# use studyname rather than study_gxe (?)
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

wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp"), vars_to_include = c(), studies_to_exclude) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, vars_to_include, eur_only = T) %>% 
    filter(!study_gxe %in% studies_to_exclude) %>% 
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



## create heavy vs moderate alcohol intake
figi_tmp <- mutate(figi_gxe, 
                   alcoholc_heavy_vs_moderate = ifelse(alcoholc == "1-28g/d", 0,
                                                ifelse(alcoholc == ">28g/d", 1, NA)),
                   alcoholc_heavy_vs_moderate = factor(alcoholc_heavy_vs_moderate, labels = c("1-28g/d", ">28g/d")))
table(figi_tmp$alcoholc_heavy_vs_moderate)

exclude <- c("ASTERISK", "ESTHER_VERDI", "PPS3", "PPS4", "PHS", "MECC_1")
wrap(figi_tmp, 'alcoholc_heavy_vs_moderate', vars_to_exclude = c(), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds")
head(x)
table(x$alcoholc_heavy_vs_moderate)
table(x$study_gxe, x$outcome)







# ------ HRT ------
# https://drive.google.com/file/d/1WoUnA_mh6oKb6UAC6-X_-MLEx47yAZ6r/view?usp=sharing
#
# HawaiiCCS has updated data but still high missingness for HRT. exclude for now

# IMPORTANT - reference group should be "NO ANY HORMONE USE"
table(figi_gxe$hrt_ref_pm, useNA = 'ifany')
table(figi_gxe$eo_ref_pm, useNA = 'ifany')
table(figi_gxe$hrt_ref_pm2, figi_gxe$eo_ref_pm, useNA = 'ifany')
table(figi_gxe$hrt_ref_pm2, figi_gxe$ep_ref_pm, useNA = 'ifany')

# to create reference groups for eo_ref_pm and ep_ref_pm: 
# "yes" from eo_ref_pm, "no", from hrt_ref_pm2
figi_tmp <- mutate(figi_gxe, 
            eo_ref_pm_gxe = ifelse(eo_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                            ifelse(hrt_ref_pm2 == "No", "No", NA)), 
            ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                   ifelse(hrt_ref_pm2 == "No", "No", NA)),
            eo_ref_pm_gxe = factor(eo_ref_pm_gxe),
            ep_ref_pm_gxe = factor(ep_ref_pm_gxe))

table(figi_tmp$eo_ref_pm_gxe)
table(figi_tmp$ep_ref_pm_gxe)
table(hrt_ref_pm2=figi_tmp$hrt_ref_pm2,eo_ref_pm_gxe=figi_tmp$eo_ref_pm_gxe, useNA = 'ifany')
table(hrt_ref_pm2=figi_tmp$hrt_ref_pm2,ep_ref_pm_gxe=figi_tmp$ep_ref_pm_gxe, useNA = 'ifany')


# hrt_ref_pm
exclude = c("HawaiiCCS_AD")
wrap(figi_gwas, 'hrt_ref_pm', vars_to_exclude = c("energytot_imp", "sex"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm_full_glm.rds")
table(x$meno, useNA = 'ifany')
summary(as.numeric(x$age_ref))
summary(hist(x$age_ref))


# hrt_ref_pm2 (definition of controls more in line with analysis plan..)
exclude = c("HawaiiCCS_AD")
wrap(figi_gwas, 'hrt_ref_pm2', vars_to_exclude = c("energytot_imp", "sex"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$hrt_ref_pm2)



# eo_ref_pm_gxe (see above)
exclude = c("HawaiiCCS_AD", "REACH_AD", "SMS_AD")
wrap(figi_tmp, 'eo_ref_pm_gxe', vars_to_exclude = c("energytot_imp", "sex"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

# ep_ref_pm_gxe (see above)
exclude = c("HawaiiCCS_AD", "REACH_AD", "SMS_AD")
wrap(figi_tmp, 'ep_ref_pm_gxe', vars_to_exclude = c("energytot_imp", "sex"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)


wrap(figi_gwas, 'ep_ref_pm_gxe', vars_to_exclude = c("energytot_imp", "sex"))







# ------ NSAIDs ------

# asp_ref
exclude <- c()

wrap(figi_gwas, 'asp_ref', studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_asp_ref_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$asp_ref)


# nsaids
# (you need to redefine the mutually adjusted variables to numeric)
tmp <- figi_gwas %>% 
  mutate(aspirin = as.numeric(aspirin) - 1)

wrap(tmp, 'nsaids', studies_to_exclude = exclude, vars_to_include = 'aspirin')
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_nsaids_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$nsaids)

# aspirin
# N should be identical to that of nsaids
tmp <- figi_gwas %>% 
  mutate(nsaids = as.numeric(nsaids) - 1)

wrap(tmp, 'aspirin', studies_to_exclude = exclude, vars_to_include = 'nsaids')
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_aspirin_basic_covars_glm.rds")
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$nsaids)





# ------ T2D ------

# diab
# recall that studies to exclude are studies that have E information but excluded for other reasons
exclude = c("MECC_1", "MECC_2", "MECC_3")
wrap(figi_gwas, 'diab', vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_diab_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

## Questions:
## 1) ATBC
## 2) MEC_2
## 3) HawaiiCCS - include or not? 

atbc <- figi_gwas %>% 
  filter(!(is.na(diab)), 
           EUR_subset == 1,
         study_gxe == "ATBC")
table(atbc$study_gxe, atbc$outcome, atbc$diab)

mec_2 <- figi_gwas %>% 
  filter(!(is.na(diab)),
         EUR_subset == 1,
         study_gxe == "MEC_2")
table(mec_2$study_gxe, mec_2$outcome, mec_2$diab)





# ------ red/processed meat ------

# redmeatqc2
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, 'redmeatqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

# some studies have only two values in the quartile variable (2,3)
# these need to be included in the analysis
asterisk <- filter(figi_gwas, 
                   !(is.na(redmeatqc2)), 
                   EUR_subset == 1,
                   study_gxe == "ASTERISK")
table(asterisk$outcome, asterisk$redmeatqc2)


colo23 <- filter(figi_gwas, 
                 !(is.na(redmeatqc2)), 
                 EUR_subset == 1,
                 study_gxe == "Colo23")
table(colo23$outcome, colo23$redmeatqc2)


# to add the above studies, let me just fix it here directly with rbinds
tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(redmeatqc2), 
                EUR_subset == 1,
                study_gxe %in% c("ASTERISK", "Colo23")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                redmeatqc2 = as.numeric(redmeatqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, redmeatqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'redmeatqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "redmeatqc2", "_basic_covars_gxescan.rds"), version = 2)

# final check
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_gxescan.rds")
head(y)


# procmeatqc2
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, 'procmeatqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_procmeatqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

# Colo23, HPFS, MECC, NHS, 
colo23 <- filter(figi_gwas, 
                 !(is.na(procmeatqc2)), 
                 EUR_subset == 1,
                 study_gxe == "Colo23")
table(colo23$outcome, colo23$procmeatqc2)

hpfs <- filter(figi_gwas, 
                 !(is.na(procmeatqc2)), 
                 EUR_subset == 1,
                 grepl("HPFS", study_gxe))
table(hpfs$outcome, hpfs$procmeatqc2)

mecc <- filter(figi_gwas, 
               !(is.na(procmeatqc2)), 
               EUR_subset == 1,
               grepl("MECC", study_gxe))
table(mecc$outcome, mecc$procmeatqc2)

nhs <- filter(figi_gwas, 
               !(is.na(procmeatqc2)), 
               EUR_subset == 1,
               grepl("NHS", study_gxe))
table(nhs$outcome, nhs$procmeatqc2)



# to add the above studies, let me just fix it here directly with rbinds
tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(procmeatqc2), 
                EUR_subset == 1,
                study_gxe %in% c("Colo23", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MECC_1", "MECC_2", "MECC_3", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                procmeatqc2 = as.numeric(procmeatqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, procmeatqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_procmeatqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'procmeatqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "procmeatqc2", "_basic_covars_gxescan.rds"), version = 2)


# final check
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_procmeatqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)







# ------ fiber/fruit/veg ------

# 
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, 'redmeatqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)









# ------ Folate ------

# folate_totqc2
exclude = c("UKB_1", "SMS_AD")
wrap(figi_gwas, 'folate_totqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_totqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

wtf <- filter(figi_gwas, study_gxe == "ASTERISK")
wtf <- filter(figi_gwas, grepl("DACH", study_gxe))
wtf <- filter(figi_gwas, grepl("PHS", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("REACH_AD", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("SMS", study_gxe)) ## **
table(wtf$folate_totqc2)

tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(folate_totqc2), 
                EUR_subset == 1,
                study_gxe %in% c("PHS", "REACH_AD")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                folate_totqc2 = as.numeric(folate_totqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, folate_totqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_totqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'folate_totqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "folate_totqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_totqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$folate_totqc2)



# folate_dietqc2
exclude = c("UKB_1", "SMS_AD")
wrap(figi_gwas, 'folate_dietqc2', vars_to_exclude = c(""), vars_to_include = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_dietqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

wtf <- filter(figi_gwas, study_gxe == "ASTERISK")
wtf <- filter(figi_gwas, grepl("DACH", study_gxe))
wtf <- filter(figi_gwas, grepl("PHS", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("REACH_AD", study_gxe)) ## **
wtf <- filter(figi_gwas, grepl("SMS", study_gxe)) ## **
table(wtf$folate_dietqc2)








# ------ Calcium ------
#(COME BACK TO THIS)
#(ASK YI ABOUT CCFR CALCIUM TOT VS DIET)


# calcium_totqc2
exclude = c("ATBC", "PPS3", "PPS4", "MCCS_1", "MCCS_2")
wrap(figi_gwas, 'calcium_totqc2', vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_totqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

ccfr <- filter(figi_gwas, 
                 !(is.na(calcium_totqc2)), 
                 EUR_subset == 1,
                 grepl("CCFR", study_gxe))
table(ccfr$outcome, ccfr$calcium_totqc2)

reach <- filter(figi_gwas, 
               !(is.na(calcium_totqc2)), 
               EUR_subset == 1,
               grepl("REACH", study_gxe))
table(reach$outcome, reach$calcium_totqc2)

sms_ad <- filter(figi_gwas, 
                !(is.na(calcium_totqc2)), 
                EUR_subset == 1,
                grepl("SMS_AD", study_gxe))
table(sms_ad$outcome, sms_ad$calcium_totqc2)

ukb <- filter(figi_gwas, 
                 !(is.na(calcium_totqc2)), 
                 EUR_subset == 1,
                 grepl("UKB", study_gxe))
table(ukb$outcome, ukb$calcium_totqc2)


# to add the above studies, let me just fix it here directly with rbinds
tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(calcium_totqc2), 
                EUR_subset == 1,
                study_gxe %in% c("CCFR_1", "CCFR_4", "CCFR_3", "REACH_AD", "UKB_1")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                calcium_totqc2 = as.numeric(calcium_totqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, calcium_totqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_totqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'calcium_totqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "calcium_totqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_totqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)





# calcium_dietqc2
exclude = c("ATBC", "PPS3", "PPS4")
wrap(figi_gwas, 'calcium_dietqc2', vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_dietqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)

reach <- filter(figi_gwas, 
              !(is.na(calcium_dietqc2)), 
              EUR_subset == 1,
              grepl("REACH_AD", study_gxe))
table(reach$outcome, reach$calcium_dietqc2)

ukb <- filter(figi_gwas, 
              !(is.na(calcium_dietqc2)), 
              EUR_subset == 1,
              grepl("UKB", study_gxe))
table(ukb$outcome, ukb$calcium_dietqc2)


tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(calcium_dietqc2), 
                EUR_subset == 1,
                study_gxe %in% c("CCFR_1", "CCFR_4", "CCFR_3", "REACH_AD", "UKB_1")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                calcium_dietqc2 = as.numeric(calcium_dietqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, calcium_dietqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_dietqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'calcium_dietqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "calcium_dietqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_calcium_dietqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)







# ------ Smoking ------

# smk_ever
exclude = c("ATBC", "PPS3", "PPS4", "NGCCS", "MECC_3")
wrap(figi_gwas, 'smk_ever', vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_smk_ever_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$smk_ever)

nshds <- filter(figi_gwas, 
               !(is.na(smk_ever)), 
               EUR_subset == 1,
               grepl("NSHDS", study_gxe))
table(nshds$outcome, nshds$smk_ever)




# smk_aveday
# this is a continuous variable.. 
# no smk_aveday should be zero among former and current smokers by definition
# ~ 62 missing smk_aveday counts
x <- figi_gxe %>% 
  filter(smoke != "Never smoked")

describeBy(x$smk_aveday, x$smoke)
summary(x$smk_aveday)
hist(x$smk_aveday)

xx <- filter(x, smk_aveday <= 0)


exclude = c("PPS3", "PPS4", "NGCCS", "MECC_3")
wrap(x, exposure = 'smk_aveday', is_e_categorical = F, vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_smk_aveday_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$smk_aveday)


# smk_pkyr
# another continuous variable
# 125 with zero pkyr
x <- figi_gxe %>% 
  filter(smoke != "Never smoked")

describeBy(x$smk_pkyr, x$smoke)
summary(x$smk_pkyr)
hist(x$smk_pkyr)

xx <- filter(x, smk_pkyr <= 0)

exclude = c("PPS3", "PPS4", "NGCCS", "MECC_3")
wrap(x, exposure = 'smk_pkyr', is_e_categorical = F, vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_smk_pkyr_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$smk_aveday)



# ----- Fruits vegetables and fiber ------ 

# fiberqc2
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, exposure = 'fiberqc2', is_e_categorical = T, vars_to_exclude = c(""), studies_to_exclude = exclude)
y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_fiberqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)
table(y$study_gxe, y$fiberqc2)
table(figi_gwas$study_gxe, figi_gwas$fiberqc2)


# vegetableqc2
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, exposure = 'vegetableqc2', is_e_categorical = T, vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_vegetableqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$vegetableqc2)
table(figi_gwas$study_gxe, figi_gwas$vegetableqc2)

dachs <- filter(figi_gwas, 
                !(is.na(vegetableqc2)), 
                EUR_subset == 1,
                grepl("DACHS", study_gxe))
table(dachs$outcome, dachs$vegetableqc2)

asterisk <- filter(figi_gwas, 
              !(is.na(vegetableqc2)), 
              EUR_subset == 1,
              grepl("ASTERISK", study_gxe))
table(asterisk$outcome, asterisk$vegetableqc2)


tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(vegetableqc2), 
                EUR_subset == 1,
                study_gxe %in% c("DACHS_1", "DACHS_2", "DACHS_3", "ASTERISK")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                vegetableqc2 = as.numeric(vegetableqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, vegetableqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_vegetableqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'vegetableqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "vegetableqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_vegetableqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)


# fruitqc2
exclude = c("PPS3", "PPS4")
wrap(figi_gwas, exposure = 'fruitqc2', is_e_categorical = T, vars_to_exclude = c(""), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$fruitqc2)
table(figi_gwas$study_gxe, figi_gwas$fruitqc2)

dachs <- filter(figi_gwas, 
                !(is.na(fruitqc2)), 
                EUR_subset == 1,
                grepl("DACHS", study_gxe))
table(dachs$outcome, dachs$fruitqc2)

asterisk <- filter(figi_gwas, 
                   !(is.na(fruitqc2)), 
                   EUR_subset == 1,
                   grepl("ASTERISK", study_gxe))
table(asterisk$outcome, asterisk$fruitqc2)

nhs5 <- filter(figi_gwas, 
                   !(is.na(fruitqc2)), 
                   EUR_subset == 1,
                   grepl("NHS_5", study_gxe))
table(nhs5$outcome, nhs5$fruitqc2)

tmp <- figi_gxe %>% 
  dplyr::filter(!is.na(fruitqc2), 
                EUR_subset == 1,
                study_gxe %in% c("DACHS_1", "DACHS_2", "DACHS_3", "ASTERISK")) %>% 
  dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
                sex = ifelse(sex == "Female", 0, 1),
                fruitqc2 = as.numeric(fruitqc2) - 1 ) %>% 
  dplyr::select(vcfid, outcome, fruitqc2, age_ref_imp, sex, energytot_imp, study_gxe, PC1, PC2, PC3)

out <- rbind(x, tmp)
saveRDS(out, "~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds")

cov_gxescan <- format_data_gxescan(out, 'fruitqc2')
saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "fruitqc2", "_basic_covars_gxescan.rds"), version = 2)

y <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds")
table(y$study_gxe, y$outcome)




# ----- BMI5 ------ 

# bmi5
exclude = c("")
wrap(figi_gwas, exposure = 'bmi5', is_e_categorical = F, vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_bmi5_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$bmi5)


# ----- Height -----

# height10
exclude = c("")
wrap(figi_gwas, exposure = 'height10', is_e_categorical = F, vars_to_exclude = c("energytot_imp"), studies_to_exclude = exclude)
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_height10_basic_covars_glm.rds")
table(x$study_gxe, x$outcome)
table(x$study_gxe, x$bmi5)









# numeric variables
wrap(figi_gwas, 'bmi', is_e_categorical = F, vars_to_exclude = c())
wrap(figi_gwas, 'bmi5', is_e_categorical = F,  vars_to_exclude = c())
wrap(figi_gwas, 'heightcm', is_e_categorical = F,  vars_to_exclude = c())










#-----------------------------------------------------------------------------#
# same as above, but keep full covariate set
#-----------------------------------------------------------------------------#

# wrap_fullcovar <- function(d, exposure) {
#   x <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"))
#   out <- d %>% 
#     filter(vcfid %in% x$vcfid)
#   saveRDS(out, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_full_glm.rds"), version = 2)
#   
# }
# 
# # NSAIDs
# wrap_fullcovar(figi_gwas, 'asp_ref')
# wrap_fullcovar(figi_gwas, 'nsaids')
# wrap_fullcovar(figi_gwas, 'aspirin')
# 
# # smoking
# wrap_fullcovar(figi_gwas, 'smoke')
# wrap_fullcovar(figi_gwas, 'smk_ever')
# wrap_fullcovar(figi_gwas, 'smk_pkyrqc2')
# 
# # Alcohol
# # split by moderate vs heavy drinking
# # note that nondrinker counts vary slightly between subsets of moderate and heavy drinking
# # because of zero cell count drops (alcoholc_moderate keeps MECC_1 whereas heavy drinking analysis drops it)
# wrap_fullcovar(figi_gwas, 'alcoholc_heavy')
# wrap_fullcovar(figi_gwas, 'alcoholc_moderate')
# 
# # diabetes
# wrap_fullcovar(figi_gwas, 'diab')
# 
# # dietary variables
# wrap_fullcovar(figi_gwas, 'calcium_totqc2')
# wrap_fullcovar(figi_gwas, 'calcium_dietqc2')
# wrap_fullcovar(figi_gwas, 'calcium_supp')
# wrap_fullcovar(figi_gwas, 'folate_totqc2')
# wrap_fullcovar(figi_gwas, 'folate_dietqc2')
# wrap_fullcovar(figi_gwas, 'folate_sup_yn')
# wrap_fullcovar(figi_gwas, 'fiberqc2')
# wrap_fullcovar(figi_gwas, 'vegetableqc2')
# wrap_fullcovar(figi_gwas, 'fruitqc2')
# wrap_fullcovar(figi_gwas, 'redmeatqc2')
# wrap_fullcovar(figi_gwas, 'procmeatqc2')
# 
# 
# # HRT
# wrap_fullcovar(figi_gwas, 'hrt_ref_pm')
# wrap_fullcovar(figi_gwas, 'eo_ref_pm')
# wrap_fullcovar(figi_gwas, 'ep_ref_pm')
# 
# # numeric variables
# wrap_fullcovar(figi_gwas, 'bmi')
# wrap_fullcovar(figi_gwas, 'bmi5')
# wrap_fullcovar(figi_gwas, 'heightcm')


