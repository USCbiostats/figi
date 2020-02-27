#=============================================================================#
# FIGI Analysis 12/10/2019
# create master analytical dataset for GWAS and GWIS
#
# Use principal components calculated 190729
# - excludes 557 TCGA samples (genotype data N/A)
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GwasSet_190729.eigenvec"
pc <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")


#-----------------------------------------------------------------------------#
# GWAS SET
#
# notes for now
# - remind me again criteria for determinig gxe == 1. some studies are
#   gxe == 0 but still have some exposure information like aspirin
# - calcium supplemental - which variable to use if they want to analyze this var
#-----------------------------------------------------------------------------#
figi_gwas <- inner_join(figi, pc, by = c("vcfid" = "IID")) %>% 
  dplyr::filter(drop == 0) %>% 
  dplyr::mutate(outcome = factor(outc),
                age_ref_imp = as.numeric(age_ref_imp),
                sex = factor(sex, labels = c("Female", "Male")), 
                famhx1 = factor(famhx1), 
                energytot = as.numeric(energytot), 
                energytot_imp = as.numeric(energytot_imp), 
                educ = factor(educ, labels = c("Less than High School", "High School/GED", "Some College", "College/Graduate School")), 
                studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname), 
                study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe), 
                # create factors for study name variables
                studyname_fct = factor(studyname), # use in gwas analysis
                study_gxe_fct = factor(study_gxe), # use in gxe  analysis
                # studyname_fct = fct_recode(studyname_fct, Colo23 = "Colo2&3"), 
                # study_gxe_fct = fct_recode(study_gxe_fct, Colo23 = "Colo2&3"), 
                # nsaids
                asp_ref = factor(asp_ref),
                nsaids = factor(nsaids),
                aspirin = factor(aspirin),
                # bmi
                bmi = as.numeric(bmi), 
                bmi5 = as.numeric(bmi5), 
                # height
                heightcm = as.numeric(heightcm),
                height10 = as.numeric(height10),
                # smoking
                smoke = factor(smoke),
                smoke = fct_relevel(smoke, "Never smoker", "Former smoker", "Smoker"), 
                smk_ever = factor(smk_ever), 
                smk_pkyrqc2 = factor(smk_pkyrqc2), 
                smk_aveday = as.numeric(smk_aveday),
                smk_pkyr = as.numeric(smk_pkyr),
                # alcohol
                # modified so moderate drinking is ref in alcoholc_moderate
                alcoholc = factor(alcoholc), # keep original labels
                alcoholc = fct_relevel(alcoholc, "nondrinker", "1-28g/d", ">28g/d"),
                alcoholc_heavy = factor(ifelse(alcoholc == "1-28g/d", NA, alcoholc), labels = c("nondrinker", ">28g/d")),
                alcoholc_heavy = fct_drop(alcoholc_heavy),
                alcoholc_moderate = factor(ifelse(alcoholc == ">28g/d", NA, alcoholc), labels = c("nondrinker", "1-28g/d")),
                alcoholc_moderate = fct_drop(alcoholc_moderate), 
                alcoholc_moderate = fct_relevel(alcoholc_moderate, "1-28g/d", "nondrinker"),
                alcoholc_heavy_vs_moderate = ifelse(alcoholc == "1-28g/d", 0,
                                             ifelse(alcoholc == ">28g/d", 1, NA)),
                alcoholc_heavy_vs_moderate = factor(alcoholc_heavy_vs_moderate, labels = c("1-28g/d", ">28g/d")),
                # HRT
                hrt_ref_pm = factor(hrt_ref_pm), 
                hrt_ref_pm2 = factor(hrt_ref_pm2), 
                eo_ref_pm = factor(eo_ref_pm),
                ep_ref_pm = factor(ep_ref_pm),
                eo_ref_pm_gxe = ifelse(eo_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA)), 
                ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                       ifelse(hrt_ref_pm2 == "No", "No", NA)),
                eo_ref_pm_gxe = factor(eo_ref_pm_gxe),
                ep_ref_pm_gxe = factor(ep_ref_pm_gxe),
                # diabetes
                diab = factor(diab), 
                # calcium
                calcium_totqc2 = factor(calcium_totqc2),
                calcium_dietqc2 = factor(calcium_dietqc2),
                calcium_supp = factor(calcium_supp),
                # folate
                folate_totqc2 = factor(folate_totqc2),
                folate_dietqc2 = factor(folate_dietqc2),
                folate_sup_yn = factor(folate_sup_yn), 
                # fiber
                fiberqc2 = factor(fiberqc2),
                # vegetables
                vegetableqc2 = factor(vegetableqc2),
                # fruits
                fruitqc2 = factor(fruitqc2), 
                # meat
                redmeatqc2 = factor(redmeatqc2),
                procmeatqc2 = factor(procmeatqc2), 
                # helper variables to generate table1 with t-test/chisq p values
                table1_outcome = factor(outcome, exclude = NULL, levels = c("Case", "Control", "P-value")),
                table1_outcome = fct_relevel(table1_outcome, "Control", "Case"), # ok to ignore "other" level, gets excluded with gxe == 1
                table1_asp_ref = factor(asp_ref, exclude = NULL, levels = c("No", "Yes", "P-value")))
 



# # test glm
# test_char <- glm(outcome ~ asp_ref + age_ref_imp, data = figi, family = 'binomial')
# summary(test)
# 
# test_factor <- glm(outcome ~ asp_ref + age_ref_imp, data = figi_gwas, family = 'binomial')
# summary(test)
# 
# # test table1
# table1(~ asp_ref + age_ref_imp | outcome, data=figi_gwas)

saveRDS(figi_gwas, "~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds", version = 2)


#-----------------------------------------------------------------------------#
# csv file for jim
#-----------------------------------------------------------------------------#
# rm(list = ls())
# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")
# 
# write.csv(figi_gwas, file = "/media/work/tmp/FIGI_v2.3_gwas_set_N138014_20191211.csv", quote = T, row.names = F)






