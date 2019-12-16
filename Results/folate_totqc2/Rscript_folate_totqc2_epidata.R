#=============================================================================#
# FIGI Analysis 12/09/2019 - HRC v2.3
#
# Important:
# MAKE SURE TO INCORPORATE THIS SCRIPT AFTER CALCULATING PCs
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")

#-----------------------------------------------------------------------------#
# GxE Set asp_ref
#-----------------------------------------------------------------------------#
cov_gxe <- figi %>%
  filter(drop == 0 & gxe == 1)

# checks
xtabs(~ cov_gxe$folate_totqc2, addNA = T)
xtabs(~ cov_gxe$sex, addNA = T)
xtabs(~ cov_gxe$study_gxe, addNA = T)
xtabs(~ cov_gxe$outc + cov_gxe$folate_totqc2, addNA = T)

# remove missing asp_ref, recode covariates
cov_full <- cov_gxe %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         folate_totqc2 = as.numeric(folate_totqc2))

cov <- cov_full %>% 
  filter(!is.na(folate_totqc2)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, folate_totqc2) %>% 
  filter(complete.cases(.))

table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))

cov_full_subset <- cov_full %>% 
  filter(vcfid %in% cov$vcfid)


# create indicator variables for GxEScanR
cov_gxescan <- cov
for(t in unique(cov_gxescan$study_gxe)) {
  cov_gxescan[paste0(t)] <- ifelse(cov_gxescan$study_gxe==t,1,0)
}

# checks (what studies are included and how many)
check <- dplyr::select(cov_gxescan, unique(cov_gxescan$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# clean up , make Kentucky reference study_gxe (arbitrary)
cov_gxescan <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -asp_ref, asp_ref)

# save files
# saveRDS(cov,         file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72269_GLM.rds', version = 2)
# saveRDS(cov_gxescan, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72269.rds",     version = 2)

# save files
# remember that you only recoded outcome, asp_ref, sex, etc. So yes it contains all covariates but it's stll specific to the exposure in question. I'll create individual files now but eventually can create a single global, cleaned E dataset
saveRDS(cov_full_subset, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_totqc2_age_ref_imp_sex_study_gxe_full.rds', version = 2)
saveRDS(cov,             file = '~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_folate_totqc2_age_ref_imp_sex_study_gxe_glm.rds', version = 2)
saveRDS(cov_gxescan,     file = "~/data/GxEScanR_PhenotypeFiles/FIGI_gxeset_asp_ref_age_ref_imp_sex_study_gxe_pc1_pc2_pc3.rds",      version = 2)

