#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
# FIGI Analysis 08/08/2019 - reimputed UKB + reach
#
# Create covariate tables for GxEScanR
#
# COVARIATE DATA.FRAME
# 
# Notes:
# use gxe set as determined in `gxe` variable (in addition to drop == 0)
# 	- gxe == 1 removes nonEUR studies
#		- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
#
# studyname = adj. for study + platform
# remove case-only studies (causes crash)
# add principal components
# 	- remember these were calculated on whole set N = 141,362
#   - updated PCs, 2 sets available: GWAS and GxE
#
# Combine NFCCR_1 NFCCR_2 (double check this with Yi)
#		- NFCCR_1 only has cases, NFCCR_2 has cases/controls 
#   - actually, DON'T DO THIS!!!!
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
pca <- "/home/rak/data/PCA/190729/FIGI_GxESet_190729.eigenvec"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")


# ------ GxE Set ------ #
pc30k <- fread(pca, skip = 1, col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>% 
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0), 
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1), 
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         bmi5 = as.numeric(bmi5)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 10), seq(1,10)), bmi5) %>% 
  filter(!is.na(bmi5)) # same as complete.cases(.) FYI

table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))


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
cov_gxescan <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -bmi5, bmi5)

# save files
saveRDS(cov,         file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_sex_age_pc10_studygxe_88324_GLM.rds', version = 2)
saveRDS(cov_gxescan, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_sex_age_pc10_studygxe_88324.rds",     version = 2)





# ------ GxE Set - BY SEX ------ #
# keep in mind studies that only have one sex

# ------- #
# females 
cov <- figi %>% 
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0), 
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1), 
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         bmi5 = as.numeric(bmi5)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 10), seq(1,10)), bmi5) %>% 
  filter(!is.na(bmi5),
         sex == 0) %>% 
  dplyr::select(-sex)

# drop case only studies
table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))

cov_gxescan <- cov
for(t in unique(cov_gxescan$study_gxe)) {
  cov_gxescan[paste0(t)] <- ifelse(cov_gxescan$study_gxe==t,1,0)
}

# checks (what studies are included and how many)
check <- dplyr::select(cov_gxescan, unique(cov_gxescan$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# clean up , make Kentucky reference study_gxe (arbitrary)
cov_gxescan_female <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -bmi5, bmi5)
saveRDS(cov_gxescan_female, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_age_pc10_studygxe_42558_female.rds', version = 2)



# ----- #
# Males 
cov <- figi %>% 
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0), 
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1), 
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         bmi5 = as.numeric(bmi5)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 10), seq(1,10)), bmi5) %>% 
  filter(!is.na(bmi5),
         sex == 1) %>% 
  dplyr::select(-sex)

# drop case only studies
table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))

cov_gxescan <- cov
for(t in unique(cov_gxescan$study_gxe)) {
  cov_gxescan[paste0(t)] <- ifelse(cov_gxescan$study_gxe==t,1,0)
}

# checks (what studies are included and how many)
check <- dplyr::select(cov_gxescan, unique(cov_gxescan$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# clean up , make Kentucky reference study_gxe (arbitrary)
cov_gxescan_male <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -bmi5, bmi5)
saveRDS(cov_gxescan_male, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_age_pc10_studygxe_45766_male.rds', version = 2)


# ----- #
# check

x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_age_pc10_studygxe_42558_female.rds")
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_age_pc10_studygxe_45766_male.rds")

system("scp ~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_bmi5_age_pc10_studygxe* hpc:/staging/dvc/andreeki/GxEScanR")
