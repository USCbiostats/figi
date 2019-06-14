#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI Analysis 01/17/2019
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
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(data.table)

rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")

# sanity
# any(duplicated(Epi$vcfid))
# any(duplicated(Epi$pooledcompassid))
# table(Epi$studyname)
# table(Epi$V2)


#------ GxE Set - asp_ref (N = 102,792 --> 72,438) ------
pc30k <- fread("~/Dropbox/code/FIGI_Results/PrincipalComponents/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, studyname, paste0(rep("PC", 10), seq(1,10)), asp_ref) %>% 
  filter(complete.cases(.))

table(cov$studyname, cov$outcome)
exclude_studies <- c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT") # case-only studies.. 

cov <- filter(cov, !studyname %in% exclude_studies)

for(t in unique(cov$studyname)) {
  cov[paste0(t)] <- ifelse(cov$studyname==t,1,0)
}

# Checks
check <- dplyr::select(cov, unique(cov$studyname)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final
cov <- cov %>% 
  dplyr::select(-UKB_1, -studyname,
                -asp_ref, asp_ref) 
  
# test run
# formula_text = paste0('outcome ~ sex + age_ref_imp + ', paste(unique(cov$studyname), collapse = "+"))
# glm(as.formula(formula_text) , data = cov, family = 'binomial')


# save cov file as rds + text
saveRDS(cov, file = '~/FIGI_GxESet_asp_ref_sex_age_pcs_studyname_72438.rds')
write.table(cov, file = "~/FIGI_GxESet_asp_ref_sex_age_pcs_studyname_72438.cov", quote = F, row.names = F, col.names = T, sep = '\t')




#------ GxE Set - aspirin (N = 102,792 --> 71,783) ------
pc30k <- fread("~/Dropbox/code/FIGI_PCA_IBD_Results/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname),
         aspirin = ifelse(aspirin == "Yes", 1, ifelse(aspirin == "No", 0, NA))) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, studyname, paste0(rep("PC", 10), seq(1,10)), aspirin) %>% 
  filter(complete.cases(.))

table(cov$studyname, cov$outcome)
exclude_studies <- c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT") # case-only studies.. 

cov <- filter(cov, !studyname %in% exclude_studies)

for(t in unique(cov$studyname)) {
  cov[paste0(t)] <- ifelse(cov$studyname==t,1,0)
}

# Checks
check <- dplyr::select(cov, unique(cov$studyname)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final
cov <- cov %>% 
  dplyr::select(-UKB_1, -studyname,
                -aspirin, aspirin) 

# test run
# formula_text = paste0('outcome ~ sex + age_ref_imp + ', paste(unique(cov$studyname), collapse = "+"))
# glm(as.formula(formula_text) , data = cov, family = 'binomial')


# save cov file as rds + text
saveRDS(cov, file = '~/FIGI_GxESet_aspirin_sex_age_pcs_studyname_71783.rds')
write.table(cov, file = "~/FIGI_GxESet_aspirin_sex_age_pcs_studyname_71783.cov", quote = F, row.names = F, col.names = T, sep = '\t')





#------ GxE Set - sex (N = 102,792 --> 72,438) ------
pc30k <- fread("~/Dropbox/code/FIGI_PCA_IBD_Results/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, studyname, paste0(rep("PC", 10), seq(1,10)), sex) %>% 
  filter(complete.cases(.))

table(cov$studyname, cov$outcome)
exclude_studies <- c("GALEON", "MOFFITT", "NGCCS") # case-only studies.. 

cov <- filter(cov, !studyname %in% exclude_studies)

for(t in unique(cov$studyname)) {
  cov[paste0(t)] <- ifelse(cov$studyname==t,1,0)
}

# Checks
check <- dplyr::select(cov, unique(cov$studyname)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final
cov <- cov %>% 
  dplyr::select(-UKB_1, -studyname,
                -sex, sex) 

# test run
# formula_text = paste0('outcome ~ sex + age_ref_imp + ', paste(unique(cov$studyname), collapse = "+"))
# glm(as.formula(formula_text) , data = cov, family = 'binomial')


# save cov file as rds + text
saveRDS(cov, file = '~/FIGI_GxESet_sex_age_pcs_studyname_101171.rds')
write.table(cov, file = "~/FIGI_GxESet_sex_age_pcs_studyname_101171.cov", quote = F, row.names = F, col.names = T, sep = '\t')
