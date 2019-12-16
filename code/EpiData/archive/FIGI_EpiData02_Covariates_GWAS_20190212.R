#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI Analysis 21/12/2019
#
# Create covariate object for GxEScanR
# 
# Notes:
# - 
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
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")


df <- figi %>% 
  filter(drop == 0, gxe == 1)

table(df$outc)


# for garrett - i'm going to give him outcomes for case/control, dropping 'others'

df <- figi %>% 
  filter(drop == 0, 
         outc != "Other") %>% 
  mutate(outc = ifelse(outc == "Case", 1, 0))

table(df$outc)

write.table(df[, c("vcfid", "outc")], file = "~/Desktop/figi.txt", quote = F, row.names = F, col.names = T)




df <- data.frame(table(Epi$studyname, Epi$study_gxe)) %>% 
  filter(Freq != 0)

#------ GWAS Set (N = 122158) ------
pc30k <- fread("~/Dropbox/code/FIGI_Results/PrincipalComponents/FIGI_GwasSet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

# for GWAS set, be mindful of outc == 'Other' ... drop?
# For this analysis, outcome = crc
cov <- Epi %>%
  filter(drop == 0) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(crc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname))

table(cov$studyname, cov$outc)
table(cov$outc, cov$adenoma_adv)
table(cov$outc, cov$adenoma)
table(cov$outc, cov$crc) # **
table(cov$crc, cov$adenoma_adv)
table(cov$outcome, cov$crc)

case_only <- c("CGN", "FIRE3", "GALEON", "HispanicCCS", "MAVERICC", "MOFFITT", "MSKCC", "NGCCS", "Taiwan", "TRIBE")

cov <- Epi %>%
  filter(drop == 0,
         crc != "",
         !studyname %in% case_only) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(crc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, studyname, paste0(rep("PC", 10), seq(1,10)), sex) %>% 
  filter(complete.cases(.))

for(t in unique(cov$studyname)) {
  cov[paste0(t)] <- ifelse(cov$studyname==t,1,0)
}

# Checks
check <- dplyr::select(cov, unique(cov$studyname)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final
cov <- cov %>% 
  dplyr::select(-UKB_1, -studyname, -sex, sex) 
names(cov)

# save cov file as rds + text
saveRDS(cov, file = '~/FIGI_GWASSet_age_sex_pc10_studyname_122158.rds')



