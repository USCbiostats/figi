#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
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

pca <- "/home/rak/data/PCA/190506/FIGI_GxESet_KGP_pc20_190430.eigenvec"
# filename <- "/home/rak/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_folate_dietqc2_sex_age_pc3_studygxe_"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")

#-----------------------------------------------------------------------------#
# GxE Set - folate_dietqc2
#-----------------------------------------------------------------------------#
pc30k <- fread(pca, skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>%
  dplyr::filter(drop == 0 & gxe == 1) %>% 
  dplyr::inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  dplyr::mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         energytot = as.numeric(energytot), 
         folate_dietqc2 = dplyr::recode(folate_dietqc2, `1`=1, `2`=2, `3`=3, `4`=4)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3)), energytot, folate_dietqc2) %>% 
  dplyr::filter(complete.cases(.))

# ------ exclude studies (case only etc) ------ #
table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)
exclude_studies <- as.vector(unique(drops$Var1))

cov <- filter(cov, !study_gxe %in% exclude_studies)

# ----- save, no dummy var ------
# saveRDS(cov, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72820_GLM.rds", version = 2)


for(t in unique(cov$study_gxe)) {
  cov[paste0(t)] <- ifelse(cov$study_gxe==t,1,0)
}

# Checks (what studies are included and how many)
check <- dplyr::select(cov, unique(cov$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final (make Kentucky the reference)
cov <- cov %>% 
  dplyr::select(-Kentucky, -study_gxe,
                -folate_dietqc2, folate_dietqc2)

# ------ save, WITH dummy var ------
saveRDS(cov, file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_folate_dietqc2_sex_age_pc3_energytot_studygxe_52447.rds', version = 2)


