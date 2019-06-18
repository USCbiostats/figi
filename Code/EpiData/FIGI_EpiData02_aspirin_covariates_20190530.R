#=============================================================================#
# FIGI Analysis 02/26/2019
# FIGI Analysis 03/13/2019
# FIGI Analysis 05/16/2019
#
# Creates the following:
# - Covariate table for GxEScanR (flat text table)
# - Covariate data.frame (same as above, but study_gxe as factor)
# - Covariate data.frame with more covariates (for post-harmonization, post-hoc analyses)
#
# Notes:
# use gxe set as determined in `gxe` variable (in addition to drop == 0)
# 	- gxe == 1 removes nonEUR studies --- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
#   - study_gxe = adj. for study + platform. Remove case-only studies when running GxEScanR (causes crash)
# principal components
#   - 2 sets available: GWAS and GxE
#   - see documentation for details
#
# maybe generalize this step as a function but for now edit manually
# also need to create datasets by tumor site
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
pca <- "/home/rak/data/PCA/190506/FIGI_GxESet_KGP_pc20_190430.eigenvec"
filename <- "/home/rak/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485"
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")


#-----------------------------------------------------------------------------#
# GxE Set aspirin
#-----------------------------------------------------------------------------#
pc30k <- fread(pca, skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

cov <- figi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         aspirin = ifelse(aspirin == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study_gxe, paste0(rep("PC", 3), seq(1,3)), aspirin) %>% 
  filter(complete.cases(.))

table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)
exclude_studies <- as.vector(unique(drops$Var1))

cov <- filter(cov, !study_gxe %in% exclude_studies)


# ------ covar data.frame (minimal)  ------
# no study_gxe indicator variables (as factors)
saveRDS(cov, file = paste0(filename, "_GLM.rds"), version = 2)


# ------ GxEScan covariate table ------
out <- cov

for(t in unique(out$study_gxe)) {
  out[paste0(t)] <- ifelse(out$study_gxe==t,1,0)
}

# Checks (what studies are included and how many)
check <- dplyr::select(out, unique(out$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# final (making Kentucky the reference)
out <- out %>% 
  dplyr::select(-Kentucky, -study_gxe,
                -aspirin, aspirin)

saveRDS(out, file = paste0(filename, '.rds'), version = 2)


