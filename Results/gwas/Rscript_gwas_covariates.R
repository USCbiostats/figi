#=============================================================================#
# FIGI Analysis 12/10/2019
#
# create master analytical dataset for GWAS
#
# ------ Notes ------ #
# Principal Components:
# - use PCs calculated 07/29/2019
#
# use gxe set as determined in `gxe` variable (in addition to drop == 0)
# 	- gxe == 1 removes nonEUR studies
#		- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
#
# studyname = adj. for study + platform in GWAS set
# study_gxe = adj. for study + platform in GxE  set
#
# exclude case-only studies
# exclude TCGA (N = 557)
#   - not included in PC calculation, genotype data not available
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
# recode variables as needed, 
#-----------------------------------------------------------------------------#





#-----------------------------------------------------------------------------#
# GWAS SET - age as E variable
#-----------------------------------------------------------------------------#

# discrepancy between 'figi' and 'pc' = 557 --- TCGA samples, which are not included in analysis due to permission reasons
# remove outc == "Other" (adenoma/adv_adenoma, and 631 samples with missing outcome, unsure what they are)
# also remove "corect_oncoarray_nonEUR_reimpute"
# 3 PCs are sufficient according to same metrics used for GxE scans

# wtf <- anti_join(figi, pc, by = c("vcfid" = "IID"))
# wtf <- anti_join(pc, figi, by = c("IID" = "vcfid")) # 427 ?? 

# let me see why anti_join is so many, e.g. what are these files
# load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
# what <- filter(figi, vcfid %in% wtf$IID)


cov_gxe <- figi %>%
  filter(drop == 0) %>% 
  inner_join(pc, by = c('vcfid' = 'IID')) 

write.csv(cov_gxe, file = "~/FIGI_v2.3_gwas_set_N137588_20191209.csv", quote = T, row.names = F)
testing <- fread("~/FIGI_v2.3_gwas_set_N137588_20191209.csv")

xtabs(~ cov_gxe$sex, addNA = T)
xtabs(~ cov_gxe$study_gxe, addNA = T)
xtabs(~ cov_gxe$outc)

cov <- cov_gxe %>% 
  filter(outc != "Other", 
         filename != "corect_oncoarray_nonEUR_reimpute") %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe)) %>% 
  dplyr::select(vcfid, outcome, sex, study_gxe, paste0(rep("PC", 3), seq(1,3)), age_ref_imp) %>% 
  filter(complete.cases(.))

table(cov$study_gxe, cov$outcome)
sort(unique(cov$study_gxe))
drops <- data.frame(table(cov$study_gxe, cov$outcome)) %>% 
  filter(Freq == 0)

cov <- filter(cov, !study_gxe %in% unique(drops$Var1))
table(cov$study_gxe, cov$outcome)

# indicator variables for GxEScanR
cov_gxescan <- cov
for(t in unique(cov_gxescan$study_gxe)) {
  cov_gxescan[paste0(t)] <- ifelse(cov_gxescan$study_gxe==t,1,0)
}

# checks (what studies are included and how many)
check <- dplyr::select(cov_gxescan, unique(cov_gxescan$study_gxe)) %>% 
  summarise_all(sum) %>% t() %>%  as.data.frame(.) %>% 
  rownames_to_column()

# clean up , make Kentucky reference study_gxe (arbitrary)
cov_gxescan <- dplyr::select(cov_gxescan, -Kentucky, -study_gxe, -age_ref_imp, age_ref_imp)

# save files
saveRDS(cov,         file = '~/data/GxEScanR_PhenotypeFiles/FIGI_GwasSet_sex_age_pc3_studygxe_124702_GLM.rds', version = 2)
saveRDS(cov_gxescan, file = "~/data/GxEScanR_PhenotypeFiles/FIGI_GwasSet_sex_age_pc3_studygxe_124702.rds",     version = 2)


