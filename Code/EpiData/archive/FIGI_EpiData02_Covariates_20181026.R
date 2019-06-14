#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI Analysis 10/26/2018
#
# Create covariate tables for GxEScanR + EPACTS
# Updated PCs also
#
#
# COVARIATE DATA.FRAME
# 
# Notes:
# use gxe set as determined in `gxe` variable 
# 	- gxe == 1 removes nonEUR studies
#		- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
# 	- gxe == 1 removes outc == "other"!
#
# studyname = adj. for study + platform
# remove case-only studies
# add principal components
# 	- remember these were calculated on whole set N = 141,362
#   - updated PCs, 2 sets available: GWAS and GxE
#
# Combine NFCCR_1 NFCCR_2
#		- NFCCR_1 only has cases, NFCCR_2 has cases/controls
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(data.table)

rm(list = ls())
# load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_10182018.RData")
# load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_10312018.RData")
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_11052018.RData")


# sanity
any(duplicated(Epi$vcfid))
any(duplicated(Epi$pooledcompassid))
table(Epi$studyname)
table(Epi$V2)

#------ GxE Set (N = 103,154) ------
pc30k <- fread("~/Dropbox/code/FIGI_EpiData/FIGI_GxESet_20181018.eigenvec", skip = 1, 
							 col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
exclude_samples <- fread("~/Dropbox/code/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
exclude_studies <- c( "MOFFITT", "NGCCS", "GALEON") # case-only studies.. 

cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies,
				 !vcfid %in% exclude_samples$IID) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0),
	       age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname),
				 studyname = as.factor(studyname)) %>% 
  dplyr::select(vcfid, outcome, age_ref, sex, study, studyname, paste0(rep("PC", 20), seq(1,20)))

saveRDS(cov, file = '~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_N_101526.rds')

write.table(cov, file = "~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_N101526.cov", quote = F, row.names = F, col.names = T, sep = '\t')



#------ GxE Set - GxEScanR (N = 103,154) ------
# for 11/05/2018 run, not adjusting for studyname, first 10 PCs only

pc30k <- fread("~/Dropbox/code/FIGI_EpiData/FIGI_GxESet_20181018.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
exclude_samples <- fread("~/Dropbox/code/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
exclude_studies <- c( "MOFFITT", "NGCCS", "GALEON", "UKB") # case-only studies + UKB

cov <- Epi %>%
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  filter(gxe == 1, 
         !study %in% exclude_studies,
         !vcfid %in% exclude_samples$IID) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref = as.numeric(age_ref),
         sex = ifelse(sex == "Female", 0, 1), 
         studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname),
         studyname = as.factor(studyname)) %>% 
  dplyr::select(vcfid, outcome, age_ref, paste0(rep("PC", 10), seq(1,10)), sex)

saveRDS(cov, file = '~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_GxEScanR_N86639_20181105.rds')



#------ GxE Set - GxEScanR (N = 102,799) ------

pc30k <- fread("~/Dropbox/code/FIGI_EpiData/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
exclude_samples <- fread("~/Dropbox/code/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
exclude_studies <- c( "MOFFITT", "NGCCS", "GALEON", "UKB") # case-only studies + UKB

cov <- Epi %>%
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  filter(gxe == 1, 
         !study %in% exclude_studies,
         !vcfid %in% exclude_samples$IID) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref = as.numeric(age_ref),
         sex = ifelse(sex == "Female", 0, 1), 
         heightcm = as.numeric(heightcm),
         studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname),
         studyname = as.factor(studyname)) %>% 
  dplyr::select(vcfid, outcome, age_ref, sex, paste0(rep("PC", 10), seq(1,10)), heightcm) %>% 
  filter(complete.cases(.))
  

# 79453

saveRDS(cov, file = '~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_GxEScanR_height_N86286_20181106.rds')
saveRDS(cov, file = '~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_GxEScanR_height_N79453_20181106.rds')




# platform indicator variable (eventually include for GxEScanR)
# gwasset <- unique(cov$gwas_set)[!unique(cov$gwas_set) %in% "corect_oncoarray"]
# for(s in gwasset) {
# 	varname = paste0("plat_ind_", s)
# 	cov[[varname]] <- ifelse(cov$gwas_set == s, 1, 0)
# }

cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), starts_with("plat_ind_"), sex )

age_pc_platform_sex <- cov_pcs[complete.cases(cov_pcs), ]


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Full Covariate File (model selection @ GLM)
# N = 87554
#
# on second thought, full vars is messy, edit as needed
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
exclude_studies <- c( "UKB", "MOFFITT", "NGCCS", "GALEON")
cat(unique(sort(names(cov))), sep = '", "')
varlist <- c("age_ref", "alcohol", "alcohol_ref", "alcoholc", "asp_ref", "bmi", "educ", "energytot", "exercise", "famhx1", "gxe", "nsaids", "nsaids_ever", "nsaidsyrs", "outc", "outcome", "PC1", "PC10", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "platform", "sex", "smoke", "study", "study_gxe", "study_site", "studyname", "V2", "vcfid", "studynamef")

cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies,
				 !vcfid %in% exclude_samples$IID) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 outcome = ifelse(outc == "Case", 1, 0),
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname),
				 studyname = as.factor(studyname))
cov_pcs <- cov[complete.cases(cov), ] %>% 
	dplyr::select(varlist)


# clean up blank cells
empty_as_na <- function(x){
	if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
	ifelse(as.character(x)!="", x, NA)
}
cov_pcs_clean <- cov_pcs %>% mutate_all(funs(empty_as_na))

write.table(cov_pcs_clean, file = "~/FIGI_Covariates_10092018_FULL.txt", quote= F, row.names = F, sep = '\t')

