#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI Analysis 10/26/2018
#
# Create covariate tables for GxEScanR + EPACTS
# Updated PCs also
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(data.table)

rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_10182018.RData")

# sanity
any(duplicated(Epi$vcfid))
any(duplicated(Epi$pooledcompassid))
table(Epi$studyname)
table(Epi$V2)


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


# GxE Set ----
pc30k <- fread("~/Dropbox/code/FIGI_EpiData/FIGI_GxESet.eigenvec", skip = 1, 
							 col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
exclude_samples <- fread("~/Dropbox/code/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
exclude_studies <- c( "UKB", "MOFFITT", "NGCCS", "GALEON")





#''''''''''''''''''''''''''''''''''''''''''''''''''''
# 09/24/2018 GxEScanR Run (no study variable)
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
exclude_samples <- fread("~/Dropbox/code/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
exclude_studies <- c( "UKB", "MOFFITT", "NGCCS", "GALEON")



cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies,
				 !vcfid %in% exclude_samples$IID) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 outcome = ifelse(outc == "Case", 1, 0),
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname))

# platform indicator variable
gwasset <- unique(cov$gwas_set)[!unique(cov$gwas_set) %in% "corect_oncoarray"]
for(s in gwasset) {
	varname = paste0("plat_ind_", s)
	cov[[varname]] <- ifelse(cov$gwas_set == s, 1, 0)
}

cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), starts_with("plat_ind_"), sex )

age_pc_platform_sex <- cov_pcs[complete.cases(cov_pcs), ]


#---------------------------------------------------
# Full Covariate File (model selection @ GLM)
# N = 87554
#
# on second thought, full vars is messy, edit as needed
#---------------------------------------------------
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
				 studynamef = as.factor(studyname))
cov_pcs <- cov[complete.cases(cov), ] %>% 
	dplyr::select(varlist)

# clean up blank cells
empty_as_na <- function(x){
	if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
	ifelse(as.character(x)!="", x, NA)
}
cov_pcs_clean <- cov_pcs %>% mutate_all(funs(empty_as_na))

write.table(cov_pcs_clean, file = "~/FIGI_Covariates_10092018_FULL.txt", quote= F, row.names = F, sep = '\t')

# wtf <- fread("~/FIGI_Covariates_10092018_FULL.txt", stringsAsFactors = F, header = T)

#=====================================================

#---------------------------------------------------
# Drop studies with imbalanced sex distributions
# Drop UKB
#---------------------------------------------------
exclude_studies <- c("MOFFITT", "NGCCS", "GALEON", "UKB", "ATBC", "HPFS", "PHS", "NHS", "SELECT", "USC_HRT_CRC", "WHI")

cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies,
				 !vcfid %in% exclude_samples$IID) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 outcome = ifelse(outc == "Case", 1, 0),
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname))

# study/platform indicator variable (SEARCH is ref)
studylist <- unique(cov$studyname)[!unique(cov$studyname) %in% "SEARCH"]
for(s in studylist) {
	varname = paste0("study_plat_ind_", s)
	cov[[varname]] <- ifelse(cov$studyname == s, 1, 0)
}

cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), starts_with("study_plat_ind"), sex )
	#dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), sex )

age_PCs_STUDY_sex <- cov_pcs[complete.cases(cov_pcs), ]

write.table(age_PCs_STUDY_sex, file = "~/FIGI_Covariates_09272018_age_PCs_STUDY_sex.txt", quote= F, row.names = F, col.names = T, sep = '\t')
system("scp ~/FIGI_Covariates_09272018_age_PCs_STUDY_sex.txt hpcc:~ ")





#---------------------------------------------------
# Adjustment by platform.. (gwas_set)
#---------------------------------------------------
exclude_studies <- c("MOFFITT", "NGCCS", "GALEON", "UKB", "ATBC", "HPFS", "PHS", "NHS", "SELECT", "USC_HRT_CRC", "WHI")
exclude_studies <- c("MOFFITT", "NGCCS", "GALEON", "UKB", "ATBC", "HPFS", "PHS", "NHS", "SELECT", "USC_HRT_CRC", "WHI")

cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies,
				 !vcfid %in% exclude_samples$IID) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 outcome = ifelse(outc == "Case", 1, 0),
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname))

# platform indicator variable
gwasset <- unique(cov$gwas_set)[!unique(cov$gwas_set) %in% "corect_oncoarray"]
for(s in gwasset) {
	varname = paste0("plat_ind_", s)
	cov[[varname]] <- ifelse(cov$gwas_set == s, 1, 0)
}

cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), starts_with("plat_ind_"), sex )

age_pc_platform_sex <- cov_pcs[complete.cases(cov_pcs), ]

write.table(age_pc_platform_sex, file = "~/FIGI_Covariates_09272018_age_pc_platform_sex.txt", quote= F, row.names = F, col.names = T, sep = '\t')
system("scp ~/FIGI_Covariates_09272018_age_pc_platform_sex.txt hpcc:~ ")


