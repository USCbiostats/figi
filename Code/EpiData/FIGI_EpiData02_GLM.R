#-----------------------------------------------------------
# FIGI Test Run ALL FILES... 
#
# Using only samplefile covariates (no genetic data yet)
#-----------------------------------------------------------

library(tidyverse)
library(data.table)

rm(list = ls())
load("~/git/FIGI_EpiData/FIGI_SampleFile_Summary_08132018.RData")

any(duplicated(All_Merged$vcfid))
table(All_Merged$studyname)
table(All_Merged$V2)


#------ Study variable ------#
# first, exclude nonEUR studies... 
# Then deal with case only studies
# for now, combining based on geographical
# for now, keep gxe == 0 (for study specific runs...)

df <- data.frame(table(All_Merged$study, All_Merged$study_site))


# some easy exclusions.. (asian)
exclude <- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
cov <- All_Merged %>%
	filter(drop == 0, 
				 outc != "Other", 
				 !is.na(outc), 
				 #gxe == 1, 
				 !studyname %in% exclude) %>% 
	mutate(study_adj = study) %>% 
	mutate(study_adj = ifelse(study == "GALEON" | 
															study == "MOFFITT" |
															study == "NGCCS" |
															study == "CRCGEN" , "GALEON_MOFFITT_NGCCS_CRCGEN", study_adj)) %>% 
	mutate(study_adj = ifelse(study == "UKB" | study == "SEARCH", "SEARCH_UKB", study_adj)) %>% 
	mutate(study_adj = ifelse(study == "SMC_COSM" | study == "SMS", "SMS_COSM_SMS_AD", study_adj))

table(cov$outc)
table(cov$studyname)
table(cov$study)
table(cov$outc, cov$study)
table(cov$outc, cov$study_adj)
table(cov$outc, cov$V2)

# Add Principal Components
#system("scp hpcc:/auto/pmd-02/figi/PCA/ALL_Merged_PCA_drops.eigenvec ~/")
pc30k_header = c("FID", "IID", paste0(rep("PC", 10), seq(1,10)))
pc30k <- fread("~/ALL_Merged_PCA_drops.eigenvec", skip = 1, col.names = pc30k_header)

#------------------------------------------------------------
# Remove samples with missing covariates 
# (note PCs were calculated on whole set N = 141,362)
#
# for study adjustment, use SEARCH_UKB as ref
# (this is still in progress, finish after some discussion..
#
#------------------------------------------------------------
cov_pcs <- inner_join(cov, pc30k, by = c('vcfid' = 'IID')) %>%
	mutate(age_ref = as.numeric(age_ref), 
				 sex = ifelse(sex == "Male", 0, 1),
				 outcome = ifelse(outc == "Case", 1, 0))



# clean and organize
cov_pcs <- cov_pcs %>% 
	filter(gxe == 1,
				 studyname != "UKB_1", 
				 studyname != "MOFFITT", 
				 studyname != "NFCCR_1",
				 studyname != "NGCCS", 
				 studyname != "GALEON",
				 study != "SMS") %>% 
	mutate(sex = factor(sex, labels = c("Male", "Female")))

table(cov_pcs$outcome, cov_pcs$studyname)
table(cov_pcs$studyname)
table(cov_pcs$study)
table(cov_pcs$outcome, useNA = 'always')
table(cov_pcs$crc, useNA = 'always')

table(cov_pcs$race_self)
y <- filter(cov_pcs, race_self == "Asian")
x <- glm(outcome ~ age_ref + PC1 + PC2 + PC3 + PC4 + sex, data = cov_pcs, family = binomial)

summary(x)
