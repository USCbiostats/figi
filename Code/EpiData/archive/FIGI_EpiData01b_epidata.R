#=======================================================
# EPI DATA
# 
# Incorporate E variables
# Process duplicates, problem samples
# Subset to gxe set
#=======================================================
library(data.table) 
library(tidyverse)
rm(list = ls())
setwd('~/git/DATA/FIGI_samplefile_epi-181120/')
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_181120.RData")


#------ Prelim Set ------#
# remove dups/IBD related samples
# remove drops
# remove special character (fix)
All_Merged_NoDups <- All_Merged %>%
	#filter(!vcfid %in% ibd_flora$drop_sample) %>% 
	mutate(pooledcompassid = as.character(pooledcompassid),
				 pooledcompassid = gsub('<c5>', "", pooledcompassid, perl = T),
				 pooledcompassid = gsub('Å', "", pooledcompassid, perl = T))
table(All_Merged_NoDups$drop)
any(duplicated(All_Merged_NoDups$vcfid))
any(duplicated(All_Merged_NoDups$pooledcompassid)) # do again below after remove drop == 1

# Latest UKB data from Yi has huge number of samples, probably sent by mistake
# merge with samplefile, rewrite CSV file into the ./original folder
# N = 27,594
# ukb <- fread("140ukb0_original.csv") %>%
# 	inner_join(All_Merged_NoDups[, c('pooledcompassid'), drop = F], by = 'pooledcompassid')
# write.csv(ukb, file = "./original/140ukb0.csv", quote = F, row.names = F)

# read original epidata
# samplefile_varlist_keep <- c("netcdfid", "V2", "vcfid", "pooledcompassid", "gwas_set", "studyname", "platform" , "compassid", "study" , "outc", "sex", "age_ref", "study_site", "race_self", "hispanic_self" , "cancer_site", "cancer_site_sum1",  "cancer_site_sum2",  "cancer_site_sum3",  "crc",  "adenoma_adv", "adenoma", "age_dxsel",  "drop", "drop_reason", "gxe", "study_gxe")
samplefile_varlist_keep <- c("pooledcompassid", "vcfid", "V2", "gwas_set", "studyname", "platform", "study" , "drop", "drop_reason", "gxe", "study_gxe", "adenoma_adv", "adenoma", "crc", "age_ref_imp", "age_dxsel_imp")

# epi_varlist <- c("pooledcompassid", "compassid", "gecco_study", "study_site", "outc", "age_ref", "sex")
epi_varlist <- scan(file = "List of variables.txt", what = character())

# Read all individual study data (CSV), rbind
tmp <- do.call(rbind, lapply(list.files("./original", full.names = T), fread, select = epi_varlist, colClasses = list(character = epi_varlist), stringsAsFactors = F)) %>% 	mutate(pooledcompassid = as.character(pooledcompassid),
				 pooledcompassid = gsub('<c5>', "", pooledcompassid, perl = T),
				 pooledcompassid = gsub('Å', "", pooledcompassid, perl = T))

# merge with samplefile after remove dropped samples
# Epi <- inner_join(All_Merged_NoDups[                                  , samplefile_varlist_keep, drop = F], tmp, by = 'pooledcompassid')
Epi <- inner_join(All_Merged_NoDups[which(All_Merged_NoDups$drop != 1), samplefile_varlist_keep, drop = F], tmp, by = 'pooledcompassid')

#------ Some Checks ------#
# Genotype YES, Epi NO - should be ZERO after removing "drop" samples
x <- anti_join(All_Merged_NoDups[                                  , , drop = F], tmp, by = 'pooledcompassid')
x <- anti_join(All_Merged_NoDups[which(All_Merged_NoDups$drop != 1), , drop = F], tmp, by = 'pooledcompassid')
# Genotype NO, Epi YES
# (i'm not too concerned since it's probably part of UW QC - can confirm with Yi when you send email about the reimputed samples)
y <- anti_join(tmp, All_Merged_NoDups[                                  , c('pooledcompassid', 'drop')], by = 'pooledcompassid')
#y <- anti_join(tmp, All_Merged_NoDups[which(All_Merged_NoDups$drop != 1), c('pooledcompassid', 'drop')], by = 'pooledcompassid')
table(y$gecco_study)

# the numbers below have probably changed after Flora IBD drops..
# 101  108  136  139  142 
# 25 1675  863  557    3
# CCFR (misc), dachs4, nshds, tcga, and the reach_ad (n=3)

# there should NOT be any duplicates
# after dropping gxe == 0, 3 pairs still remain, need to verify IBD results and ask Flora
# However, remember that for PCA calculation you want to do a GWAS + GxE set
# so --- drop all pairs of duplicates to be sure
any(duplicated(Epi$pooledcompassid)) 
dups <- Epi[duplicated(Epi$pooledcompassid),]
dups_check <- filter(Epi, pooledcompassid %in% dups$pooledcompassid)
# table(dups_check$gxe)
# dups_check <- filter(dups_check, gxe == 1)
# any(duplicated(dups_check$pooledcompassid))


# single pair remaining, drop from analysis
# 6150341003_R02C02_OmniExpress
# J1601777_OncoArray_Custom
Epi <- filter(Epi, !vcfid %in% dups_check$vcfid)


#------ SAVE ------#
rm(list=setdiff(ls(), c("samplefile", "Epi", "All_Merged")))
save.image("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_11052018.RData")


# tabulate outc with study
table(Epi$outc, Epi$study)


# write out data for jim

write.csv(Epi, file = "~/FIGI_Epi_20181022.csv", quote = T, row.names = F)
test <- read.csv("~/FIGI_Epi_20181022.csv", header = T)
