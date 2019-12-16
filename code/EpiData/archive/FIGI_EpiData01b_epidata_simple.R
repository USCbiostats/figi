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

epi_varlist <- scan(file = "List of variables.txt", what = character())
samplefile_varlist_keep <- c("pooledcompassid", "vcfid", "V2", "gwas_set", "studyname", "platform", "study" , "drop", "drop_reason", "gxe", "study_gxe", "adenoma_adv", "adenoma", "crc", "age_ref_imp", "age_dxsel_imp")

#------ GWAS SET -------#
# must fix SMC_COSM special character
gwas <- All_Merged %>% 
  #filter(drop == 0) %>% 
  mutate(pooledcompassid =  gsub('<c5>', "", pooledcompassid, perl = T),
         pooledcompassid = gsub('Å', "", pooledcompassid, perl = T))

tmp <- do.call(rbind, lapply(list.files("./original", full.names = T), fread, select = epi_varlist, colClasses = list(character = epi_varlist), stringsAsFactors = F)) %>% 
  mutate(pooledcompassid =  gsub('<c5>', "", pooledcompassid, perl = T),
         pooledcompassid = gsub('Å', "", pooledcompassid, perl = T))

Epi <- inner_join(gwas[, samplefile_varlist_keep], tmp, by = 'pooledcompassid')

#------ Some Checks ------#
# do this after filter drop == 0
x <- anti_join(gwas, tmp, by = 'pooledcompassid') # should be zero (after filtering drop == 0)
# do this without filtering drop == 0
y <- anti_join(tmp, gwas[, c('pooledcompassid', 'drop')], by = 'pooledcompassid') # many likely dropped as result of QC
table(y$gecco_study)

# the numbers below probably change after Flora IBD drops
# 101  108  136  139  142 
# 25 1675  863  557    3
# CCFR (misc), dachs4, nshds, tcga, and the reach_ad (n=3)

# make sure there are no duplicate pooledcompassids after filtering drop == 0
x <- filter(Epi, drop == 0)
any(duplicated(x$pooledcompassid)) # good

#------ SAVE ------#
rm(list=setdiff(ls(), c("samplefile", "Epi", "All_Merged")))
save.image("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")


# Write out CSV data for Jim
write.csv(Epi, file = "~/FIGI_Epi_20181126.csv", quote = T, row.names = F)
test <- read.csv("~/FIGI_Epi_20181126.csv", header = T)

# 11/29/2018 - checks for jim

# NSAIDs
# asp_ref = Yes if aspirin=Yes OR nsaids=Yes; otherwise No if aspirin=No OR nsaids=No; Otherwise NA
library(expss)
x <- filter(Epi, drop == 0 & gxe == 1)
table(x$asp_ref, x$nsaids)
cro(x$asp_ref, x$nsaids)
cro(x$aspirin, x$nsaids)


# asp_ref by studyname
library(reshape)
y <- data.frame(table(x$study, x$nsaids)) %>% 
  reshape::cast(Var1 ~ Var2) %>% dplyr::rename(Miss = V1) %>% 
  dplyr::select(Var1, Yes, No, Miss)
write.csv(y, file = "~/FIGI_study_nsaids_tab.csv", quote = F, row.names = F)