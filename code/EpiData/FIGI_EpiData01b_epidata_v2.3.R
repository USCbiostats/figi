#=============================================================================#
# EPI DATA V2.3
# 12/05/2019
# 
# Incorporate E variables
# Process duplicates, problem samples
# Subset to gxe set
#
#
# ------ NOTES ------ #
#
# Changes for HRC v2.3 vs v2.1
#
# CCFR_4
# 160000059843_CCFR - dropped
#
# HispanicCCS
# 87 samples not matching (pooledcompassid contains "HispanicCCS" instead of "HCCS")
# (check with Yi to make sure this is deliberate)
#
# Colo23
# 211 samples. Need to edit pooledcompassid, update csv file
#
# PMH-CCFR
# epi file not provided (yet?)
#
# --- UPDATE --- #
# Yi sent updated samplefile (12/10/2019), no longer an issue. 
#
#=============================================================================#
library(data.table) 
library(tidyverse)
setwd("~/data/FIGI_samplefile_epi-191205/")
#-----------------------------------------------------------------------------#
# edit colo23 csv file pooledcompassid (only run once)
#-----------------------------------------------------------------------------#

# colo <- read.csv('epi/109colo230.csv', header=T, fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
#   mutate(pooledcompassid = gsub("&", "", pooledcompassid))
# write.csv(colo, file = "epi/109colo230.csv", quote = F, row.names = F)


#-----------------------------------------------------------------------------#
# incorporate epi data
#-----------------------------------------------------------------------------#
rm(list = ls())

# start with samplefile, filter by drops/gxe
# load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_190729.RData")
load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_191205.RData")

df <- figi_samplefile %>% 
  filter(drop == 0)
any(duplicated(df$pooledcompassid))


# when reading epi data, don't include variables names found in samplefile (except for pooledcompassid)
list_of_variables <- read.table("~/data/FIGI_samplefile_epi-191205/list_of_variables.txt", stringsAsFactors = F) %>% 
  filter(!V1 %in% names(figi_samplefile)[!names(figi_samplefile) %in% 'pooledcompassid']) %>% 
  .[,1]

wrap <- function(x) {
  read.csv(x, colClasses = 'character', fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
    dplyr::select_if(names(.) %in% list_of_variables)
}
tmp <- lapply(list.files("~/data/FIGI_samplefile_epi-191205/epi", full.names = T), wrap)


# find intersection of column names for all epi datasets
# run wrap function again with this set of variable names
list_of_variables_epi <- Reduce(intersect, lapply(tmp, names))

wrap <- function(x) {
  read.csv(x, colClasses = 'character', fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
    dplyr::select_if(names(.) %in% list_of_variables_epi)
}
tmp <- do.call(rbind, lapply(list.files("~/data/FIGI_samplefile_epi-191205/epi", full.names = T), wrap))


# final figi set
figi <- inner_join(df, tmp, by = 'pooledcompassid')
  

any(duplicated(figi$vcfid))
any(duplicated(figi$pooledcompassid))



# Epi data, no genotype information
mismatches <- anti_join(tmp, df, 'pooledcompassid')
table(mismatches$gecco_study)

mismatches <- anti_join(df, tmp, 'pooledcompassid')
table(mismatches$study)

rm(list=setdiff(ls(), c("figi_samplefile", "figi")))
save.image("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")



#-----------------------------------------------------------------------------#
# Sample counts for Yi
#-----------------------------------------------------------------------------#
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")

# gwas set
counts <- data.frame(table(figi$study, figi$outc)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  rename(study = Var1)
write.csv(counts, file = "~/Dropbox/HRC_v2.3_GWAS_Set_counts.csv", quote = F, row.names = F)

# gwas set
gxe <- filter(figi, gxe == 1)
counts <- data.frame(table(gxe$study, gxe$outc)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  rename(study = Var1)
write.csv(counts, file = "~/Dropbox/HRC_v2.3_GxE_Set_counts.csv", quote = F, row.names = F)





#-----------------------------------------------------------------------------#
# Some error checking, accounting for mismatches
# 
#-----------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
library(figifs)
rm(list = ls())

# compare july and december files.. 
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
july <- figi
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")
december <- figi

# non-union for gwas set
# a single sample changed (CCFR_4) - OK!
non_matches <- anti_join(july, december, by = 'vcfid')
table(non_matches$drop)
table(non_matches$study_gxe)


# gxe set
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
july <- filter(figi, gxe == 1)
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")
december <- filter(figi, gxe == 1)

non_matches <- anti_join(july, december, by = 'vcfid')



