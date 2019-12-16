#=======================================================
# EPI DATA V2.1
# 03/13/2019
# 
# Incorporate E variables
# Process duplicates, problem samples
# Subset to gxe set
#=======================================================
library(data.table) 
library(tidyverse)
rm(list = ls())

# start with samplefile, filter by drops/gxe
load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_190313.RData")

df <- figi_samplefile %>% 
  filter(drop == 0) %>% 
any(duplicated(df$pooledcompassid))

# ------ Preliminary Set ------ #
# E available, G not available = many (no genetic data)
# G available, E not available = OSUMC only (should they be sending these data?)

wrapp <- function(x) {
  exclude <- c("compassid", "outc", "study_site", "age_ref", "sex", "race_self", "hispanic_self", "cancer_site", "age_dxsel", "cancer_site_sum1", "cancer_site_sum2", "cancer_site_sum3")
  epi_varlist <- scan(file = "~/data/FIGI_samplefile_epi-190313/list_of_variables.txt", what = character())
  epi_varlist <- epi_varlist[!epi_varlist %in% exclude]
  fread(x, select = epi_varlist, colClasses = list(character = epi_varlist), stringsAsFactors = F) %>%
    mutate(pooledcompassid = as.character(pooledcompassid))
}

tmp <- do.call(rbind, lapply(list.files("~/data/FIGI_samplefile_epi-190313/epi", full.names = T), wrapp)) 

figi <- inner_join(tmp, df, by = 'pooledcompassid') %>% 
  mutate(pooledcompassid = as.character(pooledcompassid),
         pooledcompassid = gsub('<c5>', "", pooledcompassid, perl = T),
         pooledcompassid = gsub('Å', "", pooledcompassid, perl = T))

any(duplicated(figi$vcfid))
any(duplicated(figi$pooledcompassid))

rm(list=setdiff(ls(), c("figi_samplefile", "figi")))
save.image("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_20190313.RData")

