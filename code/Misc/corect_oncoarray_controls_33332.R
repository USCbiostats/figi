#=============================================================================#
# 2109/09/17
#
# output list of vcfid (CONTROLS ONLY) for corect_oncoarray
# this will be used to subset plink files to controls only
#
# plink format (IID, FID)
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")

# GWAS set, controls only
cov_gxe <- figi %>%
  filter(filename == "corect_oncoarray",
         drop == 0,
         outc == "Control") %>% 
  mutate(IID = vcfid, FID = vcfid)
  
table(figi_samplefile$filename)

write.table(cov_gxe[, c("IID", "FID")], file = "/media/work/tmp/corect_oncoarray_controls_12684.txt", quote = F, row.names = F)




# Jeroen uses oncoarray custom ... 
# x <- filter(figi_samplefile, filename == "corect_oncoarray")
# table(x$study_gxe)
# 
# 
# x <- filter(figi_samplefile, filename == "oncoarray_to_usc")
# table(x$study_gxe)
# 
# 
# x <- filter(figi_samplefile, platform == "OncoArray+Custom")
# table(x$study_gxe, useNA = 'ifany')
# table(x$filename, useNA = 'ifany')
# 
# table(figi_samplefile$platform)
