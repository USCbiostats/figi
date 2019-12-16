#=============================================================================#
# Estimate Rsq using alt allele probability variance
# (see BinaryDosage fork on KimAE account -- master)
# 
# use GWASSet, but EXCLUDE non-EUR oncoarray batch samples
# total N = 131679
#=============================================================================#
library(tidyverse)
rm(list = ls())

load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
out <- filter(figi, filename != "corect_oncoarray_nonEUR_reimpute")[, 'vcfid']

saveRDS(out, file = "/home/rak/data/Rsq_Estimate/FIGI_GWASSet_vcfid_N_131679.rds", version = 2)






# x <- readRDS("/Users/mak/data/FIGI_EpiData_rdata/FIGI_chr22.rds")
# 
# samples <- x[[9]]
# head(samples)
# 
# 
# wtf <- inner_join(figi, samples, by = c("vcfid" = "SID"))
# 
# wtf <- anti_join(figi, samples, by = c("vcfid" = "SID"))
# table(wtf$filename)
