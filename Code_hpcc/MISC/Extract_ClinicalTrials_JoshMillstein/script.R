#=============================================================================#
#
# Extract IDs for Joshua Millstein (clinical trials)
# (from VCF files with fixed sample names)
#
#=============================================================================#
library(tidyverse)
library(data.table)
rm(list = ls())

setwd("~/Dropbox/FIGI/Code_hpcc/MISC/Extract_ClinicalTrials_JoshMillstein/")
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190313.RData")

# id list from joshua email. Contains netcdfids
idlist <- fread("IDs_FIRE3_Mavericc_Tribe.txt", header = F)
table(idlist$V3)


# merge, get vcfid, output for vcftools
x <- inner_join(figi_samplefile, idlist, by = c("netcdfid" = "V1"))
write.table(x$vcfid, file = "IDs_FIRE3_Mavericc_Tribe_vcftools.txt", quote = F, row.names = F, col.names = F)

table(x$gwas_set)
