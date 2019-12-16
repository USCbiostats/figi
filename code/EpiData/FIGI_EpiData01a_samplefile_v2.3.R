#==================================================================================================#
# FIGI Analysis 11/20/2018 Update v2.0
# FIGI Analysis 03/13/2019 Update v2.1
# FIGI Analysis 04/24/2019 Update v2.2
# FIGI Analysis 07/29/2019 Update v2.2 (small fix, accidentally excluded a few samples)
# FIGI Analysis 12/05/2019 Update v2.3
#
# Please refer to documentation Yi created
#
# Excludes the following studies:
# - BWHS, COLON, DACHS_4, NSHDS, OSUMC (no approval)
# - OFCCR - lots of overlap with CCFR_3 CCFR_4
# - TCGA - no exposure data available
#
# Process Epi Data, save RData object, output sample name list for VCF update
# Output covariates for GxEScanR
#==================================================================================================#
library("tidyverse")
library("data.table")
rm(list = ls())
setwd('~/data/FIGI_samplefile_epi-191205/')

#-----------------------------------------------------------------------------#
# DATA MATCHING ----
#
# a few samples present in vcf files but not in samplefile:
# 
# ccfr_omni
# 124015376_124015376012
# 110029094_110030029094
# 110049431_110030049431
#
# oncoarray_to_usc
# J1457396
# J8477920
# J5180634
# J6301643
# J8838797
# J8941867
# J6267972
# J7063624
#
#
# new drops for HRC v2.3
# 160000059843_CCFR (drop from 0 to 1)
#-----------------------------------------------------------------------------#
# samplefile <- fread("Samplefile_HRC1.1_epi_v2.3-20191205_usc.txt")
samplefile = read.table('Samplefile_HRC1.1_epi_v2.3-20191205_usc.txt', header=T, fileEncoding='latin1', stringsAsFactors= F, sep='\t', na.strings=c('',NA))

any(duplicated(samplefile$vcfid))
table(samplefile$gwas_set)
table(samplefile$gxe)
table(samplefile$studyname)
table(samplefile$study_gxe)


# ------ check that sample names between vcf and samplefile match ------ #
# (to avoid issues with GxEScanR)
# (reassembled vcf samples to be 100% sure, from HRC_VCF_SampleRename folder)
# (also adds variable called 'filename', which corresponds to imputation batch)
wrap <- function(x) {
  fread(x, stringsAsFactors = F, sep = "\t", col.names = "vcfid", header = F) %>% 
    mutate(filename = x)
}

vcf_samplenames <- do.call(rbind, lapply(list.files(path = "vcf_sample_list", full.names = T, pattern = "sample_list"), wrap)) %>% 
  mutate(filename = gsub("vcf_sample_list/sample_list_vcf_", "", filename),
         filename = gsub(".txt", "", filename))

# 3 samples ccfr_omni with no epi data, no samplefile entry
# 8 samples from oncoarray_to_usc that aren't present in every single chromosome, dropped all 
check <- inner_join(samplefile, vcf_samplenames, by = "vcfid") 
check <- anti_join(samplefile, vcf_samplenames, by = "vcfid") # 8186 in samplefile only 
check <- anti_join(vcf_samplenames, samplefile, by = "vcfid") # 11 in vcf files only - OK

figi_samplefile <- full_join(samplefile, vcf_samplenames, by = "vcfid") %>% 
  dplyr::filter(!vcfid %in% check$vcfid) # remove the 11 non-matches to avoid confusion

rm(vcf_samplenames, check, samplefile, wrap)

table(figi_samplefile$filename, useNA = 'ifany')

save.image("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_191205.RData")




