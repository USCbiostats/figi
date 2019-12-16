#==================================================================================================#
# FIGI Analysis 11/20/2018 Update v2.0
# FIGI Analysis 03/13/2019 Update v2.1
# FIGI Analysis 04/24/2019 Update v2.2
# FIGI Analysis 07/29/2019 Update v2.2 (small fix, accidentally excluded a few samples)
#
# Please refer to documentation Yi created
# this newest sample file cleans up some sample issues we identified previously, don't get confused
# - simply merge all vcf sample names w/ samplefile
#
# This samplefile excludes:
# - BWHS, COLON, DACHS_4, NSHDS, OSUMC (no approval)
# - OFCCR - lots of overlap with CCFR_3 CCFR_4
# - TCGA - need to obtain data manually but don't bother!
#   (marginal analyses only, no E variables available)
#
# Process Epi Data, save RData object, output sample name list for VCF update
# Output covariates for GxEScanR
#==================================================================================================#
library("tidyverse")
library("data.table")
rm(list = ls())
setwd('~/data/FIGI_samplefile_epi-190729/')

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
#-----------------------------------------------------------------------------#
samplefile <- fread("Samplefile_HRC1.1_epi_v2.1-20190306_usc.txt")

any(duplicated(samplefile$vcfid))
table(samplefile$gwas_set)
table(samplefile$studyname)
table(samplefile$study_gxe)

# ------ check that sample names between vcf and samplefile match ------ #
# (to avoid issues with GxEScanR)
# (reassembled vcf samples to be 100% sure, from HRC_VCF_SampleRename folder)
wrapp <- function(x) {
  fread(x, stringsAsFactors = F, sep = "\t", col.names = "vcfid", header = F) %>% 
    mutate(filename = x)
}

vcf_samplenames <- do.call(rbind, lapply(list.files(path = "vcf_sample_list", full.names = T, pattern = "sample_list"), wrapp)) %>% 
  mutate(filename = gsub("vcf_sample_list/sample_list_vcf_", "", filename),
         filename = gsub(".txt", "", filename))

check <- inner_join(samplefile, vcf_samplenames, by = "vcfid")
check <- anti_join(samplefile, vcf_samplenames, by = "vcfid")
check <- anti_join(vcf_samplenames, samplefile, by = "vcfid")

# 3 samples ccfr_omni with no epi data, no samplefile entry
# 8 samples from oncoarray_to_usc that aren't present in every single chromosome, dropped all 

figi_samplefile <- full_join(samplefile, vcf_samplenames, by = "vcfid") %>% 
  filter(!vcfid %in% check$vcfid)

rm(vcf_samplenames, check, samplefile, wrapp)

table(figi_samplefile$filename, useNA = 'ifany')

save.image("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_190729.RData")





# ------ old but important code ------ #
#==================================================
# OUTPUT TABLE TO UPDATE SAMPLE NAMES (VCF + PLINK)
# (VCFID)
# Don't drop "drop == 1" yet
#
# OUTPUT LIST OF SAMPLES BY GWAS_SET
#==================================================
rm(list = ls())

# still fine using 08132018 version
# but use newest one for corect_oncoarray_nonEUR_reimpute
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_10182018_HRC_FIX.RData")
table(All_Merged$V2)

# by gwas_set (use the V2 variable, see above)
batches <- unique(All_Merged$V2)
batches <- "corect_oncoarray_nonEUR"
# Output table to update samples names (VCF - bcftools)
for(batch in batches) {
	write.table(All_Merged[which(All_Merged$V2 == batch) ,c(1, 3)], file = paste0("~/FIGI_Sample_Rename_VCF_08162018_", batch, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}

for(batch in batches) {
  write.table(All_Merged[which(All_Merged$V2 == batch) ,c(1, 3)], file = paste0("~/FIGI_Sample_Rename_VCF_08162018_", batch, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}


# system("cp ~/FIGI_Sample_Rename_VCF_08162018_ccfr_1m_1mduo.txt FIGI_Sample_Rename_VCF_08162018_ccfr_1m_1mduo_reimpute.txt")
# system("cp ~/FIGI_Sample_Rename_VCF_08162018_corect_oncoarray_nonEUR.txt FIGI_Sample_Rename_VCF_08162018_corect_oncoarray_nonEUR_reimpute.txt")
# system("scp ~/FIGI_Sample_Rename_VCF_08162018_* hpcc:~")
# system("rm ~/FIGI_Sample_Rename_VCF_08162018_*")


# output all vcfid N = 145936
write.table(All_Merged[, c('vcfid', 'vcfid')], file = "~/FIGI_Sample_Rename_VCF_08162018_ALL.txt",  quote = F, row.names = F, col.names = F, sep = "\t")
system("scp ~/FIGI_Sample_Rename_VCF_08162018_ALL.txt hpcc:~")
system("rm FIGI_Sample_Rename_VCF_08162018_ALL.txt")

# output all vcfid after drops
write.table(All_Merged[which(All_Merged$drop == 0), c('vcfid', 'vcfid')], file = "~/FIGI_Sample_Rename_VCF_08162018_ALL_drops.txt",  quote = F, row.names = F, col.names = F, sep = "\t")
system("scp ~/FIGI_Sample_Rename_VCF_08162018_ALL_drops.txt hpcc:~")
system("rm FIGI_Sample_Rename_VCF_08162018_ALL_drops.txt")

# Output table to update samples names (PLINK)
# for(batch in batches) {
# 	write.table(tmp[which(tmp$gwas_set2 == batch) ,c(1, 1, 3, 3)], file = paste0("~/FIGI_Sample_Rename_PLINK_08132018_", batch, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")
# }
# system("scp ~/FIGI_Sample_Rename_PLINK_08132018_* hpcc:~")
# system("rm ~/FIGI_Sample_Rename_PLINK_08132018_*")

# Redo for corect_oncoarray_nonEUR_reimpute (10/26/2018) ----
rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_10182018_HRC_FIX.RData")

corect_oncoarray_nonEUR_reimpute <- filter(All_Merged, V2 == "corect_oncoarray_nonEUR") %>% 
  dplyr::select("netcdfid", "vcfid")

write.table(corect_oncoarray_nonEUR_reimpute, file = paste0("./working/FIGI_Sample_Rename_VCF_10262018_corect_oncoarray_nonEUR_reimpute.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

