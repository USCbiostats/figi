#==============================================================================#
# FIGI Analysis 09/04/2018 **
# FIGI Analysis 10/17/2018 Update
# FIGI Analysis 10/31/2018 Update - more drops
# FIGI Analysis 11/20/2018 Update - more drops and other alterations
#
# This samplefile excludes:
# - BWHS, COLON, DACHS_4, NSHDS, OSUMC (no approval)
# - OFCCR - lots of overlap with CCFR_3 CCFR_4
# - TCGA - need to obtain data manually but don't bother!
#   (marginal analyses only, no E variables available)
#
# Process Epi Data, save RData object, output sample name list for VCF update
# Output covariates for GxEScanR
#==============================================================================#
library("tidyverse")
library("data.table")

setwd('~/data/FIGI_samplefile_epi-181120/')
rm(list = ls())

exclude <- c("BWHS", "COLON", "DACHS_4", "NSHDS", "OSUMC", "TCGA", "OFCCR")
#==================================================#
# DATA MATCHING ----
#==================================================#
#------ verify sample names match, output sample names for bcftools ------#
samplefile <- fread("Samplefile_HRC1.1_epi_v2-20181126_usc.txt") %>% 
  filter(!studyname %in% exclude)
table(samplefile$gwas_set)
table(samplefile$studyname)


# VCF FILE SAMPLE LIST
# USC files have an extra platform column from previous code. Create same columns for GECCO to match
axiom_acs <- fread("./samplenames_vcf_ak/axiom_acs_aus_nf_samplelist.txt", header = F)
axiom_mecc <- fread("./samplenames_vcf_ak/axiom_mecc_cfr_ky_samplelist.txt", header = F)
corect_oncoarray <- fread("./samplenames_vcf_ak/corect_oncoarray_samplelist.txt", header = F)
mecc <- fread("./samplenames_vcf_ak/mecc_samplelist.txt", header = F) %>%	mutate(V1 = as.character(V1))

ccfr_1m <- fread("./samplenames_vcf_ak/ccfr_1m_1mduo_samplelist.txt", header = F) %>%
  separate(V1, into = c("V4", "trash"), sep = "_", remove = F) %>%
  dplyr::rename(tmp = V1, V1 = V4)  %>% dplyr::select(-trash)
ccfr_omni <- fread("./samplenames_vcf_ak/ccfr_omni_samplelist.txt", header = F) %>%
  separate(V1, into = c("trash", "V4"), sep = "_", remove = F) %>%
  dplyr::rename(tmp = V1, V1 = V4)  %>% dplyr::select(-trash) # FYI no dups between 1m & omni

corsa_axiom <- fread("./samplenames_vcf_ak/corsa_axiom.txt", header = F) %>%
  mutate(V2 = "corsa_axiom")
cytosnp_comb <- fread("./samplenames_vcf_ak/cytosnp_comb.txt", header = F) %>%
  mutate(V2 = "cytosnp_comb")
initial_comb_datasets <- fread("./samplenames_vcf_ak/initial_comb_dataset.txt", header = F) %>%
  mutate(V2 = "initial_comb_datasets")
# initial_comb <- fread("./samplenames_vcf_ak/initial_comb.txt", header = F) %>%
# 	mutate(V2 = "initial_comb") # overlaps completely with initial_comb_dataset
newfoundland_omniquad <- fread("./samplenames_vcf_ak/newfoundland_omniquad.txt", header = F) %>%
  mutate(V2 = "newfoundland_omniquad")
omni_comb <- fread("./samplenames_vcf_ak/omni_comb.txt", header = F) %>%
  mutate(V2 = "omni_comb")
omni_express_cone <- fread("./samplenames_vcf_ak/omni_express_cone.txt", header = F) %>%
  mutate(V2 = "omniexpress_exomechip")
oncoarray_to_usc <- fread("./samplenames_vcf_ak/oncoarray_to_usc.txt", header = F) %>%
  mutate(V2 = "oncoarray_to_usc")
plco3 <- fread("./samplenames_vcf_ak/plco3.txt", header = F) %>%
  mutate(V2 = "plco_3")
# plco4 <- fread("./samplenames_vcf_ak/plco4.txt", header = F) %>%
# 	mutate(V2 = "plco4") # overlaps completely with oncoarray_to_usc samples
ukbiobank <- fread("./samplenames_vcf_ak/ukbiobank_samplelist.txt", header = F) %>%
  mutate(V1 = as.character(V1), V2 = 'ukbiobank')

# new files since
corect_oncoarray_nonEUR <- fread("~/Dropbox/code/corect_oncoarray_reimpute/corect_oncoarray_release_v3_impute_nonEUR.fam", header = F) %>%	mutate(V1 = as.character(V1), V2 = 'corect_oncoarray_nonEUR') %>% select(V1, V2)
dachs3 <- fread("./samplenames_vcf_ak/dachs3.txt", header = F) %>%
  mutate(V1 = as.character(V1), V2 = 'dachs3')

reach <- fread("./samplenames_vcf_ak/reach.txt", header = F) %>%
  separate(V1, into = c("V4", "trash"), sep = "_", remove = F) %>%
  dplyr::rename(tmp = V1, V1 = V4)  %>% dplyr::select(-trash) %>%
  mutate(V2 = "reach")


#------ merges ------#
figi_axiom_acs <- inner_join(axiom_acs, samplefile, by = c("V1" = "netcdfid")); rm(axiom_acs)
figi_axiom_mecc <- inner_join(axiom_mecc, samplefile, by = c("V1" = "netcdfid")) %>% filter(gwas_set == "axiom_mecc_cfr_ky") ; rm(axiom_mecc)
figi_corect_oncoarray <- inner_join(corect_oncoarray, samplefile, by = c("V1" = "netcdfid")); rm(corect_oncoarray)
figi_mecc <- inner_join(mecc, samplefile, by = c("V1" = "netcdfid")) %>% filter(gwas_set == "mecc"); rm(mecc)
figi_corsa_axiom <- inner_join(corsa_axiom, samplefile, by = c("V1" = "netcdfid")); rm(corsa_axiom)
figi_cytosnp_comb <- inner_join(cytosnp_comb, samplefile, by = c("V1" = "netcdfid")); rm(cytosnp_comb)
figi_initial_comb_datasets <- inner_join(initial_comb_datasets, samplefile[which(samplefile$gwas_set == "initial_comb"),], by = c("V1" = "netcdfid")); rm(initial_comb_datasets)
figi_newfoundland_omniquad <- inner_join(newfoundland_omniquad, samplefile, by = c("V1" = "netcdfid")); rm(newfoundland_omniquad)
figi_omni_comb <- inner_join(omni_comb, samplefile, by = c("V1" = "netcdfid")); rm(omni_comb)
figi_omni_express_cone <- inner_join(omni_express_cone, samplefile, by = c("V1" = "netcdfid")); rm(omni_express_cone)
figi_plco3 <- inner_join(plco3, samplefile, by = c("V1" = "netcdfid")); rm(plco3)
figi_ukbiobank <- inner_join(ukbiobank, samplefile, by = c("V1" = "netcdfid")); rm(ukbiobank)
figi_dachs3 <- inner_join(dachs3, samplefile, by = c("V1" = "netcdfid")); rm(dachs3)

### REMOVE REACH_AD, adding with new PKG from keith (REACH SPECIFICALLY)
# 20912 vs. 20904
figi_oncoarray_to_usc <- inner_join(oncoarray_to_usc, samplefile, by = c("V1" = "netcdfid")) %>% filter(studyname != "REACH_AD"); rm(oncoarray_to_usc) 

figi_ccfr_1m <- inner_join(ccfr_1m, samplefile, by = c("V1" = "netcdfid")) %>% dplyr::select(-V1) %>% dplyr::rename(V1 = tmp); rm(ccfr_1m) #  replace netcdfid value on samplefile to reflect VCF

### ccfr_omni (N = 3 not matched) ** bring up in discussion in detail ###
# 1600 vs. 1597
figi_ccfr_omni <- inner_join(ccfr_omni, samplefile[which(samplefile$gwas_set == "101ccfr_usc2"),], by = c("V1" = "netcdfid")) %>%	dplyr::select(-V1) %>% dplyr::rename(V1 = tmp); rm(ccfr_omni) #  replace netcdfid value on samplefile to reflect VCF

### corect_oncoarray_nonEUR (N = 2 not matched) ** bring up in discussion in detail ###
# 6901 vs. 6899
figi_corect_oncoarray_nonEUR <- inner_join(corect_oncoarray_nonEUR, samplefile, by = c("V1" = "netcdfid")); rm(corect_oncoarray_nonEUR)

figi_reach <- inner_join(reach, samplefile, by = c("V1" = "netcdfid")) %>% dplyr::select(-V1) %>% dplyr::rename(V1 = tmp); rm(reach) #  replace netcdfid value on samplefile to reflect VCF



#------ Concatenate all ------#
All_Merged <- do.call(rbind, mget(apropos("figi"))) %>% dplyr::rename(netcdfid = V1)
Anti_samplefile <- anti_join(samplefile, All_Merged, by = 'vcfid')
table(Anti_samplefile$studyname) # drop == 1

any(duplicated(All_Merged$vcfid))
any(duplicated(All_Merged$pooledcompassid)) # drop == 1

rm(list=setdiff(ls(), c("samplefile", "All_Merged")))
save.image("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_181120.RData")


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