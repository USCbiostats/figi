#-----------------------------------------------------------
# FIGI Analysis 09/04/2018 **
#
# Create covariate tables for GxEScanR + EPACTS etc
#-----------------------------------------------------------
library(tidyverse)
library(data.table)

rm(list = ls())
load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_09242018.RData")

# sanity
any(duplicated(All_Merged$vcfid))
table(All_Merged$studyname)
table(All_Merged$V2)

table(Epi$studyname)
table(Epi$V2)


#----------------------------------------------------------------------------------------
# COVARIATE Table
# 
# Notes:
# use gxe set as determined in `gxe` variable 
# 	- gxe == 1 removes nonEUR studies
# 	- gxe == 1 removes outc == "other"
# use studyname as study/platform adjustment variable (use oncoarray as ref for now)
# remove case-only studies
# add principal components
# temporarily combined NFCCR_1 NFCCR_2 studynames.. 
#----------------------------------------------------------------------------------------

# (note PCs were calculated on whole set N = 141,362)
#system("scp hpcc:/auto/pmd-02/figi/PCA/ALL_Merged_PCA_drops.eigenvec ~/")
pc30k <- fread("~/git/FIGI_EpiData/ALL_Merged_PCA_drops.eigenvec", skip = 1, 
							 col.names = c("FID", "IID", paste0(rep("PC", 10), seq(1,10))))

#exclude_studies <- c("HispanicCCS", "Taiwan", "SMHS", "SWHS", "HCES-CRC", "PURIFICAR")
exclude_studies <- c("MOFFITT", "NGCCS", "GALEON")
	

cov <- Epi %>%
	inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
	filter(gxe == 1, 
				 !study %in% exclude_studies) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1), 
				 outcome = ifelse(outc == "Case", 1, 0),
				 studyname = ifelse(grepl("NFCCR", studyname), "NFCCR", studyname))

table(cov$outc, cov$study)
table(cov$outc, cov$studyname)

# STUDY INDICATOR VARIABLE (UKB as reference)
studylist <- unique(cov$studyname)[!unique(cov$studyname) %in% "UKB_1"]
for(s in studylist) {
	varname = paste0("study_plat_ind_", s)
	cov[[varname]] <- ifelse(cov$studyname == s, 1, 0)
}

cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), starts_with("study_plat_ind"), sex )

cov_pcs_noNA <- cov_pcs[complete.cases(cov_pcs), ]
#-------------------------------------------
# write out file, send to HPC
write.table(cov_pcs_noNA, file = "~/FIGI_Covariates_09042018.txt", quote= F, row.names = F, col.names = T, sep = '\t')
system("scp ~/FIGI_Covariates_09042018.txt hpcc:~ ")
system("scp ~/FIGI_Covariates_09042018.txt rak:~ ")



#------------------------------------------
# remove 'problem' samples from oncoarray_to_usc
# (some IDs not present in all chromosome files
oncoarray_to_usc <- fread("~/git/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt")
cov_pcs <- cov %>% 
	filter(studyname != "UKB_1") %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), sex )
cov_pcs_noNA <- cov_pcs[complete.cases(cov_pcs), ]
cov_pcs_noNA_oncoarrayfix <- filter(cov_pcs_noNA, ! vcfid %in% oncoarray_to_usc$IID) # 7 samples excluded as a result.. (the eighth could have been removed b/c of covars)

write.table(cov_pcs_noNA_oncoarrayfix, file = "~/FIGI_Covariates_NOSTUDY_09242018.txt", quote= F, row.names = F, col.names = T, sep = '\t')
# system("scp ~/FIGI_Covariates_NOSTUDY_09212018.txt hpcc:~ ")
# system("scp ~/FIGI_Covariates_NOSTUDY_09212018.txt rak:~ ")
system("scp ~/FIGI_Covariates_NOSTUDY_09242018.txt hpcc:~ ")
system("scp ~/FIGI_Covariates_NOSTUDY_09242018.txt rak:~ ")

# WITHOUT STUDY INDICATOR VARIABLES (just why not)
cov_pcs <- cov %>% 
	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), sex )
cov_pcs_noNA <- cov_pcs[complete.cases(cov_pcs), ]
write.table(cov_pcs_noNA, file = "~/FIGI_Covariates_NOSTUDY_09042018.txt", quote= F, row.names = F, col.names = T, sep = '\t')
system("scp ~/FIGI_Covariates_NOSTUDY_09042018.txt hpcc:~ ")
system("scp ~/FIGI_Covariates_NOSTUDY_09042018.txt rak:~ ")



# Fixing issue with oncoarray_to_usc

#--------------------------------------------------------
# write file - EPACTS
# cov_pcs_noNA_EPACTS <- cov_pcs_noNA %>%
# 	mutate(FAM_ID = vcfid, 
# 				 IND_ID = vcfid, 
# 				 FAT_ID = 0, 
# 				 MOT_ID = 0, 
# 				 SEX = sex, 
# 				 DISEASE = outcome) %>%
# 	dplyr::select(FAM_ID, IND_ID, FAT_ID, MOT_ID, SEX, DISEASE, age_ref, paste0('PC', seq(1,4)))
# 
# 
# cat("#", file = "~/FIGI_Covariates_08272018_EPACTS.ped")
# write.table(cov_pcs_noNA_EPACTS, file = "~/FIGI_Covariates_08272018_EPACTS.ped", quote= F, row.names = F, col.names = T, sep = '\t', append = T)
# system("scp ~/FIGI_Covariates_08272018_EPACTS.ped hpcc:~ ")
# system("scp ~/FIGI_Covariates_08272018_EPACTS.ped rak:~ ")



########################## OLD

#=====================================================
# PROCESS AND OUTPUT COVARIATE INFORMATION for GXESCAN
#=====================================================
# rm(list = ls())
# load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_09042018.RData")


# 
# # OUTPUT covariate file to test GxEScanR
# # Need to add PCs - do it two ways, pc20k and pc100k as a test
# system("scp hpcc:/staging/dvc/andreeki/TestRun_100KSNPS/*.eigenvec ~/work/FIGI_TestRun_85k")
# system("scp hpcc:/staging/dvc/andreeki/TestRun_20KSNPS/*.eigenvec ~/work/FIGI_TestRun_85k")
# 
# pc20k <- fread("~/work/FIGI_TestRun_85k/ALL_Merged_20KSNPS_KGP.eigenvec", stringsAsFactors = F, col.names = c("FID", "IID", paste0(rep("PC", 10), seq(1,10))))
# pc100k <- fread("~/work/FIGI_TestRun_85k/ALL_Merged_100KSNPS_KGP.eigenvec", stringsAsFactors = F, col.names = c("FID", "IID", paste0(rep("PC", 10), seq(1,10))))
# 	
# 
# ALL_Merged_pc20k_cov <- inner_join(All_Merged, pc20k, by = c('vcfid' = 'IID')) %>%
# 	mutate(sex = ifelse(sex == "Female", 0, 1),
# 				 outcome = ifelse(outc == "Case", 1, 0)) %>%
# 	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,10)), sex)
# 
# ALL_Merged_pc100k_cov <- inner_join(All_Merged, pc100k, by = c('vcfid' = 'IID')) %>%
# 	mutate(sex = ifelse(sex == "Female", 0, 1),
# 				 outcome = ifelse(outc == "Case", 1, 0)) %>%
# 	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,10)), sex)
# 
# 
# write.table(ALL_Merged_pc20k_cov, file = "~/Desktop/ALL_Merged_pc20k_cov.txt", quote = F, row.names = F, col.names = T)
# write.table(ALL_Merged_pc100k_cov, file = "~/Desktop/ALL_Merged_pc100k_cov.txt", quote = F, row.names = F, col.names = T)



#----------------------------------------
## Covariate file for the Chr8 test run
## exclude UKBIOBANK and REACH
# rm(list = ls())
# load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_08132018.RData")
# 
# cov <- mutate(FIGI_ImpGeno_EpiData,
# 							gwas_set2 = ifelse(gwas_set == "101ccfr_usc", "ccfr_1m_1mduo",
# 																 ifelse(gwas_set == "101ccfr_usc2", "ccfr_omni",
# 																 			 ifelse(gwas_set == "initial_comb", "initial_comb_datasets",
# 																 			 			 ifelse(gwas_set == "OmniExpress+ExomeChip", "omniexpress_exomechip", 			 			 			 			              ifelse(gwas_set == "OncoArray+Custom", "oncoarray_to_usc",
# 																 			 			 																																																				 ifelse(gwas_set == "UKB_28K_HRC1.1", "ukbiobank", gwas_set)))))),
# 							gwas_set2 = ifelse(studyname == "DACHS_3", "dachs3",
# 																 ifelse(studyname == "REACH_AD", "reach", gwas_set2))) %>%
# 	filter(V2 != "reach",
# 				 V2 != "ukbiobank",
# 				 V2 != "corect_oncoarray_nonEUR")
# 
# table(cov$gwas_set2)
# table(cov$studyname)
# table(cov$V2)
# 
# # output list of sample IDs to subset merged plink file prior PCA calculation
# write.table(cov[, c('vcfid', 'vcfid')], file = "~/FIGI_chr8_TestRun_Samplelist.txt", quote = F, col.names = F, row.names = F)
# system("scp ~/FIGI_chr8_TestRun_Samplelist.txt hpcc:/staging/dvc/andreeki/PCA/") 
# system("rm ~/FIGI_chr8_TestRun_Samplelist.txt")
# 
# 
# # read in PCs, output covariate table!
# eigenval <- c(551.062, 296.726, 183.206, 69.4328, 52.918, 50.8268, 50.0444, 43.0176, 37.8236, 33.6117)
# eigenval_cumsum <- 100*cumsum(eigenval)/sum(eigenval) # this will be helpful for RMarkdown
# 
# system("scp hpcc:/staging/dvc/andreeki/PCA/ALL_Merged_PCA.eigenvec ~/")
# 
# pc30k_header = c("FID", "IID", paste0(rep("PC", 10), seq(1,10)))
# pc30k <- fread("~/ALL_Merged_PCA.eigenvec", skip = 1, col.names = pc30k_header)
# 
# 
# table(cov$gxe)

# 1629 missing age-ref... 
# REMOVE FOR CONSISTENCY... 
# BUT KEEP IN MIND PCs were calculated on whole set N = 105689
# cov_pc30k <- inner_join(cov, pc30k, by = c('vcfid' = 'IID')) %>%
# 	mutate(age_ref = as.numeric(age_ref), 
# 				 sex = ifelse(sex == "Female", 0, 1),
# 				 outcome = ifelse(outc == "Case", 1, 0)) %>%
# 	dplyr::select(vcfid, outcome, age_ref, paste0('PC', seq(1,4)), sex)


# write file - GxEScanR
# cov_pc30k_noNA <- cov_pc30k[complete.cases(cov_pc30k), ]
# write.table(cov_pc30k_noNA, file = "~/FIGI_chr8_Covariates.txt", quote= F, row.names = F, col.names = T, sep = '\t')
# system("scp ~/FIGI_chr8_Covariates.txt hpcc:~ ")



# write file - EPACTS
# cov_pc30k_ped <- cov_pc30k %>%
# 	mutate(FAM_ID = vcfid, 
# 				 IND_ID = vcfid, 
# 				 FAT_ID = 0, 
# 				 MOT_ID = 0, 
# 				 SEX = sex, 
# 				 DISEASE = outcome) %>%
# 	dplyr::select(FAM_ID, IND_ID, FAT_ID, MOT_ID, SEX, DISEASE, age_ref, paste0('PC', seq(1,4)))
# cov_pc30k_ped_noNA <- cov_pc30k_ped[complete.cases(cov_pc30k_ped), ]
# 
# cat("#", file = "~/FIGI_chr8_Covariates_EPACTS.ped")
# write.table(cov_pc30k_ped_noNA, file = "~/FIGI_chr8_Covariates_EPACTS.ped", quote= F, row.names = F, col.names = T, sep = '\t', append = T)
# system("scp ~/FIGI_chr8_Covariates_EPACTS.ped hpcc:~ ")
# 












# 
# ##### ----- results quick process ------#####
# library(qqman)
# 
# # EPACTS
# #system("scp hpcc:/staging/dvc/andreeki/figi_chr8_epacts/figi_chr8_epacts_results.epacts.gz ~/")
# 
# epacts <- fread("~/figi_chr8_epacts_results.epacts", header = T) %>%
# 	mutate(CHR = `#CHROM`,
# 				 logp = -log10(PVALUE)) %>%
# 	filter(!is.na(PVALUE))
# qq(epacts$PVALUE)
# manhattan(epacts, chr = "CHR", bp = "END", p = "PVALUE", logp = T)
# 
# # sig result
# epacts_sig <- filter(epacts, logp > 5)
# 
# # GxEScan
# #system("scp hpcc:/staging/dvc/andreeki/test/merge/cytosnp_comb_merge_chr8.recode.results.out ~/")
# gxescan <- fread("~/cytosnp_comb_merge_chr8.recode.results.out", header = T) %>%
# 	mutate(pvalue = pnorm(zG, lower.tail = T),
# 				 logp = -log10(pnorm(zG, lower.tail = T)),
# 				 pvalue_gxe = pnorm(z_GE, lower.tail = T), 
# 				 logp_gxe = -log10(pnorm(z_GE, lower.tail = T)))
# 
# gxescan_gxe <- fread("~/cytosnp_comb_merge_chr8.recode.results.out", header = T) %>%
# 	mutate(pvalue = pnorm(zG, lower.tail = T),
# 				 logp = -log10(pnorm(zG, lower.tail = T)),
# 				 pvalue_gxe = pnorm(z_GE, lower.tail = T), 
# 				 logp_gxe = -log10(pnorm(z_GE, lower.tail = T))) %>%
# 	filter(!is.na(pvalue_gxe))
# 
# qq(gxescan$pvalue)
# manhattan(gxescan, chr = "CHR", bp = "BP", p = "pvalue", logp = T)
# 
# qq(gxescan_gxe$pvalue_gxe)
# manhattan(gxescan_gxe, chr = "CHR", bp = "BP", p = "pvalue_gxe", logp = T)
# 
# gxescan_sig <- filter(gxescan, logp > 5)
# 
# 
# # create df for talk
# 
# e <- epacts_sig %>%
# 	mutate(id = paste(CHR, END, sep = ":")) %>%
# 	dplyr::rename(epacts_logP = logp) %>%
# 	dplyr::select(id, epacts_logP, MAF)
# 
# g <- gxescan_sig %>%
# 	mutate(id = paste(CHR, BP, sep = ":")) %>%
# 	dplyr::rename(gxescan_logP = logp) %>%
# 	dplyr::select(id, gxescan_logP, Cases, Controls, A1, A2)
# 
# sigs <- full_join(e, g, by = 'id') %>% arrange(id) %>%
# 	dplyr::select(id, A1, A2, MAF, Cases, Controls, epacts_logP, gxescan_logP)
