#=============================================================================#
# NSAIDS - asp_ref
# replicate Nan 2015 analysis and findings
# 11/25/2019
# 
# Using compassids provided by Yi
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(animation)
rm(list = ls())

global_N <- unique(gxe$Subjects)
global_covs <- c("age_ref_imp", "sex", "study_gxe", "PC1-3")
global_E <- "asp_ref"


# Start with asp_ref variable (main findings were for combined aspirin and/or NSAIDs)
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")
pc30k <- fread("~/data/PCA/190729/FIGI_GxESet_190729.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

# Yi ID List for Nan analysis
idlist <- fread("~/data/Nan_etal/GxNSAIDs_IDs.txt")
any(duplicated(idlist$compassid))

# create studyname_compassid var
idlist <- fread("~/data/Nan_etal/GxNSAIDs_IDs.txt") %>% 
  mutate(tmpid = paste0(studyname, "_", compassid))

dat <- figi %>% 
  mutate(tmpid = paste0(studyname, "_", compassid)) %>% 
  filter(tmpid %in% idlist$tmpid)

non_matches <- anti_join(idlist, dat, 'tmpid')

wtf <- figi_samplefile %>%
  mutate(tmpid = paste0(studyname, "_", compassid)) %>%
  filter(tmpid %in% non_matches$tmpid)

table(wtf$drop) # all dropped from FIGI.. 


#-----------------------------------------------------------------------------#
# meta-analysis of asp_ref in this subset ----
#-----------------------------------------------------------------------------#
tmp <- pheno_cc %>%
  mutate(study_gxe = study)

x <- getCounts_byOutcome(tmp, outcome, asp_ref)
y <- getGLM_byGroup(tmp, outcome, asp_ref, c("age_ref_imp", "sex"))
run_meta_analysis_create_forestplot(x, y, "aspirin use @ ref time \nModel: outcome~aspirin+age_ref_imp+sex")


#-----------------------------------------------------------------------------#
# Significant SNPs ----
# rs2965667 chr12:17444733
# rs10505806 chr12:17488764
# rs16973225 (case-only) chr15:82229999
#-----------------------------------------------------------------------------#

# I'm using a function I wrote previously..

# sample subset for GetSNPValues
saveRDS(pheno_cc$vcfid, file = paste0(write_binarydosage_vcfid_filename(), ".rds"), version = 2)

# snps (vector of index positions)
rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds")

x <- c()
x <- c(x, filter(rsq_filter, grepl("12:17444733", ID))[1,1])
x <- c(x, filter(rsq_filter, grepl("12:17488764", ID))[1,1])

y <- filter(rsq_filter, grepl("15:82229999", ID))[1,1]

get_binarydosage_index(x, 12)
get_binarydosage_index(y, 15)



#-----------------------------------------------------------------------------#
# Merge Dosages, run GLM ----
#-----------------------------------------------------------------------------#
tmp <- data.frame(readRDS("files/GetSNPValues_asp_ref_age_ref_imp_sex_study_gxe_PC1-3_N_72820_chr12_out.rds")) %>% 
  rownames_to_column('vcfid')

pheno_cc_dosages <- inner_join(pheno_cc, tmp, by = 'vcfid')

summary(glm(outcome ~ X12.17444733*asp_ref + age_ref_imp + sex + study_gxe + study_site + PC1 + PC2 + PC3, data = pheno_cc_dosages, family = 'binomial'))

summary(glm(outcome ~ X12.17488764*asp_ref + age_ref_imp + sex + study_gxe + study_site +  PC1 + PC2 + PC3, data = pheno_cc_dosages, family = 'binomial'))
