#=============================================================================#
# NSAIDS - asp_ref
# replicate Nan 2015 analysis and findings
# 06/25/2019
# 
# Handful of significant findings, are they also significant in our data? 
# recreate some of the meta analyses
# run GxE
# consider running a full scan?
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(animation)
# rm(list = ls())

global_N <- unique(gxe$Subjects)
global_covs <- c("age_ref_imp", "sex", "study_gxe", "PC1-3")
global_E <- "asp_ref"



# Start with asp_ref variable (main findings were for combined aspirin and/or NSAIDs)
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")
pc30k <- fread("~/data/PCA/190506/FIGI_GxESet_KGP_pc20_190430.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))


    # Investigate PMH-CCFR
    x <- data.frame(table(figi_samplefile$study, figi_samplefile$studyname)) %>% 
      filter(Freq != 0)
    
    y <- data.frame(table(figi$study, figi$studyname)) %>% 
      filter(Freq != 0)
    
    z <- filter(figi_samplefile, grepl("PMH", study_gxe))
    # many of them are drop == 1 (reason: "Related Genotyping Platform". I think it was from IBD analysis)
    # let's just not include them in this exploratory analysis
    
    # Investigate CCFR, OFCCR
    # analysis keeps Australia, Toronto (Ontario), Seattle, and OFCCR (Ontario)
    z <- filter(figi_samplefile, grepl("OFCCR", studyname)) # most of these are dropped since large overlap with ccfr ontario
    
    z <- filter(figi_samplefile, grepl("CCFR", studyname))
    table(z$studyname, z$study_site) # study_site captures the different CCFR studies


# ------ Recreate Nan Subset ------ #
# Exclude PMH-CCFR for this run... small N anyway

study_keep <- c("CCFR_1", "CCFR_2", "CCFR_3", "CCFR_4", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PLCO_1_Rematch", "PLCO_2", "PLCO_3", "PLCO_4_AD", "VITAL", "WHI_1", "WHI_2", "WHI_3")

# first, explore counts for studies with aspirin information available
pheno <- figi %>% 
  filter(drop == 0 & gxe == 1) %>%
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>%
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) %>% 
  filter(!is.na(asp_ref),
         study_gxe %in% study_keep,
         !study_site %in% c("Hawaii", "Los Angeles", "Mayo Foundation"))

table(pheno$study, pheno$outcome)  
x <- data.frame(table(pheno$study, pheno$year_ref)) %>% filter(!Freq == 0)


# DACHS - remove "OmniExpress_ExomeChip"
dachs <- pheno %>% 
  filter(grepl("DACHS", study_gxe),
         race_self == "White", 
         platform != "OmniExpress+ExomeChip")
table(dachs$year_ref)


# DALS - use all? 
dals <- pheno %>% 
  filter(grepl("DALS", study_gxe),
         race_self == "White")


# HPFS - omniexpress only?
hpfs <- pheno %>% 
  filter(grepl("HPFS", study_gxe),
         race_self == "White",
         platform == "OmniExpress")
table(hpfs$outcome)


# NHS - omniexpress only?
nhs <- pheno %>% 
  filter(grepl("NHS", study_gxe),
         race_self == "White",
         platform == "OmniExpress")
table(nhs$outcome)


# PLCO - cytoSNP
plco <- pheno %>% 
  filter(grepl("PLCO", study_gxe), 
         race_self == "White",
         platform == "CytoSNP")
table(plco$outcome)


# VITAL - use all?
vital <- pheno %>% 
  filter(grepl("VITAL", study_gxe), 
         race_self == "White")


# WHI - 
whi <- pheno %>% 
  filter(grepl("WHI", study_gxe), 
         race_self == "White",
         platform != "OncoArray+Custom")
table(whi$platform)

# CCFR 
ccfr <- pheno %>% 
  filter(grepl("CCFR", study_gxe),
         race_self == "White",
         platform == "Illumina_1M/1Mduo")

# ------ merge all ------ #

pheno_cc <- rbind(ccfr, dachs, dals, hpfs, nhs, plco, vital, whi) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, study, study_gxe, study_site, paste0(rep("PC", 3), seq(1,3)), asp_ref)





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
