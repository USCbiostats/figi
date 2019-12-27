#=============================================================================#
# GxEScanR Results pre-process
# 
# 1) filter Rsq > 0.8
# 2) create 3df stats, general clean up 
# 3) output statistics for LD CLumping
# 
# save object as .rds file. More convenient since you create plots repeatedly
#
# output a smaller file to merge with GxE scans for two-step methods
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# Arguments (run this script externally)
# keep in mind some paths are coded here, if you ever change files around
args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
filename <- args[2] # ex: FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820


#-----------------------------------------------------------------------------#
# Read results, filter by Rsq
#-----------------------------------------------------------------------------#
# Rsq estimate based on alt allele probability, maf > 0.001, rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")

# Results
gxe_all <- do.call(rbind, lapply(list.files(path = paste0("~/data/Results/", exposure), full.names = T, pattern = filename), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID))

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)

gxe_short <- dplyr::select(gxe, ID, betaG, chiSqG, betaGE, chiSqGE)

saveRDS(gxe_short, file = paste0("~/data/Results/", exposure, "/processed/", filename, "_results_short.rds"), version = 2)

