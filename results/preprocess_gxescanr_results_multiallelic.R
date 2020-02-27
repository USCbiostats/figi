#=============================================================================#
# GxEScanR Results pre-process
# 
# 1) filter Rsq > 0.8
# 2) create 3df stats, general clean up 
# 3) output statistics for LD CLumping
# 
# save object as .rds file. More convenient since you create plots repeatedly
#
# alternate script to fix issues with multi-allelic SNPs
# this is not an issue going forward, just for alcohol/HRT
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

exposure = 'gwas'
filename = 'FIGI_v2.3_gwasset_basic_covars_gxescan'

#-----------------------------------------------------------------------------#
# Read results, filter by Rsq
#-----------------------------------------------------------------------------#
# Rsq estimate based on alt allele probability, maf > 0.001, rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")

# Results
gxe_all <- do.call(rbind, lapply(list.files(path = paste0("~/data/results/", exposure), full.names = T, pattern = filename), fread, stringsAsFactors = F)) %>% 
  mutate(SNP = paste(SNP, Reference, Alternate, sep = ":"))

  #        chiSqEDGE = chiSqG + chiSqGE,
  #        chiSq3df = chiSqG + chiSqGxE + chiSqGE)

gxe_multi <- do.call(rbind, lapply(list.files(path = paste0("~/data/results/", exposure, "/multiallelic"), full.names = T, pattern = filename), fread, stringsAsFactor = F))

gxe_nomulti <- filter(gxe_all, !SNP %in% gxe_multi$SNP)

gxe <- rbind(gxe_nomulti, gxe_multi) %>% 
  filter(SNP %in% rsq_filter$id) %>% 
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  arrange(Chromosome, Location)


saveRDS(gxe, file = paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"), version = 2)


# # for gwas only - create smaller file for convenience
# gxe_short <- dplyr::select(gxe, SNP, betaG, chiSqG)
# saveRDS(gxe_short, file = paste0("~/data/results/", exposure, "/processed/", filename, "_results_short.rds"), version = 2)


#-----------------------------------------------------------------------------#
# Output chiSqG and chiSqGxE for LD Clumping
#-----------------------------------------------------------------------------#
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}

for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
    filter(Chromosome == chr) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/", filename, "_chr", chr, "_chiSqGxE_ldclump.txt"), quote = F, row.names = F, sep = '\t')
}

for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
    filter(Chromosome == chr) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/", filename,  "_chr", chr, "_chiSqG_ldclump.txt"), quote = F, row.names = F, sep = '\t')
}


#-----------------------------------------------------------------------------#
# Output chiSqG and chiSqGxE for LocusZoom Plots
# 
# how does locuszoom handle duplicate positions?? 
#-----------------------------------------------------------------------------#
locuszoom <- gxe %>%
  dplyr::mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
                MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
  dplyr::select(MarkerName, `P-value`)

write.table(locuszoom, file = paste0("/media/work/tmp/", filename, "_chiSqGxE_locuszoom.txt"), quote = F, row.names = F, sep = "\t")


locuszoom <- gxe %>%
  dplyr::mutate(`P-value` = pchisq(chiSqG, df = 1, lower.tail = F),
                MarkerName = paste0("chr", Chromosome, ":", Location)) %>% 
  dplyr::select(MarkerName, `P-value`)

write.table(locuszoom, file = paste0("/media/work/tmp/", filename, "_chiSqG_locuszoom.txt"), quote = F, row.names = F, sep = "\t")
