# =============================================================================#
# use allele probability based Rsq estimates
# (function is written in binarydosage forked package in my account)
#
# create vector of SNP IDs to keep in analysis
# - maf > 0.001
# - Rsq >= 0.8
# =============================================================================#
library(tidyverse)
library(data.table)

#-----------------------------------------------------------------------------#
# re-estimated Rsq using latest merged FIGI bdose files 
# (dated August 2019)
#-----------------------------------------------------------------------------#

# wrapper to read-in results and filter by rsq > 0.8
wrap <- function(x) {
  readRDS(x) %>% 
    filter(rsq > 0.8)
}

results <- do.call(rbind, lapply(list.files("~/data/Rsq_Estimate/", pattern = "FIGI_RsqEstimate_chr", full.names = T), wrap))

results2 <- mutate(results, maf = 0.5 - abs(aaf - 0.5)) %>% 
  filter(maf > 0.001)

saveRDS(results2, file = "~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds", version = 2)

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # loop over chromosome, write temporary files
# wrapper_function <- function(chr) {
#   x <- readRDS(paste0("/home/rak/data/Rsq_Estimate/FIGI_RsqEstimate_chr", chr, ".rds"))
# 
#   y <- readRDS(paste0("~/data/HRC_InfoFile/FIGI/FIGI_chr", chr, ".rds"))$SNPs %>%
#     mutate(id = paste(SNPID, Reference, Alternate, sep = ":"))
# 
#   maf <- readRDS(paste0("~/data/HRC_InfoFile/FIGI/FIGI_chr", chr, ".rds"))$SNPInfo %>%
#     mutate(maf = 0.5 - abs(AAF - 0.5))
# 
#   z <- bind_cols(y, maf)
# 
#   if (length(x) == nrow(y)) {
#     z <- z %>%
#       mutate(Rsq_estimate = x) %>%
#       filter(
#         maf > 0.001,
#         Rsq_estimate >= 0.8
#       ) %>%
#       dplyr::select(id, Rsq_estimate)
#   } else {
#     z <- z %>%
#       mutate(Rsq_estimate = 0) %>%
#       filter(
#         maf > 0.001,
#         Rsq_estimate >= 0.8
#       ) %>%
#       dplyr::select(id, Rsq_estimate)
#   }
#   saveRDS(z, file = paste0("/home/rak/TEMPORARY_chr", chr, ".rds"))
# }
# 
# for(chr in 1:22) {
#   wrapper_function(chr)
# }
# 
# 
# # read in and concatenate
# zz <- do.call(rbind, lapply(list.files("~/", pattern = "TEMPORARY", full.names = T), readRDS))
# saveRDS(zz, file = "~/data/Rsq_Estimate/FIGI_RsqEsimate_chrALL.rds", version = 2)
# 
# 
# 
# 
# #-----------------------------------------------------------------------------#
# # worth creating a non-filtered file for use by others
# #-----------------------------------------------------------------------------#
# 
# wrapper_function_two <- function(chr) {
#   x <- readRDS(paste0("/home/rak/data/Rsq_Estimate/FIGI_RsqEstimate_chr", chr, ".rds"))
#   
#   y <- readRDS(paste0("~/data/HRC_InfoFile/FIGI/FIGI_chr", chr, ".rds"))$SNPs %>%
#     mutate(id = paste(SNPID, Reference, Alternate, sep = ":"))
#   
#   maf <- readRDS(paste0("~/data/HRC_InfoFile/FIGI/FIGI_chr", chr, ".rds"))$SNPInfo %>%
#     mutate(maf = 0.5 - abs(AAF - 0.5))
#   
#   z <- bind_cols(y, maf)
#   
#   if (length(x) == nrow(y)) {
#     z <- z %>%
#       mutate(Rsq_estimate = x) %>%
#       dplyr::select(id, maf, Rsq_estimate)
#   } else {
#     z <- z %>%
#       mutate(Rsq_estimate = 0) %>%
#       dplyr::select(id, maf, Rsq_estimate)
#   }
#   saveRDS(z, file = paste0("/home/rak/FIGI_GWASSet_138015_MAF_RsqEstimate_chr", chr, ".rds"))
# }
# 
# for(chr in 1:22) {
#   wrapper_function_two(chr)
# }



