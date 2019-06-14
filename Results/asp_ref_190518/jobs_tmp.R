library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
rm(list = ls())

# annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_190518/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovT_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

# calculate lambdas ------
getlambda <- function(pvals) {
  chisq <- qchisq(1-pvals, 1)
  lambda <- round(median(chisq)/qchisq(0.5,1),4)
  lambda
}
getlambda2df <- function(pvals) {
  chisq <- qchisq(1-pvals, 2)
  lambda <- round(median(chisq)/qchisq(0.5,2),4)
  lambda
}
getlambda3df <- function(pvals) {
  chisq <- qchisq(1-pvals, 3)
  lambda <- round(median(chisq)/qchisq(0.5,3),4)
  lambda
}
getlambda1000 <- function(lambda) {
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  lambda1000
}


#--------------------------------------------------------#
# GxE results ----
#--------------------------------------------------------#
results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

# Get SNP list for bdose.. 

x <- filter(results_GxE, P < 5e-8, Chromosome == '9')

