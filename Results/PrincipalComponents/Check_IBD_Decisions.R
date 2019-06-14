#============================================
# FIGI 11/02/2018
# 
# After UW made drop decisions, some still 
# remaining in dataset?
#
# Investigate, email Flora
#============================================
library(tidyverse)
library(data.table)
library(reshape)

rm(list = ls())
setwd("~/Dropbox/code/FIGI_PCA_IBD_Results/")
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_08162018.RData")

# compare drops with original results to see if some decisions were missed.
# easier to see by melting dataset
ibd_original <- readRDS("IBD_Results_Aspera_20180816.rds") %>% 
  dplyr::select(vcfid_ID1, vcfid_ID2, family) %>% 
  mutate(group = c(.$family[1:684], seq(212,3082))) %>% # just a simple group var for dups/relateds
  dplyr::select(-family) %>% 
  melt(id = 'group') %>% 
  dplyr::select(-variable) %>%filter(!duplicated(.)) %>% arrange(group) %>% 
  dplyr::rename(vcfid = value)

ibd_flora <- fread("~/git/DATA/FIGI_samplefile_epi-103118/FIGI_IBD_pairs_drop_decision_single_entry.txt") %>% 
  mutate(IBD_drop = "yes")

ibd_original_decisions <- full_join(ibd_original, ibd_flora, by = c("vcfid" = "drop_sample"))

write.csv(ibd_original_decisions, file = "IBD_Original_Summary.csv", quote = F, row.names = F)

# check newest analysis
ibd_new <- fread("FIGI_GxESet.kin0") %>% 
  dplyr::select(ID1, ID2) %>% 
  mutate(issue = 'yes')


ibd_final <- full_join(ibd_original_decisions, ibd_new, by = c("value" = "ID1"))



check_drop_stat <- inner_join(ibd_new, Epi, by = c('ID2' = 'vcfid'))
table(check_drop_stat$drop)
  