# Compare sample files 181105 vs 181120
library(tidyverse)
library(data.table)

rm(list = ls())

exclude <- c("BWHS", "COLON", "DACHS_4", "NSHDS", "OSUMC", "TCGA", "OFCCR")

older <- fread("~/git/DATA/FIGI_samplefile_epi-181105/Samplefile_HRC1.1_epi_v2-20181105.txt") %>% 
  filter(!studyname %in% exclude)

newer <- fread("~/git/DATA/FIGI_samplefile_epi-181120/Samplefile_HRC1.1_epi_v2-20181120_usc.txt") %>% 
  filter(!studyname %in% exclude)

x <- anti_join(older, newer, by = "vcfid")

table(x$studyname)




# G set

older.clean <- filter(older, !vcfid %in% x$vcfid)


a <- dplyr::select(older.clean, vcfid, studyname, drop, gxe, drop_reason) 
b <- dplyr::select(newer, vcfid, studyname, drop, gxe, drop_reason) 

c <- inner_join(a, b, by='vcfid') %>% 
  #filter(drop.x != drop.y) %>% 
  filter(gxe.x != gxe.y) %>% 
  #dplyr::select(vcfid, studyname.x, studyname.y, drop.x, drop.y, drop_reason.y) %>% 
  dplyr::select(vcfid, studyname.x, studyname.y, gxe.x, gxe.y, drop_reason.y) %>% 
  filter(studyname.x %in% c("NFCCR_1", "NFCCR_2"))
