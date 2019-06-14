#============================================
# FIGI August 2018
# 
# Create Summary kin0 File for Flora
# to make IBD drop decisions
# (sent via aspera)
#============================================
library(tidyverse)
library(data.table)

rm(list = ls())
setwd("~/Dropbox/code/FIGI_PCA_IBD_Results/")
# use this version for consistency
# newer samplefiles have different sample DROPS
# HOPEFULLY - there are no cases where they were initially DROPPED but are now KEPT
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_08162018.RData")

pid <- dplyr::select(All_Merged, vcfid, pooledcompassid, drop, gwas_set, studyname)

ibd <- read_table2("ALL_Merged_PCA_ibd.kin0") %>%
  dplyr::select(-FID1, -FID2)

ibd <- read_table2("ALL_Merged_PCA_drops_ibd.kin0") %>%
  dplyr::select(-FID1, -FID2)

all <- ibd %>%
  inner_join(pid, by = c("ID1" = "vcfid")) %>%
  dplyr::rename(pooledcompassid_ID1 = pooledcompassid, 
         drop_ID1 = drop,
         gwas_set_ID1 = gwas_set,
         studyname_ID1 = studyname) %>% 
  inner_join(pid, by = c("ID2" = "vcfid")) %>% 
  dplyr::rename(pooledcompassid_ID2 = pooledcompassid, 
         drop_ID2 = drop,
         gwas_set_ID2 = gwas_set,
         studyname_ID2 = studyname)


# separate dataframe into IDs that occur multiple times (ID1/ID2) vs only ONCE
x <- c(all$ID1, all$ID2)
all_table <- data.frame(table(x))

mult <- filter(all_table, Freq > 1)
all_multiple <- filter(all, ID1 %in% mult$x | ID2 %in% mult$x)
all_single <- filter(all, !ID1 %in% mult$x & !ID2 %in% mult$x)



#---------------------------------
# add information from samplefile

# add a 'family' variable
# why are there less samples after 'cleaning' though.. 
ids <- unique(c(all_multiple$ID1, all_multiple$ID2))
all_multiple_clean <- data.frame()
test <- ids
n <- 1
df <- all_multiple

for(id in test) {
  x1 <- filter(df, ID1 %in% id | ID2 %in% id) 
  xt <- unique(c(x1$ID1, x1$ID2))
  x2 <- filter(df, ID1 %in% xt | ID2 %in% xt)
  xt <- unique(c(x2$ID1, x2$ID2))
  x3 <- filter(df, ID1 %in% xt | ID2 %in% xt)
  xt <- unique(c(x3$ID1, x3$ID2))
  x4 <- filter(df, ID1 %in% xt | ID2 %in% xt)
  xt <- unique(c(x4$ID1, x4$ID2))
  x5 <- filter(df, ID1 %in% xt | ID2 %in% xt)
  xt <- unique(c(x5$ID1, x5$ID2))
  
  if(nrow(x5) == 0) {
    next
  } else {
    x5$family = n
    n <- n + 1
    df <- df[-which(df$ID1 %in% xt | df$ID2 %in% xt), ]
    all_multiple_clean <- bind_rows(all_multiple_clean, x5)
  }
  
}



# bind rows

out <- bind_rows(all_multiple_clean, all_single) %>% 
  dplyr::rename(vcfid_ID1 = ID1, 
                vcfid_ID2 = ID2) %>% 
  dplyr::select(InfType, vcfid_ID1, vcfid_ID2, pooledcompassid_ID1, pooledcompassid_ID2, gwas_set_ID1, gwas_set_ID2, studyname_ID1, studyname_ID2, N_SNP, HetHet, IBS0, HetConc, HomIBS0, Kinship, IBD1Seg, IBD2Seg, PropIBD, drop_ID1, drop_ID2, family) 


write.table(out, file = "~/FIGI_IBD_ALL_N145936.txt", quote = F, row.names = F, col.names = T, sep = '\t')
write.table(out, file = "~/FIGI_IBD_Exclude_Drops_N141362.txt", quote = F, row.names = F, col.names = T, sep = '\t')


# NEW FOR 11/02/2018
# write rds for decision processing
saveRDS(out, file = "IBD_Results_Aspera_20180816.rds")



# 
# table(All_Merged$drop)
# 
# 
# #=====================================================
# 
# 
# UKB <- filter(out, studyname_ID1 == "UKB_1" | studyname_ID2 == "UKB_1") %>% 
#   dplyr::select(-vcfid_ID1, -vcfid_ID2) %>% 
#   filter(studyname_ID2 != "UKB_1")
# 
# 
# 
# # adjust for gwas_set, and then studyname
# table(All_Merged$gwas_set)
# table(All_Merged$study)
# 
# # use this as indicator variables 
# table(All_Merged$studyname) # this sort of takes into account for study + genotyping platform. 

