#========================================================#
# get rs numbers based on chr:bp information
# (in this case, to use locuszoom)
#
# using snp138 sql database on desktop
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(broom)
rm(list = ls())

# results
results <- do.call(rbind, lapply(list.files(path = "./results/", full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438_"), fread, stringsAsFactors = F))
names(results)

results_G <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaG, P)

# identify 'multiallelic' markers
dups <- results_G[duplicated(results_G$SNP), ]
dupped <- filter(results_G, SNP %in% dups$SNP) %>% arrange(desc(logP))

# Example locuszoom chromosome 8 peak
# generate sql script for results in chromosome 8
results_G_filter <- results_G %>% 
  mutate(logP = -log10(P)) %>% 
  filter(CHR == 8,
         !SNP %in% dupped$SNP) %>%  arrange(BP) %>% 
  mutate(sql = paste0('select chrom,chromEnd,name from snp138 where chrom = "chr8" and chromEnd = ', BP, ";"))
write.table(results_G_filter[,c('sql')], file = "get_rsnumbers.sql", quote = F, row.names = F, col.names = F)
#(run mysql script to get rsnumbers)
# mysql -u figi -pgecco -A -D hg19 --skip-column-names < get_rsnumbers.sql > results_rsnumbers_chr8.txt

#------------------------------------------------------------------------#

# final processing step - keep latest RS number if duplicated positions. 
rsnumbers <- fread("results_rsnumbers_chr8.txt", col.names = c("CHR", "BP", "RS")) %>% 
  arrange(desc(row_number())) %>% 
  filter(!duplicated(BP)) %>% 
  arrange(desc(row_number()))

z <- inner_join(results_G_filter, rsnumbers, by = "BP")

write.table(z[, c("SNP", "RS", "P")], file = "results_rsnumbers_chr8_locuszoom_Gmarginal.txt", quote = F, row.names = F, col.names = T, sep = '\t')

topP <- arrange(z, P)







#==========================================================================#
# GxE results
# plot all the loci in locuszoom.. 

results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_sig_GxE <- filter(results_GxE, P < 5e-8) # chromosomes 2, 3, 6, 7, 10, 16, 17, 18, 21, 22

results_GxE_filter <- results_GxE %>% 
  mutate(logP = -log10(P)) %>% 
  filter(!SNP %in% dupped$SNP) %>% arrange(BP)

# identify 'multiallelic' markers (or duplicates rows - maybe small gxescanR bug that captures rows twice)
dups <- results_GxE[duplicated(results_GxE$SNP), ]
dupped <- filter(results_GxE, SNP %in% dups$SNP)

# quick wrapper
sql_script <- function(chr) {
  tmp <- results_GxE %>% 
    mutate(logP = -log10(P)) %>% 
    filter(CHR == chr, 
           !SNP %in% dupped$SNP) %>% arrange(BP) %>% 
    mutate(sql = paste0('select chrom,chromEnd,name from snp138 where chrom = "chr', chr, '" and chromEnd = ', BP, ";"))
  write.table(tmp[,c('sql')], file = paste0("get_rsnumbers_chr", chr, ".sql"), quote = F, row.names = F, col.names = F)
}

# CHR2
sql_script(2)
sql_script(3)
sql_script(6)
sql_script(7)
sql_script(10)
sql_script(16)
sql_script(17)
sql_script(18)
sql_script(21)
sql_script(22)

#(run mysql script to get rsnumbers)
# mysql -u figi -pgecco -A -D hg19 --skip-column-names < get_rsnumbers_chr2.sql > results_rsnumbers_chr2.txt

#------------------------------------------------------------------------#
# final processing step - keep latest RS number if duplicated positions/rsnumbers

fix_rsdups <- function(chr) {
  rsnumbers <- fread(paste0("results_rsnumbers_chr", chr, ".txt"), col.names = c("CHR", "BP", "RS")) %>% 
    arrange(desc(row_number())) %>% 
    filter(!duplicated(BP)) %>% 
    arrange(desc(row_number()))
  
  pvals <- filter(results_GxE_filter, CHR == chr)
  
  z <- inner_join(pvals, rsnumbers, by = "BP")
  
  write.table(z[, c("SNP", "RS", "P")], file = paste0("results_rsnumbers_chr", chr, "_locuszoom.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
  
  # topP <- arrange(z, P)
  
}

for(x in c(2,3,6,7,10,16,17,18,21,22)) {
  fix_rsdups(x)
}



#==========================================================================#
# LD LINK
#==========================================================================#

# GxE results
# plot all the loci in locuszoom.. 

results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_sig_GxE <- filter(results_GxE, P < 5e-8) # chromosomes 2, 3, 6, 7, 10, 16, 17, 18, 21, 22

# identify 'multiallelic' markers (or duplicates rows - maybe small gxescanR bug that captures rows twice)
dups <- results_GxE[duplicated(results_GxE$SNP), ]
dupped <- filter(results_GxE, SNP %in% dups$SNP)


results_GxE_filter <- results_GxE %>% 
  mutate(logP = -log10(P)) %>% 
  filter(!SNP %in% dupped$SNP) %>% arrange(BP)

chr22 <- filter(results_GxE_filter, CHR == 22) # %>%
  # dplyr::select(CHR, BP, P)

write.table(chr22, file = "~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/results_rsnumbers_chr22_ldlink.txt", quote = F, row.names = F)
