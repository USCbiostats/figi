#=============================================================================#
# PriorityPrunner
# outline:
# - subset vcf files
# - convert to plink
# - merge/concatenate
# - create p value table (EDGxE statistic)
# - use priority prunner (by chromosome might be more feasible)
#
# if you can edit SNP names in plink, might be useful to handle multi-allelics
#
# It would make things more manageable to take a random sample of 10,000 samples...
# i might do that for the in-person meeting
#
#
#
# Cleanest possible - use results that exclude UKB entirely
#=============================================================================#
library(tidyverse)
library(data.table)
setwd("~/Dropbox/FIGI/Results/NSAIDS/")


# ------ Random sample of individuals ------
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds") %>% 
  filter(study_gxe != "UKB_1")

set.seed(2019)
cov_samp <- sample(cov$vcfid, 10000, replace = F)
cov_filter <- cov %>% 
  filter(vcfid %in% cov_samp)
write.table(cov_filter[, 'vcfid'], file = "./prune/PriorityPruner_noUKB_10000_samplelist.txt", quote = F, row.names = F, col.names = F)

table(cov$study_gxe)
table(cov_filter$study_gxe)


# ------ GxE Results SNP List ------
# make sure you remove the 'bad' SNPs from UKB (~ 62000)
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_noUKB"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

table(gxe$Chromosome)

write.table(gxe[,'ID'], file = "./prune/PriorityPruner_noUKB_ALL_snplist.txt", quote = F, row.names = F, col.names = F)


for(chr in 1:22) {
  tmp <- filter(gxe, Chromosome == chr)
  write.table(tmp[, c("Chromosome", "Location")], file = paste0("./prune/PriorityPruner_noUKB_chr", chr, "_snplist.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}




#-----------------------------------------------------------------------------#
# ------ test run (successful) ------#
# managed to get priority prunner working with the list below..
# Begin with Chromosome 22, single batch (axiom_acs_aus_nf)

# subset vcf files
# output list of variants + list of vcfids 
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_72145"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

pp <- filter(gxe, Chromosome == 22)

# marker list 
write.table(pp[, c("Chromosome", "Location")], file = "~/data/PriorityPruner/results_GxE_asp_ref_chr22.txt", quote = F, row.names = F, col.names = F, sep = '\t')

# sample list
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

write.table(cov[, 'vcfid'], file = "~/data/PriorityPruner/results_GxE_asp_ref_chr22_samplelist.txt", quote = F, col.names = F, row.names = F)





#--------------------------------------------------------#
# pvalue list for PriorityPruner ----
#--------------------------------------------------------#
# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_noUKB"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

# original results
colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
ge <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/ge_only/", full.names = T, pattern = "results_GEOnly_asp_ref_sex_age_pc10_st"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID)) %>% dplyr::select(-drop)


# Required Columns:
#   name - Name of the SNP (e.g., rs222).
# chr - Name of the chromosome (e.g., 1, chr1). Chromosome X must be denoted as either 'X', 'chrX' or '23'. Chromosome Y and MT are not supported.
# pos - Physical position of the SNP.
# a1 - First allele of the SNP.
# a2 - Second allele of the SNP.
# p - P-value or other prioritization metric between 0 and 1. This is used for prioritizing the selection of SNPs, where lower numbers are prioritized.
# forceSelect - Flag indicating if the SNP should be selected (kept) regardless of its LD with other selected SNPs or other filtering criteria specified, such as MAF or design score (1=true, 0=false).
# designScore - Design score of the SNP (any positive real number). Can be filled in with a constant value (e.g., 1) if unknown.


# let's start with p values based on EDGxE statistic.. 

pppval <- inner_join(gxe, ge[, c("ID", "chiSqGE")], by = "ID") %>% 
  mutate(name = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate),
         chiSqG_chiSqGE = chiSqG + chiSqGE,
         p = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
         forceSelect = 0,
         designScore = 1) %>% 
  rename(chr = Chromosome, 
         pos = Location,
         a1 = Reference,
         a2 = Alternate) %>% 
  dplyr::select(name, chr, pos, a1, a2, p, forceSelect, designScore)

wtf <- pppval %>% 
  filter(chr == 22)

for(x in 1:22) {
  tmp <- pppval %>% 
    dplyr::filter(chr == x)
  write.table(tmp, file = paste0("~/priorityPruner_pvals_chr", x, ".txt"), quote = F, row.names = F)
}

head(pppval)

write.table(pppval, file = "~/test.txt", quote = F, row.names = F)




######### - on GxE statistic


pppval <- inner_join(gxe, ge[, c("ID", "chiSqGE")], by = "ID") %>% 
  mutate(name = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate),
         chiSqG_chiSqGE = chiSqG + chiSqGE,
         p = pchisq(chiSqGxE, df = 1, lower.tail = F),
         forceSelect = 0,
         designScore = 1) %>% 
  rename(chr = Chromosome, 
         pos = Location,
         a1 = Reference,
         a2 = Alternate) %>% 
  dplyr::select(name, chr, pos, a1, a2, p, forceSelect, designScore)

wtf <- pppval %>% 
  filter(chr == 22)

for(x in 1:22) {
  tmp <- pppval %>% 
    dplyr::filter(chr == x)
  write.table(tmp, file = paste0("~/priorityPruner_pvals_chr", x, ".txt"), quote = F, row.names = F)
}



# need to add sex information to tfam file
# remember in this file i coded female = 0
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

fam <- fread("~/test.tfam") %>% 
  inner_join(cov, by = c("V1" = "vcfid")) %>% 
  mutate(newsex = ifelse(sex == 0, 2, sex)) %>% 
  dplyr::select(V1, V2, V3, V4, newsex, V6)
head(fam)
table(fam$newsex)

write.table(fam, file = "~/test_sex.tfam", quote = F, row.names = F, col.names = F)




