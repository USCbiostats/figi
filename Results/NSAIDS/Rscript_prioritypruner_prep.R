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
#=============================================================================#
library(tidyverse)
library(data.table)
setwd("~/Dropbox/FIGI/Results/NSAIDS/")




# ------ Random sample of individuals ------
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

set.seed(2019)
cov_samp <- sample(cov$vcfid, 10000, replace = F)
cov_filter <- cov %>% 
  filter(vcfid %in% cov_samp)
write.table(cov_filter[, 'vcfid'], file = "./prune/PriorityPruner_10000_samplelist.txt", quote = F, row.names = F, col.names = F)

table(cov$study_gxe)
table(cov_filter$study_gxe)

# ------ GxE Results SNP List ------
# make sure you remove the 'bad' SNPs from UKB (~ 62000)
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

ukb_allele_check <- do.call(rbind, lapply(list.files(path = "~/data/UKB_AlleleFreq_Check/", full.names = T, pattern = "results_GWAS_ukb_allelefreq_check"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_"),
         P = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
  filter(!duplicated(ID))

# make sure you don't accidentally remove any of these (wanna keep them i guess... no matter what?)
x <- unique(fh_annotations$SNP)
tmp1 <- filter(ukb_allele_check, SNP %in% x,
               P < 5e-8)
tmp2 <- filter(ukb_allele_check, !SNP %in% x,
               P < 5e-8)
ukb_filter <- rbind(tmp1, tmp2)
any(duplicated(ukb_filter$SNP))
rm(tmp1); rm(tmp2); rm(x); rm(ukb_allele_check); rm(fh_annotations)


# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_72145"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID),
         !ID %in% ukb_filter$ID)

colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
ge <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref/ge_only/", full.names = T, pattern = "results_GEOnly_asp_ref_sex_age_pc10_studygxe_72145"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID),
         !ID %in% ukb_filter$ID) %>% 
  dplyr::select(-drop)


edgxe <- inner_join(gxe, ge[, c("ID", "chiSqGE")], by = "ID") %>% 
  mutate(pval.GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
         chiSqG_chiSqGE = chiSqG + chiSqGE,
         pval.edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F))

# pp <- filter(edgxe, Chromosome == 22)
# 
# write.table(pp[, c("Chromosome", "Location")], file = "./prune/PriorityPruner_chr22_snplist.txt", quote = F, row.names = F, col.names = F, sep = "\t")

for(chr in 1:21) {
  tmp <- filter(edgxe, Chromosome == chr)
  write.table(tmp[, c("Chromosome", "Location")], file = paste0("./prune/PriorityPruner_chr", chr, "_snplist.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
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






