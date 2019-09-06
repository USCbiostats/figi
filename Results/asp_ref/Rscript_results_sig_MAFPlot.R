#=======================================================================================#
# NSAIDS - asp_ref Results
# 05/18/2019
# 
# Extract G*NSAIDs GxE hits
# (output vector of index positions for binarydosage files)
#
# also perform some downstream analyses
#=======================================================================================#
library(tidyverse)
library(data.table)
library(purrr)
library(lmtest)
rm(list = ls())

# annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_sig
                        nals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_190518/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovT_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))


#--------------------------------------------------------#
# GxE results ----
# MAF information for 126 significant results
#--------------------------------------------------------#
# results_GxE <- gxe %>%
#   mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
#   dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_GxE_sig <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F),
         ID = paste0(SNP, ":", Reference, ":", Alternate)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P) %>% 
  filter(P < 5e-8)

info <- readRDS(file = "files/GxEScanR_GxE_asp_ref_Rsq_MAF.rds")

# for reference, if you want the Rsq values for the 126 markers for each imputation batch:
# info_filter <- info %>% 
#   dplyr::select(contains("Rsq")) %>% 
#   mutate_if(is.character, as.numeric)
# 
# info_final <- data.frame(mapply(`*`, info_filter, sample_sizes)) %>% 
#   mutate(Rsq_tot = rowSums(.[,1:16]),
#          Rsq_avg = Rsq_tot / sum(sample_sizes),
#          ID = info$ID) %>%
#   filter(Rsq_avg > 0.8)


# calculate MAFs, create plot to highlight differences between minimac3/minimac4
maf <- filter(info, ID %in% results_GxE_sig$ID ) %>%
  dplyr::select(contains("ALT_Frq")) %>%
  mutate(SNP = seq(1:nrow(.)))

mafplot <- gather(maf, key = 'study', value = 'aaf', -SNP) %>% 
  dplyr::mutate(study = c(rep("Others", 378), rep("ccfr_omni", 126), rep("Others", 1260), rep("Reach", 126), rep("UKB", 126)),
                MAF = 0.5 - abs(aaf - 0.5))

ggplot(mafplot, aes(SNP, MAF, colour = study)) + 
  geom_point(size = 0.5) + geom_line(size = 0.3) + ggtitle("MAFs for GxNSAID significant results") + xlab("SNPs 1-126")
ggsave("figures/GxEScanR_asp_ref_chiSqGxE_sig_MAFs.png", width = 7, height = 4)



# should calculate Rsq too 
test <- filter(info, ID %in% results_GxE_sig$ID ) %>%
  dplyr::select(contains("Rsq")) %>%
  mutate(SNP = seq(1:nrow(.)))

test2 <- gather(test, key = 'study', value = 'Rsq', -SNP) %>% 
  mutate(Rsqcat = ifelse(Rsq < 0.3, "Rsq < 0.3", ifelse(Rsq >=0.3 & Rsq < 0.8, "0.3 <= Rsq < 0.8", "Rsq >= 0.8")))
table(test2$Rsqcat)

tt <- test[1,] %>% 
  gather(key = "study", value = "Rsq", -SNP)

ggplot(test2, aes(SNP, Rsqcat, colour = study)) + 
  geom_point(size = 0.5)


ggsave("figures/GxEScanR_asp_ref_chiSqGxE_sig_MAFs.png", width = 7, height = 4)


