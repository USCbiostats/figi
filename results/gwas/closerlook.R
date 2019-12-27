#=============================================================================#
# FIGI GWAS
#
# 1) Manhattan plots for each chromosome
# 2) followup with locuszoom plots (output a shell script here why not)
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

gxe <- readRDS("~/data/Results/gwas/processed/FIGI_GwasSet_age_sex_pc3_studygxe_124702_results.rds")

annotation <- fread("~/data/Annotations/temp_annotation_ver2.txt") %>% 
  dplyr::mutate(SNP = paste0(Chr, ":", Pos))
  
#-----------------------------------------------------------------------------#
# Manhattan plots by chromosome
#-----------------------------------------------------------------------------#
for(chr in 1:22) {
  tmp <- filter(gxe, Chromosome == chr)
  create_manhattanplot(tmp, 'gwas', c('age_ref_imp', 'sex', 'study_gxe', 'PC1', 'PC2', 'PC3'), stat = 'chiSqG', annotation_file = 'temp_annotation_ver2.txt', df = 1, filename_suffix = paste0("_Chr", chr))
}

#-----------------------------------------------------------------------------#
# Create locuszoom script using 'interesting' loci
# We only would care about  hits that persist after removing gwas hits + LD snps
#-----------------------------------------------------------------------------#
gxe <- gxe %>% 
  mutate(P = pchisq(chiSqG, df = 1, lower.tail = F))

wrapper <- function(chr, centers) {
  dplyr::filter(gxe, Chromosome == chr, 
                !SNP %in% annotation$SNP, 
                P < 5e-8) %>%
    dplyr::select(SNP, Chromosome, Location, Reference, Alternate, P) %>% 
    dplyr::mutate(group = kmeans(log(Location), centers = centers)[[1]]) %>% 
    arrange(group, P)
  }



chr1 <- wrapper(1, 4)
chr2 <- wrapper(2, 1)
chr4 <- wrapper(4, 1)
chr5 <- wrapper(5, 4)
chr6 <- wrapper(6, 4)
chr7 <- wrapper(7, 1)
chr8 <- wrapper(8, 4)
chr10 <- wrapper(10, 2)
chr11 <- wrapper(11, 4)
chr12 <- wrapper(12, 2)
chr14 <- wrapper(14, 4)

chr15 <- wrapper(15, 1)
chr17 <- wrapper(17, 1)

chr18 <- wrapper(18, 1)
chr19 <- wrapper(19, 1)
chr20 <- wrapper(20, 4)
chr21 <- wrapper(21, 1)


