#=============================================================================#
# 09/18/2019
#
# create LD based annotation for the all CRC top hits
# keep in mind you need this for each chromosome
#
# output - create list of IDs to calculate LD by chromosome
# format: 1 variant ID (CHR:BP:REF:ALT) per line
#=============================================================================#
library(tidyverse)
library(data.table)


# table reports SNPs with risk/non-risk alleles. Need REF/ALT
hits <-fread("~/Dropbox/FIGI/Documentation/Annotations_GWASHits_20190913/crc_gwas_indep_signals_050719_temp.tsv")[,-c(16,17)] %>% 
  mutate(SNP = paste(CHROM, POSITION, sep = ":"))

table(hits$CHROM)
any(duplicated(unique(hits$SNP)))

# not a perfect solution, but if you get multi-allelics it should still be ok, remove those from the output

for(chr in 1:20){
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))[[11]]
  tmp <- filter(figi, SNPID %in% hits$SNP) %>% 
    mutate(ID = paste(SNPID, Reference, Alternate, sep = ":"))
  write.table(tmp$ID, file = paste0("/media/work/tmp/plink_gwas_hit_ld_calculation_chr", chr, ".txt"), quote = F, row.names = F, col.names = F)
}


# Process the results, identify issues, redo if needed
ld <- do.call(rbind, lapply(list.files("~/data/Annotations/gwas_141/", full.names = T), fread))

tmp <- unique(paste(ld$CHR_A, ld$BP_A, sep = ":"))
miss <- hits[which(!hits$SNP %in% tmp), ]

# looks like there's a single SNP that is NOT found in HRC imputed markers - 6:105966894
# omit for now, it's still a hit in the process of being reported, probably from a different imputation panel, maybe based on fred hutch internal sequencing-based reference panel


#-----------------------------------------------------------------------------#
# re-create annotation file similar to the one provided by Stephanie
# i think the most important parts are the Chr, Pos, Colour
fh_annotations_old <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

summary(ld$R2)
summary(fh_annotations_old$r2)

ld <- do.call(rbind, lapply(list.files("~/data/Annotations/gwas_141/", full.names = T), fread)) %>% 
  rename(Chr = CHR_B,
         Pos = BP_B) %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"), 
         Colour = cut(R2, breaks = c(0.2, 0.4, 0.6, 0.8, 0.99, Inf), labels = c("lightblue", "green", "orange", "red", "purple")))
  
table(ld$Colour)

write.table(ld, file = "~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv", quote = F, row.names = F)

