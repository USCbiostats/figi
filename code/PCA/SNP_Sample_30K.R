#======================================================================================
# 08/01/2018 (for full analysis)
#
# Get sample of SNPs for PCA, IBD 
# Ensure SNPs overlap with 1KGP reference panel 
# 
# * In imputation setting, considering using 'high quality' SNPs (Rsq > 0.99, MAF > 0.05)
# Since dealing with multiple batches, use corect_oncoarray (largest imputation set) 
# to subset markers
#
# Run this on HPC - SHOULD ONLY NEED TO BE RAN ONCE
#======================================================================================
library(dplyr)
library(data.table)
library(tidyr)
rm(list = ls())

#------ Directories ------#
setwd("/staging/dvc/andreeki/PCA") # working directory.. 
oncoarray_infofiles <- "/auto/pmd-02/figi/HRC_Unzipped/corect_oncoarray/"
kgp_plinkfiles <- "/auto/pmd-02/figi/andreeki/Refs/1KGP/"


#------ Read files, filter ------#
x <- do.call(rbind, lapply(paste0("zcat < ", list.files(path=oncoarray_infofiles, pattern = "info.gz", full.names = T)), fread, stringsAsFactors = F))

x.multiallelic <- x[duplicated(x$SNP),] # TO REMOVE MULTIALLELIC SNPS

x.edit <- filter(x, !SNP %in% x.multiallelic$SNP,
                 Genotyped != "Typed_Only", 
                 Rsq > 0.99 & MAF > 0.05)

table(x.edit$Genotyped)


#------ Use HRC tab file to convert SNP (CHR:BP) into rsID ------#
# (wrayner tools website)
hrc <- fread("grep -v ^# /auto/pmd-02/figi/andreeki/Refs/HRC.r1-1.GRCh37.wgs.mac5.sites.RSONLY.tab") %>%
  mutate(SNP = paste(V1, V2, sep = ":"))

x.hrc <- inner_join(x.edit, hrc, by = 'SNP') %>%
    dplyr::select(-V1, -V2) %>% dplyr::rename(ID = V3)

#------ Find overlapping SNPs with 1KGP ------#
# (see Folate_GWAS repository for converting 1KGP into plink files)
# (plink files stored in .../figi/andreeki/Refs/1KGP)
for(chr in 1:22) {
    tmp <- fread(paste0(kgp_plinkfiles, "kgp.chr", chr, ".biallelic.bim")) %>%
        dplyr::rename(ID = V2)
    overlap <- inner_join(x.hrc, tmp, by = 'ID')
    write.table(overlap$ID, file = paste0('overlap.chr', chr, '.rs.txt'), quote = F, row.names = F, col.names = F)
  }
kgp <- do.call(bind_rows, lapply(list.files(pattern = "overlap"), fread, stringsAsFactors = F, header = F))

# remove duplicates after joining hrc_kgp (kgp has duplicate IDs with different alleles)
dups <- kgp[duplicated(kgp$V1), ]

kgp.nodups <- filter(kgp, !V1 %in% dups$V1) %>%
  dplyr::rename(ID = V1)

# OVERLAP WITH RSIDs
hrc_kgp <- inner_join(x.hrc, kgp.nodups, by = "ID")
any(duplicated(hrc_kgp$ID))
write.table(hrc_kgp, file = "~/hrc_kgp_innerjoin_rsid.txt", quote = F, row.names = F) # takes a long time to get here, save as file


#----------------------------------------------------
# Take sample of SNPs
#----------------------------------------------------

# (this is after getting high quality SNPs that overlap with 1KGP)
set.seed(2018)
n.sample = 30000 # choosing 30k snps for final sample

hrc_kgp_sample <- hrc_kgp[sample(nrow(hrc_kgp), n.sample), ] %>%
	separate(SNP, into = c("Chr", "BP"), remove = F) %>%
	arrange(Chr, BP)

table(hrc_kgp_sample$Chr)

# output to text file
write.table(hrc_kgp_sample[,c("Chr", "BP"), ], file = "FIGI_PC_Backbone_Sample_30K.txt", quote = F, row.names = F, col.names = F)
write.table(hrc_kgp_sample[,c("ID"), ], file = "FIGI_PC_Backbone_Sample_30K_rsID.txt", quote = F, row.names = F, col.names = F)