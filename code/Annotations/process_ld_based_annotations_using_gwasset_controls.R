#=============================================================================#
# LD Based Annotations (Huygue 2019) 
# Does a few things related to GWAS annotations and expectation based two-step
# 
# 1) calculate r2 for all markers flanking reported hits +- 1MB
#    - 1MB is necessary because some r2 > 0.2 are actually pretty far away from
#      main hit. 
#    - using FIGI controls (N ~ 73,601)
# 2) get lists of SNPs related to GWAS top hits to remove from the analysis
#    and treat separately as markers already discovered
#    you can do LD vs regional markers. LD misses some SNPs, regional is
#    arbitrary. Maybe I can try anything within an LD region (flanked by the 
#    outer SNPs in LD with main hit with r2 > 0.2
#    
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# FIGI Controls vcfid
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
figi_gwasset_controls <- figi %>% 
  filter(outc == 'Control') %>% 
  dplyr::select(vcfid)

# Jeroen GWAS hits
# important - EXCUDE 6:105966894 (does not appear to be in HRC reference panel)
gwas <-fread("~/Dropbox/FIGI/Documentation/Annotations_GWASHits_20190913/crc_gwas_indep_signals_050719_temp.tsv")[,-c(16,17)] %>% mutate(SNP = paste(CHROM, POSITION, sep = ":")) %>% 
  dplyr::filter(SNP != "6:105966894") %>% 
  arrange(SNP)

# output BinaryDosage index position for 140 hits above, +- 1MB flanking region
# extract dosage values from BinaryDisage files, calculate r2 among controls
  # just a bit of useful regex
  # gsub( ":.*$", "", gwas$SNP)
  # gsub(".*:","",gwas$SNP)

# important - note that I picked a region of 2MB. Sometimes markers in correlation with the top hit are actually very distant in terms of base pair location. That's what happened to chromosome 12 initially

snps <- gwas$SNP
for(snp in snps){
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", gsub( ":.*$", "", snp), ".rds"))[[11]] %>% 
    mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
  
  hit_location <- as.numeric(gsub(".*:","",snp))
  
  tmp <- figi %>% 
    filter(Location >= hit_location-1000000 & Location <= hit_location+1000000)
  
  tmp_index <- which(figi$id %in% tmp$id)
  
  saveRDS(tmp_index, file = paste0("~/huygue_gwas_140_", snp, ".rds"), version = 2)
}

# NOTE - make sure you rename the .rds files to remove colons because it causes problems with file transfers etc 
# rename  's/:/_/' huygue_gwas_140_* 




#--------------------------------------------------#
# get dosage for the 140 top hits only
#--------------------------------------------------#
snps <- gwas$SNP
# snps <- "17:10707241"
for(snp in snps){
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", gsub( ":.*$", "", snp), ".rds"))[[11]] %>% 
    mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
  
  tmp <- figi %>% 
    filter(SNPID == snp)
  
  tmp_index <- which(figi$id %in% tmp$id)
  
  saveRDS(tmp_index, file = paste0("~/huygue_gwas_140_tophitonly_", snp, ".rds"), version = 2)
}


# rename  's/:/_/' huygue_gwas_140_tophitonly_* 





#--------------------------------------------------#
# GetSNPValues Results
# 1) clean up output
# 2) create annotation file for EasyStrata
#--------------------------------------------------#
file_list <- list.files("~/data/Annotations/gwas_140_r2_figi_controls/", pattern = "r2", full.names = T)
ld <- do.call(c, lapply(file_list, function(x) readRDS(x)))
ld_filter <- ld[which(ld > 0.2)]
names(ld_filter) <- gsub("\\.1$", "", names(ld_filter)) # some names have a ".1" at the end.. remove and deal with duplicates

# duplicates - keep the ones with the highest r2
x <- ld_filter[order(names(ld_filter), -ld_filter)] # keep duplicate with the highest r2 value..
# (although.. ultimately it doesn't matter since you're filtering only on chr:bp)
ld_filter_nodups <- x[!duplicated(names(x))]

temp_annotation <- data.frame(r2 = ld_filter_nodups, snp = names(ld_filter_nodups)) %>% 
  dplyr::filter(!grepl('\\.1$', snp)) %>%
  separate(snp, into = c("Chr", "Pos")) %>% 
  mutate(Chr = gsub("X", "", Chr), 
         Colour = cut(r2, breaks = c(0.2, 0.5, 0.8, 0.999, Inf), labels = c("green", "orange", "red", "purple")))
write.table(temp_annotation, file = "~/data/Annotations/temp_annotation_ver2.txt", quote = F, row.names = F, sep = '\t')



# some tests...






# ---------------------------------------------------- #
# try to create a manhattan plot with these annotations
# using asp_ref as example
gxe <- readRDS("~/data/Results/asp_ref/processed/FIGI_GxESet_asp_ref_age_sex_pc3_studygxe_72820_results.rds")
plot_exposure <- "asp_ref"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', annotation_file = "temp_annotation_ver2.txt", df = 1)

# manhattan plot after removing these SNPs
# to check for left over significant hits, and to make sure these are independent, not part of a previously reported GWAS hit
temp_annotation <- temp_annotation %>% 
  mutate(SNP = paste0(Chr, ":", Pos))
gxe_remove_annotations <- filter(gxe, !SNP %in% temp_annotation$SNP)
create_manhattanplot(gxe_remove_annotations, plot_exposure, plot_covariates, stat = 'chiSqG', annotation_file = 'temp_annotation_ver2.txt', df = 1, filename_suffix = "_NO_GWAS_HITS")

# some left over peaks after removing "top hits" from the results file. 
# these need further investigation - create locuszoom plots for each one of those
# (Although i'd wait until main effects GWAS is finished to scrutinize them)
# the point is, you have top hits to remove now for separate analysis

# but first, implement LD + region based filters (could be within region but not in LD with top hit)


#--------------------------------------------------#
# using LD information, get any SNP that falls inside LD region.. 
# (because sometimes SNPs that are inside an LD region are actually not 
# in any correlation with the top hit)
#--------------------------------------------------#
gxe <- readRDS("~/data/Results/gwas/processed/FIGI_GwasSet_age_sex_pc3_studygxe_124702_results.rds")
file_list <- list.files("~/data/Annotations/gwas_140_r2_figi_controls/", pattern = "r2", full.names = T)

# read each rds file, then filter r2 > 0.2
# output a list
read_filter_wrapper <- function(x) {
  tmp <- readRDS(x)
  tmp_filter <- tmp[which(tmp >= 0.2)]
  names(tmp_filter) <- gsub("\\.1$", "", names(tmp_filter)) # remove ".1", deal with duplicates
  out <- tmp_filter[order(names(tmp_filter), -tmp_filter)]
  out_filter_nodups <- out[!duplicated(names(out))]
  # keep in mind some chromosomes have overlapping regions
}

ld_filter <- lapply(file_list, read_filter_wrapper)


# takes each list item, get lower and upper basepair bound
# output a list too
range_wrapper <- function(x) {
  names_tmp <- as.numeric(gsub(".*\\.","",names(x)))
  
  tmp <- names(x)[order(names_tmp)][c(1, length(x))]
  tmp_noX <- gsub("X", "", tmp)
  tmp_noC <- gsub("\\.", ":", tmp_noX)
}
ld_filter_window <- lapply(ld_filter, range_wrapper)

# get a list of SNPs that you can exclude from results
# these SNPs removed because they fall within a region
# that contains SNPs in LD with the gwas top hit
gxe_tmp <- gxe %>% 
  dplyr::select(SNP, Chromosome, Location)

range_get_wrapper <- function(x) {
  chr <-  as.numeric(gsub( ":.*$", "", x))
  position <- as.numeric(gsub(".*:","",x))
  
  tmp <- gxe_tmp %>% 
    filter(Chromosome == chr[1], 
           between(Location, position[1], position[2]))
}

test <- lapply(ld_filter_window, range_get_wrapper)
out <- do.call(rbind, test)

saveRDS(out, file = "~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds", version = 2)


# some SNPs in the annotation file aren't in the results
# maybe it's due to MAF, not being in HRC? Rsq maybe. 
# wtf <- anti_join(temp_annotation, out, by = 'SNP')


# create manhattan plot to check
gxe_filter <- gxe %>% 
  filter(!SNP %in% out$SNP)

create_manhattanplot(gxe_filter, plot_exposure, plot_covariates, stat = 'chiSqG', annotation_file = 'temp_annotation_ver2.txt', df = 1, filename_suffix = "_NO_GWAS_HITS_REGIONS")



