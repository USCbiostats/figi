#=============================================================================#
# Clean up annotations
# use FIGI Controls as source for LD
# it'll be fine
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# Rsq estimate based on alt allele probability, maf > 0.001, rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")

# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/asp_ref/", full.names = T, pattern = "FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID))

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)


load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")

figi_gwasset_controls <- figi %>% 
  filter(outc == 'Control') %>% 
  dplyr::select(vcfid)

head(figi_gwasset_controls)
saveRDS(figi_gwasset_controls, file = "~/FIGI_controls_vcfid_73601.rds", version = 2)
#---------------------------------------------------#
# scrutinize the annotations ------
# peaks that are seemingly in LD with a GWAS top hit not being appropriately annotated
#---------------------------------------------------#


# create manhattan plots for every chromosome
for(chr in 1:22) {
  tmp <- filter(gxe, Chromosome == chr)
  create_manhattanplot(tmp, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = paste0("_BY_CHROMOSOME_chr", chr))
}


# chromosome 1
# my question is - for the third peak, by aren't markers physically near the GWAS hit not being
# labeled as high LD with that hit? 
chr1 <- filter(gxe, Chromosome == 1)

create_manhattanplot(chr1, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_TEST_chr1")

chr1_3 <- filter(chr1, Location >= 55261752-500000 & Location <= 55261752+500000)

create_manhattanplot(chr1_3, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_TEST_chr1_3")



# proof of concept..
# try to improve LD calculation by using the entire FIGI dataset (controls) 
# 1) find range of markers +- 500kb from GWAS hit
# 2) extract snps from binarydosage, subset on GWAS set CONTROLS
# 3) calculate LD, note r2 values, put them in a data frame that you can use for annotation EasyStrata

# 55261752 is the weird chromosome 1 hit with low LD but lots of elevated markers

chr1_bdinfo <- readRDS("~/data/BinaryDosage_InfoFile/FIGI_chr1.rds")[[11]] %>% 
  mutate(id = paste0(Chromosome, Location, Reference, Alternate))

chr1_bdinfo_hit <- chr1_bdinfo %>% 
  filter(Location >= 55261752-100000 & Location <= 55261752+100000)

chr1_bdinfo_index <- which(chr1_bdinfo$id %in% chr1_bdinfo_hit$id)
saveRDS(chr1_bdinfo_index, file = "~/huygue_gwas_140_chr1.rds", version = 2)
# saveRDS(chr1_bdinfo_index, file = "~/figi_asp_ref_chr1_55261752_index.rds", version = 2)

# After extracting dosages.... 
# note that gwasset without asians N ~ 131679
gwas_set <- readRDS("/home/rak/data/GxEScanR_PhenotypeFiles/FIGI_GWASSet_vcfid_N_138572.rds")

# this top hit should be the GWAS top hit for this locus
tophit <- readRDS("~/tmp/test4_out.rds") %>% 
  filter(vcfid %in% figi_gwasset_controls$vcfid) %>% 
  dplyr::select(X1.55246035)


# loop/apply whatever
file_list <- list.files("~/tmp", pattern = "test", full.names = T)

calculate_r2 <- function(x) {
  a <- readRDS(x) %>% 
    filter(vcfid %in% figi_gwasset_controls$vcfid) %>% 
    dplyr::select(-vcfid)
  
  b <- apply(a, 2, function(x) cor(x, tophit)^2)
  b
}

result <- do.call(c, lapply(file_list, calculate_r2) )

result_filter <- result[which(result > 0.2)]


# let's create a temporary data.frame to use as annotation for a single chromosome (1)
temp_annotation <- data.frame(r2 = result_filter, snp = names(result_filter)) %>% 
  separate(snp, into = c("Chr", "Pos")) %>% 
  mutate(Chr = gsub("X", "", Chr), 
         Colour = cut(r2, breaks = c(0.2, 0.5, 0.8, 0.999, Inf), labels = c("green", "orange", "red", "purple")))
write.table(temp_annotation, file = "~/data/Annotations/temp_annotation.txt", quote = F, row.names = F, sep = '\t')




# start with chromosome 1
plot_exposure <- "asp_ref"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")
chr1 <- filter(gxe, Chromosome == 1) %>% 
  filter(Location <= 55261752+500000 & Location >= 55261752-500000)

create_manhattanplot(chr1, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "TEST_chr1")


# chromosome 18
# because chromosome 6 is really weird
chr18 <- filter(gxe, Chromosome == 18)

chr6 <- filter(gxe, Chromosome == 6) %>% 
  filter(Location <= 31312538+1000000 & Location >= 31312538-1000000)


create_manhattanplot(chr18, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "TEST_chr18")




#----------------------------------------------------#
# old annotation file  ------
#----------------------------------------------------#
gwas95_old <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no",
         PREVIOUSLY_REPORTED == "yes") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "green", 
         SNP = paste(Chr, Pos, sep = ":"),
         ID = gsub("/", "_", VARIANT))

gwas95 <- readRDS("~/data/Annotations/Known_95_Loci_NatGenList_Jan2019.rds")

gwas95_df <- do.call(rbind, gwas95) %>% 
  dplyr::rename(Pos = end) %>% 
  dplyr::mutate(Chr = gsub("chr", "", chr),
                Colour = case_when(0 <= r2 & r2 < 0.2 ~ "lightblue",
                                   0.2 <= r2 & r2 < 0.4 ~ "lightblue", 
                                   0.4 <= r2 & r2 < 0.6 ~ "green", 
                                   0.6 <= r2 & r2 < 0.8 ~ "orange", 
                                   0.8 <= r2 & r2 <= 1.0 ~ "red"),
                marker = wtf2)

#rs16878812
gwas95_chr6 <- gwas95$`rs16878812_CEU-FIN-GBR-IBS-TSI.bed`

gwas_test <- gwas95[[35]]
