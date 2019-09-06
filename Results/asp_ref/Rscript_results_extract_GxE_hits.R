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
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_190518/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc3_studygxe_72820_binCovT_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))



#-----------------------------------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------------------------------#
results_G <- gxe %>%
  mutate(`P-value` = pchisq(chiSqG, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)

write.table(results_G, file = "~/Desktop/GxEScanR_aspref_age_sex_pc3_studygxe_N72820.txt", quote = F, row.names = F, sep = "\t")

tmp <- gxe %>%
  mutate(P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
  filter(Chromosome == 8, P < 5e-8) %>% 
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P)


#--------------------------------------------------------#
# GxE results ----
# use binarydosageinfo rdsfiles to get index positions
# for top hits (By Chromosome)
#--------------------------------------------------------#
results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P) %>% 
  filter(P < 5e-8)

table(results_GxE$Chromosome)
chroms <- unique(names(table(results_GxE$Chromosome)))


# loop, output index positions
for(chr in chroms) {
  figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))
  figi_snps <- figi[[11]] %>% 
    mutate(ID = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate))
  saveRDS(which(figi_snps$ID %in% results_GxE$ID), file = paste0("files/gxescanr_aspref_age_sex_pc3_studygxe_GxEresults_chr", chr, ".rds"), version = 2)
}


# output sample names too
saveRDS(cov$vcfid, file = "files/aspref_samples_N72820.rds", version = 2)


#--------------------------------------------------------#
# GxE results Rsq Values ------
#--------------------------------------------------------#
batch_list <- c("axiom_acs_aus_nf",
                "axiom_mecc_cfr_ky",
                "ccfr_1m_1mduo_reimpute",
                "ccfr_omni",
                "corect_oncoarray",
                "corsa_axiom",
                "cytosnp_comb",
                "initial_comb_datasets",
                "mecc",
                "newfoundland_omniquad",
                "omni_comb",
                "omniexpress_exomechip",
                "oncoarray_to_usc",
                "plco_3",
                "reach",
                "ukbiobank")
sample_sizes <- c(2766, 7501, 2180, 1600, 36621, 2467, 10908, 5924, 983, 637, 5986, 5503, 20912, 4864, 750, 27594)
chroms <- names(table(results_GxE$Chromosome))

# get info from each chromosome info file, create matrix
get_rsq <- function() {
  do.call(rbind, lapply(chroms, function(chr)
    readRDS(paste0("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr", chr, ".rds")) %>%
      filter(ID %in% results_GxE$ID)))
}

x <- get_rsq()
saveRDS(x, file = "files/GxEScanR_GxE_asp_ref_Rsq_MAF.rds", version = 2)

#--------------------------------------------------------#
# Garrett big lasso ------
# single hit on chr22, not related to NSAID scan
#--------------------------------------------------------#

# figi <- readRDS(paste0("~/data/BinaryDosage_InfoFile/FIGI_chr", 22, ".rds"))
# figi_snps <- figi[[11]] %>% 
#   mutate(ID = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate))
# 
# figi_snps_gar <- figi_snps %>% 
#   filter(Location == "48049700")
# 
# which(figi_snps$ID %in% figi_snps_gar$ID)
# saveRDS(which(figi_snps$ID %in% figi_snps_gar$ID), file = "~/garrett_tmp.rds", version = 2)
# 
# x<-readRDS("~/garrett_tmp.rds")


#--------------------------------------------------------#
# GLM ------
# (extracted files start with "getsnpvalues")
# (remember to perform lrtest)
#--------------------------------------------------------#
library(purrr)
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72820_GLM.rds")

# read in files (start with chr22)
# dose <- data.frame(readRDS("files/getsnpvalues_aspref_age_sex_pc3_studygxe_GxEresults_chr1.rds")) %>%
#   rownames_to_column("vcfid")
# z <- inner_join(cov, dose, by = 'vcfid'); names(z)
# m1 <- glm(outcome ~ `X22.16852708`*asp_ref + age_ref_imp + sex + study_gxe + PC1 + PC2 + PC3, data = z, family = 'binomial')
# summary(m1)


chr22 <- readRDS("files/getsnpvalues_aspref_age_sex_pc3_studygxe_GxEresults_chr22.rds")


tmp_function <- function(x) {
  # data.frame(readRDS(x)) %>% rownames_to_column('vcfid')
  data.frame(readRDS(x))
}

tmp <- list.files(path = "files", pattern = "getsnpvalues_aspref_age_sex_pc3_studygxe_GxEresults_chr", full.names = T) %>% 
  map_dfc(tmp_function) %>% 
  mutate(vcfid = rownames(chr22))

df <- inner_join(cov, tmp, by = 'vcfid')

  # correlation matrix
  names(df)
  cov_df <- data.frame(cor(df[, 100:105]))
  corrplot(cor(df[, 100:105]), method = 'color')

glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = df, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = df, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)


# run GLM + lrtest
basemodel <- map(df[ , 17:142], function(x) glm_func_base(x))
intmodel <- map(df[,17:142], function(x) glm_func(x))
results_gxe <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))

        # results match. Now to explore those loci further



#--------------------------------------------------------#
# locuszoom ------
# there are several regions, we should focus on one at 
# a time
#--------------------------------------------------------#

# Marginal G 
results_G_locuszoom <- results_G %>%
  mutate(`P-value` = P,
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)
write.table(results_G_locuszoom, file = "~/locuszoom/examples/GxEScanR_aspref_age_sex_pc3_studygxe_N72820_RsqFilter.txt", quote = F, row.names = F, sep = '\t')

tmp <- results_G %>% 
  filter(P < 5e-8,
         Chromosome == 1)

# GxE
results_GxE <- gxe %>%
  mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)

write.table(results_GxE, file = "~/Desktop/GxEScanR_GxE_aspref_age_sex_pc3_studygxe_N72820.txt", quote = F, row.names = F, sep = "\t")

tmp <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  filter(P < 5e-8,
         Chromosome == 8) %>% 
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)



 # 2DF
results_2df_locuszoom <- results_2df %>%
  mutate(`P-value` = P,
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)

write.table(results_2df_locuszoom, file = "~/locuszoom/examples/GxEScanR_2df_aspref_age_sex_pc3_studygxe_N72820_RsqFilter.txt", quote = F, row.names = F, sep = '\t')

tmp <- results_2df %>% 
  filter(P < 5e-8,
         Chromosome == 6)






# one time thing - create objects tomake getting rs numbers easier
# wtf <- fread("/media/work/FIGI/IC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
# 
# for(chr in 1:21) {
#   x <- filter(wtf, `#CHROM` == chr) %>% 
#     mutate(SNP = paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT))
#   saveRDS(x, paste0("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr", chr, ".rds"))
# }



library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_SNP", host = "grch37.ensembl.org", dataset = "hsapiens_snp")
rsid <- getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
      filters = c('chr_name','start','end'), 
      values = list(9, results_GxE$Location, results_GxE$Location), 
      mart = mart)


temp <- getBM(attributes = c('refsnp_id', 'allele', 'chrom_start', 'chrom_strand'), filters = c('chr_name', 'start', 'end'), values = list(chrom, df[snp, 2] , df[snp, 3]), mart = snp_mart)



# HRC reference panel info
hrc_chr9 <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr9.rds")

# chr9
# hits were between 9:65596318 -- 9:69780874

results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P) %>% 
  filter(P < 5e-8, 
         Chromosome == 9)



results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P) %>% 
  filter(Location > (65596318 - 500000) & Location < (69780874 + 500000), 
         Chromosome == 9)

z <- inner_join(hrc_chr9, results_GxE, by = c("SNP"="ID")) %>% 
  dplyr::rename(CHROM = `#CHROM`) %>% 
  mutate(logp = -log10(P),
         SNPnew = paste0("chr", CHROM, ":", Location))
write.table(z, file = "~/Dropbox/locuszoom_tmp.txt", quote = F, row.names = F, sep = '\t')



z <- inner_join(hrc_chr9, results_GxE, by = c("SNP"="ID")) %>% 
  rename(CHROM = `#CHROM`) %>% 
  mutate(logp = -log10(P)) %>% 
  dplyr::select(CHROM, Location, Reference, Alternate, P, ID)
write.table(z, file = "~/Dropbox/locuszoom_tmp.txt", quote = F, row.names = F, sep = '\t')



#-----------------------------------------------#
# GxE Results Table ------
# see rds object created above (with Rsq and Alt Freqs)
#-----------------------------------------------#
results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F),
         ID = paste(SNP, Reference, Alternate, sep = ":")) %>%
  dplyr::select(ID, Chromosome, Location, Cases, betaGxE, P)

info <- readRDS(file = "files/GxEScanR_GxE_asp_ref_Rsq_MAF.rds")

info_filter <- info %>% 
  dplyr::select(contains("Rsq")) %>% 
  mutate_if(is.character, as.numeric)

info_final <- data.frame(mapply(`*`, info_filter, sample_sizes)) %>% 
  mutate(Rsq_tot = rowSums(.[,1:16]),
         Rsq_avg = Rsq_tot / sum(sample_sizes),
         ID = info$ID) %>%
  filter(ID %in% results_GxE_sig$ID)



test <- info_filter %>% 
  mutate(ID = info$ID) %>% 
  filter(ID %in% results_GxE_sig$ID) %>% 
  separate(ID, into = c("chr", NA, NA, NA), sep = ":") %>% 
  group_by(chr) %>% 
  summarise_all(mean)

write.table(test, file = "~/Dropbox/GxEScanR_GxE_aspref_126SigHits_Rsq_By_ImpBatch.csv", sep = ",", quote = F, row.names = F)





# saveRDS(info_final[, c('ID', 'Rsq_tot', 'Rsq_avg')], file = "files/GxEScanR_GxE_asp_ref_Rsq_MAF_IDsOnly.rds", version = 2 )



## create the table of 'sig hits'
info_final <- data.frame(mapply(`*`, info_filter, sample_sizes)) %>% 
  mutate(Rsq_tot = rowSums(.[,1:16]),
         Rsq_avg = Rsq_tot / sum(sample_sizes),
         ID = info$ID) %>% 
  dplyr::select(ID, Rsq_avg)

results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F),
         ID = paste(SNP, Reference, Alternate, sep = ":")) %>%
  dplyr::select(ID, Chromosome, Location, Cases, betaGxE, P) %>% 
  filter(P < 5e-8)

z <- inner_join(results_GxE, info_final, by = 'ID')
saveRDS(z, file = "files/GxEScanR_GxE_asp_ref_sex_age_pc_studygxe_N72820_sigresults_table.rds", version = 2)
