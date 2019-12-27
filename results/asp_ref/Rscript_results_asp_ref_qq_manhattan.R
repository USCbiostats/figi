#=============================================================================#
# NSAIDS - asp_ref Results
# 05/18/2019
# 
# Generate Plots
# Implement 2-step methods
#
# let's see what things look like after filtering by Rsq
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


# you want to make sure NOT to drop annotations from manhattan plots. let's make sure there's good overlap...
# 90 OF THE MARKERS ARE IN MY RESULTS, MISSING 5, 2 of them were previously reported
# annotation_original <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
#   filter(SUGGESTIVE=="no") %>% 
#   separate(VARIANT, into = c("SNP", "alleles"), sep = "_")
# 
# check <- inner_join(gxe, annotation_original, by = "SNP")
# check <- anti_join(annotation_original, gxe, by = "SNP")



# output chiSqGxE results for LD clumping
# calculate_pval <- function(data, statistic, df) {
#   data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
#   data
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_asp_ref_ldclump_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_asp_ref_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }

# LD Clump Results
tmp <- do.call(rbind, lapply(list.files("~/data/Results/asp_ref/clump_controls/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)






#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "asp_ref"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")

# Marginal G Results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)

# GxE results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)

# 2DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)

# 3DF results ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)

# GE, Case, Control ----
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1)
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)


# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')


#-----------------------------------------------------------------------------#
# QQ and Manhattan Clumped ----
#-----------------------------------------------------------------------------#

# manhattan plot
gxe_plot <- gxe_chiSqGxE_ld %>% 
  mutate(p = pchisq(chiSqGxE, df=1, lower.tail=F)) %>% 
  filter(-log10(p) >=6 )
names(gxe_plot)

manhattan(gxe_plot, chr = "Chromosome", bp = "Location", p = "p")



#-----------------------------------------------------------------------------#
# Two-Step Clumped ----
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = '_LDclump')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_LDclump")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_LDclump")










#-----------------------------------------------------------------------------#
# Extract dosage values in significant loci from binary dosage files ----
# (needs hand holding - look at manhattan plots to identify significant hits)
# (for GetSNPValues - need sample list and SNP index positions)
#-----------------------------------------------------------------------------#

# Significant Results
sig1 <- filter(gxe, chiSqGxE > 30) #GxE
sig2 <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p <= wt) #two step method, G and EDGE same answer

sig_all <- unique(c(sig1$ID, sig2$ID))

# get SNP index positions (one file per chromosome)
chroms <- unique(sig2$Chromosome)
for(chr in chroms) {
  figi <- readRDS(paste0("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr", chr, ".rds"))
  figi_snps <- figi[[11]] %>% 
    mutate(ID = paste0(Chromosome, ":", Location, ":", Reference, ":", Alternate))
  saveRDS(which(figi_snps$ID %in% sig_all), file = paste0("files/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr", chr, ".rds"), version = 2)
}

# write sample list (vector)
vcfid <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820_GLM.rds")[, 'vcfid']
saveRDS(vcfid, file = "files/GetSNPValues_asp_ref_sex_age_pc3_studygxe_vcfid.rds", version = 2)



# ... 
# ... 
# ...


# confirm results with GLM 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/asp_ref/dosage/", pattern = "GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820_GLM.rds"))

glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(dat[ , 2:3], function(x) glm_func_base(x))
intmodel <- map(dat[,2:3], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))



# for X6.32560631, might be driven by UKB1.. 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/asp_ref/dosage/", pattern = "GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820_GLM.rds"))

dat <- filter(dat, study_gxe != "UKB_1")

glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

basemodel <- map(dat[ , 2:3], function(x) glm_func_base(x))
intmodel <- map(dat[,2:3], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))





#-----------------------------------------------------------------------------#
# Curious - are flanking SNPs in chromosome 6 also different MAFs for UKB? 
#-----------------------------------------------------------------------------#

# get chromosome 6 +- 10 SNPs in flanking region
test <- readRDS("files/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr6.rds")
testout <- seq(test-10, test+10)
saveRDS(testout, file = paste0("files/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_TEST_chr", 6, ".rds"), version = 2)


# calculate mafs by study for each of these markers
dat_test <- data.frame(readRDS("~/data/Results/asp_ref/dosage/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_TEST_chr6_out.rds")) %>% 
  rownames_to_column("vcfid") %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc3_studygxe_72820_GLM.rds"))

dat_test_maf <- dat_test %>% 
  group_by(study_gxe) %>% 
  summarise_at(vars(X6.32560327:X6.32560839), function(x) 0.5 - abs( (sum(x) / (2*length(x))) - 0.5))


# plot a handful of plots

p <- list("X6.32560501", "X6.32560534" , "X6.32560670" ,"X6.32560695")
p <- list("X6.32560483", "X6.32560501", "X6.32560534" ,   "X6.32560631",   "X6.32560670" ,  "X6.32560695"  , "X6.32560711")
wrap_plot <- function(x) {
  vars <- sym(x)
  ggplot(dat_test_maf) +
    geom_point(aes(y = study_gxe, x = !!vars)) +
    theme_bw() +
    xlim(0,0.5)
}

pp <- lapply(p, wrap_plot)
do.call(grid.arrange, pp)



# script for GLM

# confirm results with GLM 

glm_func <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat_test, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat_test, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(dat_test[ , 2:22], function(x) glm_func_base(x))
intmodel <- map(dat_test[,2:22], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))








#-----------------------------------------------------------------------------#
# Plot MAFs ----
# GENERALIZE THIS WITH FUNCTION
#-----------------------------------------------------------------------------#
posthoc_df_maf <- dat %>%
  group_by(study_gxe) %>%
  summarise_at(vars(X5.40252294, X6.32560631), function(x) 0.5 - abs( (sum(x) / (2*length(x))) - 0.5))


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X5.40252294)) +
  theme_bw() +
  xlim(0,0.5)


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X6.32560631)) +
  theme_bw() +
  xlim(0,0.5)


#-----------------------------------------------------------------------------#
# Post Hoc Analyses ----
# prioritize SNPs that were significant only (any method)
#-----------------------------------------------------------------------------#

# Output GxE p values for LocusZoom
locuszoom_gxe <- gxe %>%
  mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)
write.table(locuszoom_gxe, file = "/media/work/tmp/LocusZoom_GxE_asp_ref_sex_age_pc3_studygxe_72820.txt", quote = F, row.names = F, sep = "\t")


sig1 <- filter(gxe, chiSqGxE > 30) #GxE
sig2 <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p <= wt) #two step method, G and EDGE same answer

# (manually note chr:bp and use local version of LocusZoom)




#-----------------------------------------------------------------------------#
# Output for genome browser track ----
# output top hit + LD > 0.8.. 
#-----------------------------------------------------------------------------#

# significant results
sig <- gxe %>% 
  mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         pval = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>% 
  filter(pval < 5e-8) 

gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p <= wt)


# first, make sure you run locuszoom to create list of LD snps
load("~/Dropbox/LocusZoom_GxE_asp_ref_sex_age_pc3_studygxe_72820_191002_5_40252294/chr5_39752294-40752294.Rdata")
chr5_40252294 <- metal %>% 
  filter(ldcut == "(0.8,1]")

load("~/Dropbox/LocusZoom_GxE_asp_ref_sex_age_pc3_studygxe_72820_191002_6_32560631/chr6_32060631-33060631.Rdata")
chr6_32560631 <- metal %>% 
  filter(ldcut == "(0.8,1]")

out_annovar <- rbind(chr5_40252294, chr6_32560631) %>% 
  mutate(SNP = paste0(chr, ":", pos_int)) %>% 
  dplyr::select(SNP, chr, pos_int) %>% unique() %>% 
  inner_join(gxe[, c("SNP", "Reference", "Alternate")], by = 'SNP') # get ref/alt information from results.. 

write.table(out_annovar[, c("chr", "pos_int", "pos_int", "Reference", "Alternate")], file = "/media/work/annovar/example/figi_asp_ref_sig.txt", quote = F, row.names = F, col.names = F)


# read annovar output, let's prepare output for genome browser custom track.. 
out_gb <- fread("/media/work/annovar/figi_asp_ref_sig_annot.hg19_multianno.csv") %>% 
  mutate(chrout = paste0("chr", Chr))

write.table(out_gb[, c("chrout", "Start", "End", "avsnp147")], file = "~/Dropbox/figi_asp_ref_sig_genomebrowser.txt", quote = F, row.names = F, col.names = F)


paste0("chr5:", 
       min(out_gb[which(out_gb$Chr == 5), "End"]),
       "-", 
       max(out_gb[which(out_gb$Chr == 5), "End"]))


paste0("chr6:", 
       min(out_gb[which(out_gb$Chr == 6), "End"]),
       "-", 
       max(out_gb[which(out_gb$Chr == 6), "End"]))



# just a record for vcftools (create vcf file of markers for CADD etc)
write.table(out_annovar[, c('SNP')], file = "~/Dropbox/figi_asp_ref_sig_vcftools.txt", quote = F, row.names = F, col.names = F, sep = '\t')

vcftools --gzvcf /media/work/FIGI/newfoundland_omniquad/newfoundland_omniquad_chr5.vcf.gz --indv s_1_Omniquad --snps ~/Dropbox/figi_asp_ref_sig_vcftools.txt --recode --out ~/figi_asp_ref_sig_vcfout





# get index positions for getSNPvalues
out_bd <- rbind(chr5_40252294, chr6_32560631) %>% 
  dplyr::mutate(ID = paste0(chr, ":", pos_int))
figi <- readRDS("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr5.rds")[[11]]
saveRDS(which(figi$SNPID %in% unique(out_bd$ID)), file = paste0("~/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr", 5, ".rds"), version = 2)


out_bd <- rbind(chr5_40252294, chr6_32560631) %>% 
  dplyr::mutate(ID = paste0(chr, ":", pos_int))
figi <- readRDS("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr6.rds")[[11]]
saveRDS(which(figi$SNPID %in% unique(out_bd$ID)), file = paste0("~/GetSNPValues_asp_ref_sex_age_pc3_studygxe_index_positions_chr", 6, ".rds"), version = 2)





# get SNP index positions (one file per chromosome)
for(chr in unique(out_bd$chr)) {
  figi <- readRDS(paste0("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr", chr, ".rds"))
  figi_snps <- figi[[11]] %>% 
    mutate(ID = paste0(Chromosome, ":", Location, ":", Reference, ":", Alternate))
  saveRDS(which(figi_snps$ID %in% unique(sig_all$ID)), file = paste0("files/GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", chr, ".rds"), version = 2)
}




#-------------------------------------#
# GxE
#-------------------------------------#



sig <- gxe %>% 
  mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         pval = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>% 
  filter(pval < 5e-8) 

# single hit on chromosome 6, let's grab this + SNPs in LD @ 0.5 based on LD structure among corect_oncoarray controls
chr6_32560631 <- c("6:32560631:C:T","6:32345891:T:G","6:32349086:T:C","6:32349772:C:T","6:32371915:C:T","6:32375973:C:A","6:32376348:G:A","6:32381374:G:C","6:32381443:A:T","6:32384801:A:T","6:32387809:T:C","6:32396930:G:A","6:32402778:T:G","6:32405026:G:A","6:32405821:T:C","6:32433067:A:G","6:32451822:C:T","6:32559825:A:G","6:32560695:G:C","6:32560839:C:T","6:32561201:A:G","6:32561334:C:G","6:32561465:C:T","6:32561466:A:G","6:32561638:C:T","6:32561681:G:A","6:32566398:G:A","6:32567256:C:T","6:32570311:T:G","6:32571122:C:T","6:32571961:T:C","6:32573415:A:G","6:32573562:T:C","6:32574603:G:A","6:32577222:G:A","6:32577380:A:G","6:32577973:T:G","6:32578590:A:G","6:32578632:C:T","6:32578772:C:A","6:32580617:T:A","6:32581973:G:A","6:32582189:A:C","6:32582650:C:G","6:32583099:A:G","6:32583357:A:T","6:32583529:G:A","6:32584693:C:G","6:32587350:G:A","6:32591588:A:G","6:32600101:G:C","6:32603798:A:T","6:32605078:A:G","6:32605609:G:A","6:32607853:A:T","6:32608077:T:C","6:32608269:G:A","6:32608478:T:C","6:32608998:T:C","6:32609545:C:T","6:32610401:G:A","6:32628712:C:T","6:32648987:C:T","6:32649088:C:G","6:32649126:A:G","6:32649161:C:T","6:32652196:T:C","6:32653263:A:G","6:32654714:T:C","6:32656947:G:A","6:32657710:A:G","6:32658335:T:C","6:32658665:G:A","6:32659099:T:C","6:32659839:G:A","6:32662607:G:A","6:32663308:T:C","6:32663431:T:C","6:32664126:A:G","6:32664990:A:T","6:32665255:T:A","6:32665319:T:G","6:32665640:G:T","6:32665759:G:T","6:32665909:C:A","6:32665912:C:A","6:32666635:C:T","6:32666651:G:T","6:32666660:C:T","6:32667107:T:C","6:32667850:C:G","6:32668411:G:A","6:32668831:T:C","6:32668846:G:T","6:32669013:T:C","6:32669454:A:C","6:32669483:A:G","6:32669761:C:G","6:32670548:G:A","6:32670912:T:C","6:32671057:G:A","6:32671086:C:T","6:32671247:C:T","6:32671332:G:C","6:32671412:C:T","6:32671755:G:A","6:32672624:T:C","6:32672641:A:C","6:32673099:T:A","6:32673385:C:G","6:32673574:A:T","6:32673931:A:G","6:32675523:G:C","6:32675634:C:G","6:32676159:C:G","6:32678491:T:C","6:32680299:G:A","6:32680620:G:T","6:32681483:T:C","6:32681530:T:C","6:32682019:A:C","6:32682137:A:G","6:32682429:C:T")

sig_chr6_32560631 <- gxe %>% 
  mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         pval = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>% 
  dplyr::filter(ID %in% chr6_32560631)

paste0("chr6:", min(sig_chr6_32560631$Location), "-", max(sig_chr6_32560631$Location))


#-------------------------------------#
# two-step chiSqG
#-------------------------------------#
gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p <= wt)

# two results: same locus as above, and chr5:40252294
chr5_40252294 <- c("5:40252294:C:T","5:40233337:G:A","5:40233662:T:G","5:40234224:A:G","5:40238192:C:T","5:40247015:G:C","5:40249579:G:A","5:40253503:C:T","5:40254195:A:G","5:40256326:A:G","5:40257771:T:C","5:40257953:C:T","5:40260962:G:A","5:40261340:C:T","5:40263604:A:G","5:40274458:A:G","5:40284614:A:G","5:40284625:T:G","5:40286497:C:T")

sig_chr5_40252294 <- gxe %>% 
  mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         pval = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>% 
  dplyr::filter(ID %in% chr5_40252294)

paste0("chr5:", min(sig_chr5_40252294$Location), "-", max(sig_chr5_40252294$Location))


#-------------------------------------#
# output csv file
#-------------------------------------#
out <- rbind(sig_chr5_40252294, sig_chr6_32560631) %>% 
  dplyr::select(-SNP) %>% 
  dplyr::rename(SNP = ID)

rsid_5 <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr5.rds") %>% 
  mutate(SNP = paste(`#CHROM`, POS, REF, ALT, sep = ":")) %>% 
  filter(SNP %in% out$ID)

rsid_6 <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr6.rds") %>% 
  mutate(SNP = paste(`#CHROM`, POS, REF, ALT, sep = ":")) %>% 
  filter(SNP %in% out$ID)

rsid <- rbind(rsid_5, rsid_6) %>% 
  dplyr::select(ID, SNP)

outf <- inner_join(out, rsid, 'SNP') %>% 
  mutate(chr = paste0("chr", Chromosome))

write.csv(outf, file = "~/Dropbox/test.csv", quote = F, row.names = F)



# output same information for annovar

outf_annovar <- inner_join(out, rsid, 'SNP') %>% 
  dplyr::mutate(start = Location,
                end = Location) %>% 
  dplyr::select(Chromosome, start, end, Reference, Alternate)
write.table(outf_annovar, file = "/media/work/annovar/example/figi_asp_ref_sig.txt", quote = F, row.names = F, col.names = F, sep = '\t')

# get a range of significant hits. 
## 5:40252294
chr <- 5
bp <- 40252294
paste(bp - 500000, bp + 500000)

sig_func <- gxe %>% 
  filter(Chromosome == chr & between(Location, bp-500000, bp+500000)) %>% 
  mutate(p = pchisq(chiSqGxE, df = 1, lower.tail = F),
         logp = -log10(p),
         start = Location - 1)

sig_func_range <- filter(sig_func, logp > 4)

paste0("chr", chr, ":", min(sig_func_range$Location), "-", max(sig_func_range$Location))

func_sig_annovar <- filter(func, logp > 4) %>% 
  mutate(start = Location) %>% 
  dplyr::select(Chromosome, start, Location, Location, Reference, Alternate)




# genome graph?
chr <- 5
bp <- 40252294
paste(bp - 500000, bp + 500000)

sig_func <- gxe %>% 
  mutate(p = pchisq(chiSqGxE, df = 1, lower.tail = F),
         logp = -log10(p),
         chr = paste0("chr", Chromosome)) %>% 
  filter(logp > 3) %>% 
  dplyr::select(chr, Location, logp)



write.csv(sig_func, file = "~/Dropbox/test.csv", quote = F, row.names = F)
