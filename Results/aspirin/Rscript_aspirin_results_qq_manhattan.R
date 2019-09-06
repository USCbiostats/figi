#=============================================================================#
# NSAIDS - aspirin results
# 06/01/2019
# 
# Filter by Weighted Average Rsq (until another solution is devised e.g. reimputation)
# Generate Plots
# Implement 2-step methods
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())


# annotations
# (right now, I only use this df to ensure that these markers stay in plot after filtering by P < 0.05)
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))


# Rsq Filter (vector of IDs to keep)

# (OLD)
# rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds") %>% 
#   filter(Rsq_avg > 0.8)

# Rsq estimate based on alt allele probability
# filtered maf > 0.001, Rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")


# results
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/aspirin/", full.names = T, pattern = "FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)

# gxe_noGWAS <- filter(gxe, !SNP %in% fh_annotations$SNP)

# output chiSqGxE results for LD clumping
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}

for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
    filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
}

for(chr in 1:22) {
  out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
    filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
    dplyr::select(SNP, P)
  write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
}




# LD Clump Results

## chiSqGxE
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump/", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqGxE_"), fread, stringsAsFactors = F))

gxe_clump_chiSqGxE <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)



## chiSqG
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqG_"), fread, stringsAsFactors = F))

gxe_clump_chiSqG <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)

# ----- testing ------ #
# # output results for LD clumping chiSqG (to compare manhattan plots)
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }

# gxe_chiSqG_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/NSAIDS/aspirin_190528/aspirin_chiSqG_ldclump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
# 
# gxe_chiSqG_ld <- gxe %>% 
#   filter(ID %in% gxe_chiSqG_ldclumped$SNP)



#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "aspirin"
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



### QUICK TESTS TO MAKE PLOTS MORE CONSISTENT
xxx <- sample(7267852, 100000, replace = F)
zzz <- gxe[xxx, ]

create_qqplot(zzz, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_DELETE")
create_manhattanplot(zzz, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_DELETE")


gxe_twostep <- format_twostep_data(dat = zzz, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_DELETE")


### Separate analysis

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_noGWAS, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "TEST")

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')




#-----------------------------------------------------------------------------#
# Two-Step after LD Clump ----
# chiSqGxE
#-----------------------------------------------------------------------------#

# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = '_clump_chiSqGxE')

# G|E 2-step Murcray ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_clump_chiSqGxE")

# EDGE 2-step Gauderman ----
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_clump_chiSqGxE")


gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p < wt)

zzz <- filter(fh_annotations, chr == "chr5")


# what proportion of twostep markers are in the annotation file..
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  separate(rsid, into = c('a', 'b', 'ref', 'alt'), sep = ":") %>% 
  mutate(ID = paste(Chr, Pos, ref, alt, sep = ":"))

wtf <- inner_join(gxe_twostep, fh_annotations, by = 'ID')

hhh <- filter(gxe_twostep, grepl("117630683", SNP))
hhh <- filter(fh_annotations, grepl("117630683", ID))


## Test, remove SNPs that are annotated...
yyy <- filter(fh_annotations, Colour == "red")

zzz <- filter(gxe_chiSqGxE_ld, !SNP %in% yyy$SNP)
zzz <- filter(gxe_chiSqGxE_ld, !SNP %in% fh_annotations$SNP)


gxe_twostep <- format_twostep_data(dat = zzz, 'chiSqG', 5, 0.05)

wtf <- inner_join(gxe_twostep, fh_annotations, by = "SNP")

create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = '_LDclump_NOGWASHITS')



#-----------------------------------------------------------------------------#
# Two-Step after LD Clump ----
# chiSqG
#-----------------------------------------------------------------------------#
plot_exposure <- "aspirin"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")

# Marginal G Results ----
create_qqplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_clump_chiSqG")
create_manhattanplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_clump_chiSqG")

# GxE results ----
create_qqplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1, filename_suffix = "_clump_chiSqG")


# D|G 2-step Kooperberg ----
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_clump_chiSqG")



# D|G 2-step Kooperberg ----
# Try removing the GWAS hits from the two-step plot... 
x <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, c(1:12, 18:21)] %>%
  mutate(SNP = paste0(CHROM, ":", POSITION)) 

tmp <- gxe_clump_chiSqG %>% 
  filter(!SNP %in% x$SNP)

gxe_twostep <- format_twostep_data(dat = tmp, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_DELETEDELETEDELETE")

xx <- filter(gxe_twostep, step2p <= wt)

#-----------------------------------------------------------------------------#
# Extract dosage values in significant loci from binary dosage files ----
# (needs hand holding - look at manhattan plots to identify significant hits)
# (for GetSNPValues - need sample list and SNP index positions)
#-----------------------------------------------------------------------------#

# Significant Results
sig1 <- filter(gxe, chiSqGxE > 30) #GxE
sig2 <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p <= wt) #two step method, G and EDGE same answer

sig_all <- bind_rows(sig1, sig2)

# get SNP index positions (one file per chromosome)
for(chr in unique(sig_all$Chromosome)) {
  figi <- readRDS(paste0("/home/rak/data/FIGI_BDose_IndexFiles/FIGI_chr", chr, ".rds"))
  figi_snps <- figi[[11]] %>% 
    mutate(ID = paste0(Chromosome, ":", Location, ":", Reference, ":", Alternate))
  saveRDS(which(figi_snps$ID %in% unique(sig_all$ID)), file = paste0("files/GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", chr, ".rds"), version = 2)
}

# write sample list (vector)
vcfid <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")[, 'vcfid']
saveRDS(vcfid, file = "files/GetSNPValues_aspirin_sex_age_pc3_studygxe_vcfid.rds", version = 2)





# ... 
# ... 
# ...


# confirm results with GLM 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/aspirin/dosage/", pattern = "GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds"))

glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(dat[ , 2:4], function(x) glm_func_base(x))
intmodel <- map(dat[,2:4], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))



# for X6.32560631, might be driven by UKB1.. 
dat <- do.call(inner_join, lapply(list.files("~/data/Results/aspirin/dosage/", pattern = "GetSNPValues_aspirin_sex_age_pc3_studygxe_index_positions_chr", full.names = T), function(x) data.frame(readRDS(x)) %>% rownames_to_column("vcfid"))) %>% 
  merge(readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds"))

dat <- filter(dat, study_gxe != "UKB_1")

glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = dat, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

basemodel <- map(dat[ , 2:4], function(x) glm_func_base(x))
intmodel <- map(dat[,2:4], function(x) glm_func(x))
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
  summarise_at(vars(X5.40252294:X6.12577203), function(x) 0.5 - abs( (sum(x) / (2*length(x))) - 0.5))


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X5.40252294)) +
  theme_bw() +
  xlim(0,0.5)


ggplot(posthoc_df_maf) +
  geom_point(aes(y = study_gxe, x = X6.12577203)) +
  theme_bw() +
  xlim(0,0.5)



#-----------------------------------------------------------------------------#
# Followup interesting hits ----
#-----------------------------------------------------------------------------#
gxe_twostep_ld <- format_2step_data(data = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>% 
  filter(step2p < wt)



#-----------------------------------------------------------------------------#
# locuszoom ------
# there are several regions, we should focus on one at 
# a time
#-----------------------------------------------------------------------------#

# Marginal G 

# GxE
locuszoom_gxe <- gxe %>%
  mutate(`P-value` = pchisq(chiSqGxE, df = 1, lower.tail = F),
         MarkerName = paste0("chr", SNP)) %>%
  dplyr::select(MarkerName, `P-value`)
write.table(locuszoom_gxe, file = "/media/work/tmp/LocusZoom_GxE_aspirin_sex_age_pc3_studygxe_66485.txt", quote = F, row.names = F, sep = "\t")


#-----------------------------------------------------------------------------#
# Extract top hits dosages from binarydosage 
#-----------------------------------------------------------------------------#
# first take care of sample list (vcfid)
# might want to move that to the top of the file.. 
sample_list <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")[, 'vcfid']
saveRDS(sample_list, file = paste0(write_binarydosage_vcfid_filename(), ".rds"), version = 2)


# get SNPs of interest based on plots above - in this case, two-step data.table
# (don't forget in this case the sig results was from clumped file)
# do it by chromosome for now
gxe_twostep <- format_2step_data(data = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) 
gxe_twostep_sig <- gxe_twostep[step2p < wt, ]

get_binarydosage_index(gxe_twostep_sig$ID, 5)


#-----------------------------------------------------------------------------#
# Extract top hits dosages from binarydosage 
#-----------------------------------------------------------------------------#
covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")
dosages <- data.frame(readRDS("files/GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr5_out.rds")) %>% 
  rownames_to_column(var = 'vcfid')

dosages <- data.frame(readRDS("/home/rak/Dropbox/FIGI/Results/aspirin_190605/files/GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr5_out.rds")) %>% 
  rownames_to_column(var = 'vcfid')


# need to merge in dosage info to the EpiData
posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid')


# then apply gxe only to dosage columns, with the same covariates used for gxescan



# then compile information, create forest plot of betas CIs. 






















