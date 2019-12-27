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

# LD Clump Results
tmp <- do.call(rbind, lapply(list.files("~/data/Results/asp_ref/clump_controls/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)




#-----------------------------------------------------------------------------#
# Expectation based weighted hypothesis testing ----
#-----------------------------------------------------------------------------#

# first, let's remove the GWAS top hits

new_annotation <- fread("~/data/Annotations/temp_annotation_ver2.txt") %>% 
  mutate(SNP = paste0(Chr, ":", Pos))

new_annotation <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste0(Chr, ":", Pos))

tmp <- filter(gxe, !SNP %in% new_annotation$SNP)
create_manhattanplot(tmp, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_NO_GWAS_HITS_OLD")


# for now, let's remove anything < 5e-8

tmp <- gxe %>% 
  mutate(p = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
  filter(p > 5e-8)

create_manhattanplot(tmp, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_NO_GWAS_HITS_SIG")


# first need to process the results, create bins, assign bins 

x <- data.table(tmp)[, step1p := pchisq(chiSqG, df = 1, lower.tail = F)
                     ][
                     , step2p := pchisq(chiSqGxE, df = 1, lower.tail = F)
                     ][
                     , logstep2p := -log10(step2p)
                     ][
                     order(step1p)
                     ][
                     , MapInfo := Location
                      ]

m = nrow(x)
nbins = floor(log2(m/5+1))
sizeBin = c(5 * 2^(0:(nbins)))
sizeBin = c(5 * 2^(0:(nbins-2)), m - 5 * (2^(nbins-1) - 1) )
endpointsBin = cumsum(sizeBin)
alphaBin = 0.05 * 2 ^ -(1:nbins) / sizeBin 

rk.pv <- c(1:m)
grp <- ceiling(log(rk.pv/5+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
x[,grp:=grp] # assigning group to the p values..
setkey(x,grp)
for(i in 1:max(grp))
{
  x[J(i),wt:=0.05*2^(-i)/nrow(x[J(i)])] # data.table syntax, create threshold value
}

# confirm that in the original results, after removing top GWAS there are no sig results
sig_results <- filter(x , step2p < wt) # fat zero





#-------------------------------------#
# now, let's try different bin assignment
# supposing 1 million SNPs, adjustment levels in expectation

# let's keep same size bins.. assuming total of 1 million
m = 1000000
nbins = floor(log2(m/5)+1) # 17 bins
sizeBin = c(5 * 2^(0:(nbins-2)), m - 5 * (2^(nbins-1) - 1) ) #bins + the last one
alphaBin <- sizeBin / 1000000
alphaBinCut <- c(-Inf, alphaBin, Inf)

rk.pv <- c(1:m)
# grp <- ceiling(log(rk.pv/5+1,base=2))
grp <- as.numeric(cut(x$step1p, breaks = alphaBinCut ))      
grp_help <- table(grp)
grp_help2 <- as.vector(as.integer(grp_help))

#rep(alphaBin, grp_help2)


y <- x[ ,grp:=grp]


setkey(y,grp)
for(i in 1:max(grp))
{
  y[J(i),wt:=0.05*2^(-i)/sizeBin[i]] # data.table syntax, create threshold value
}


# so then what would be significant here...

sig_results <- filter(y, step2p < wt)






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


x <- filter(gxe_twostep, step2p <= wt)


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




