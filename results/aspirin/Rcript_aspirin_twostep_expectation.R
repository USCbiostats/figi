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

fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv")


# fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv")

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
gxe_chiSqGxE_ldclumped <- do.call(rbind, lapply(list.files("~/data/Results/asp_ref/clump/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))

gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% gxe_chiSqGxE_ldclumped$SNP)

rm(gxe_chiSqGxE_ldclumped)



#-----------------------------------------------------------------------------#
# Expectation based weighted hypothesis testing ----
#-----------------------------------------------------------------------------#

# first, let's remove the GWAS top hits

tmp <- filter(gxe, !SNP %in% fh_annotations$SNP)


# first need to process the results, create bins, assign bins 

x <- data.table(tmp)[, step1p_g := pchisq(chiSqG, df = 1, lower.tail = F)
                     ][ 
                       , step1p_ge := pchisq(chiSqGE, df = 1, lower.tail = F)
                       ][
                         , step1p_edge := pchisq(chiSqEDGE, df = 2, lower.tail = F)
                       ][
                       , step2p := pchisq(chiSqGxE, df = 1, lower.tail = F)
                       ][
                         , logstep2p := -log10(step2p)
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
# start with marginal G step1 statistic

xx <- x[order(step1p_g)]

# let's keep same size bins.. assuming total of 1 million
m = 1000000
nbins = floor(log2(m/5)+1) # 17 bins
sizeBin = c(5 * 2^(0:(nbins-2)), m - 5 * (2^(nbins-1) - 1) ) #bins + the last one
alphaBin <- sizeBin / 1000000
alphaBinCut <- c(-Inf, alphaBin, Inf)

rk.pv <- c(1:m)
# grp <- ceiling(log(rk.pv/5+1,base=2))
grp <- as.numeric(cut(xx$step1p_g, breaks = alphaBinCut ))      
grp_help <- table(grp)
grp_help2 <- as.vector(as.integer(grp_help))

#rep(alphaBin, grp_help2)


wtf <- xx[ ,grp:=grp]


setkey(wtf,grp)
for(i in 1:max(grp))
{
  wtf[J(i),wt:=0.05*2^(-i)/sizeBin[i]] # data.table syntax, create threshold value
}


# so then what would be significant here...

sig_results <- filter(wtf, step2p < wt)


# plot the results
endpointsBin = cumsum(sizeBin) # endpoints of the bins
binsToPlot = 10
# create list and vars for plotting
min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
last.sig = alphaBin[10]

# create list where each component contains a bin
# log transform bin alpha value
# create 'x', normalized position information for each bin
glist<-list()
for(i in 1:10){
  t <- wtf[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}

# CREATE PLOT
head(glist[[1]]) # for reference


# png(write_weighted_plot_filename(deparse(substitute(dat)), statistic), width = 1280, height = 720)
color <- rep(c("blue","olivedrab4"),100)
plot(glist[[1]][,x], glist[[1]][,logstep2p],
     col = "blue1",
     xlab="Bin number for step1 p value",
     ylab="-log10(step2 chiSqGxE p value)",
     xlim=c(0,binsToPlot),
     ylim=c(0,min.p),
     axes=F, pch=19,
     cex.main = 1.6,
     cex.axis = 1.3,
     cex.lab = 1.3,
     cex.sub = 1.3,
     cex = 1)
lines(glist[[1]][,x], glist[[1]][,ref],
      col = "black",lwd=2)

# the rest of the points ...  =|
# (adding to current plot..)
for(i in 2:binsToPlot){
  points(glist[[i]][,x], glist[[i]][,logstep2p],
         col = color[i], pch = 19,
         cex.main = 1.6,
         cex.axis = 1.3,
         cex.lab = 1.3,
         cex.sub = 1.3,
         cex = 1)
  lines(glist[[i]][,x], glist[[i]][,ref],
        col = "black",lwd = 2)
}






#### let's try the above but using the edge statistic

xx <- x[order(step1p_edge)]

# let's keep same size bins.. assuming total of 1 million
m = 1000000
nbins = floor(log2(m/5)+1) # 17 bins
sizeBin = c(5 * 2^(0:(nbins-2)), m - 5 * (2^(nbins-1) - 1) ) #bins + the last one
alphaBin <- sizeBin / 1000000
alphaBinCut <- c(-Inf, alphaBin, Inf)

rk.pv <- c(1:m)
# grp <- ceiling(log(rk.pv/5+1,base=2))
grp <- as.numeric(cut(xx$step1p_edge, breaks = alphaBinCut ))      
grp_help <- table(grp)
grp_help2 <- as.vector(as.integer(grp_help))

#rep(alphaBin, grp_help2)


wtf <- xx[ ,grp:=grp]


setkey(wtf,grp)
for(i in 1:max(grp))
{
  wtf[J(i),wt:=0.05*2^(-i)/sizeBin[i]] # data.table syntax, create threshold value
}


# so then what would be significant here...

sig_results <- filter(wtf, step2p < wt)


# plot the results
endpointsBin = cumsum(sizeBin) # endpoints of the bins
binsToPlot = 10
# create list and vars for plotting
min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
last.sig = alphaBin[10]

# create list where each component contains a bin
# log transform bin alpha value
# create 'x', normalized position information for each bin
glist<-list()
for(i in 1:10){
  t <- wtf[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}

# CREATE PLOT
head(glist[[1]]) # for reference


# png(write_weighted_plot_filename(deparse(substitute(dat)), statistic), width = 1280, height = 720)
color <- rep(c("blue","olivedrab4"),100)
plot(glist[[1]][,x], glist[[1]][,logstep2p],
     col = "blue1",
     xlab="Bin number for step1 p value",
     ylab="-log10(step2 chiSqGxE p value)",
     xlim=c(0,binsToPlot),
     ylim=c(0,min.p),
     axes=F, pch=19,
     cex.main = 1.6,
     cex.axis = 1.3,
     cex.lab = 1.3,
     cex.sub = 1.3,
     cex = 1)
lines(glist[[1]][,x], glist[[1]][,ref],
      col = "black",lwd=2)

# the rest of the points ...  =|
# (adding to current plot..)
for(i in 2:binsToPlot){
  points(glist[[i]][,x], glist[[i]][,logstep2p],
         col = color[i], pch = 19,
         cex.main = 1.6,
         cex.axis = 1.3,
         cex.lab = 1.3,
         cex.sub = 1.3,
         cex = 1)
  lines(glist[[i]][,x], glist[[i]][,ref],
        col = "black",lwd = 2)
}












#-----------------------------------------------------------------------------#
# GWAS hits adjusted at 0.05/140
#-----------------------------------------------------------------------------#


fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv")


hits <- filter(gxe, ID %in% fh_annotations$SNP_B ) %>% 
  mutate(pval =  pchisq(chiSqGxE, df = 1, lower.tail = F),
         sig = ifelse(pval <= 0.05/140, 1, 0)) %>% 
  filter(sig == 1)




fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

hits <- filter(gxe, SNP %in% fh_annotations$SNP ) %>% 
  mutate(pval =  pchisq(chiSqGxE, df = 1, lower.tail = F),
         sig = ifelse(pval <= 0.05/140, 1, 0)) %>% 
  filter(sig == 1)


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




