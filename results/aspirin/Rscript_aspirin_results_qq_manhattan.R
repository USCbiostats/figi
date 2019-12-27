#=============================================================================#
# NSAIDS - aspirin results
# 06/01/2019
# 09/20/2019
# 
# Notes:
# - Rsq filter > 0.8, estimated on alternate allele probability variance
#   (see GxEScanR fork on github)
# - Generate results plots
#   all data
#   ld clumped
#   expectation based (two-step methods)
#
# Annotations are based on 140 top hits, with LD SNPs calculated from
# corect_oncoarray control samples (N ~ 12,000)
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
rm(list = ls())

# LD based annotations
# fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
#   mutate(SNP = paste(Chr, Pos, sep = ":"))
fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv")


# Rsq Filter (vector of IDs to keep)
# maf > 0.001, Rsq >= 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")


# RESULTS
gxe_all <- do.call(rbind, lapply(list.files(path = "~/data/Results/aspirin/", full.names = T, pattern = "FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_chr"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = ":"),
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE) %>% 
  filter(!duplicated(ID)) # small issue with gxescan package

gxe <- gxe_all %>% 
  filter(ID %in% rsq_filter$id)

# # --------------------------------------- #
# # output chiSqGxE results for LD clumping
# calculate_pval <- function(data, statistic, df) {
#   data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
#   data
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqGxE', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_ldclump_aspirin_age_sex_pc3_studygxe_72269_chiSqGxE_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }
# 
# for(chr in 1:22) {
#   out <- calculate_pval(gxe, 'chiSqG', df = 1) %>%
#     filter(Chromosome == chr) %>% mutate(SNP = ID) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/tmp/Plink_aspirin_ldclump_chiSqG_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }
# # --------------------------------------- #



# LD Clump Results

## chiSqGxE
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump/", full.names = T, pattern = "Plink_ldclump_aspirin_age_sex_pc3_studygxe_72269_chiSqGxE_chr"), fread, stringsAsFactors = F))
gxe_clump_chiSqGxE <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)

## chiSqGxE - Controls only
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump_controls/", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqGxE_controls"), fread, stringsAsFactors = F))
gxe_clump_chiSqGxE <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)

## chiSqG
tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump", full.names = T, pattern = "Plink_aspirin_ldclump_chiSqG_"), fread, stringsAsFactors = F))
gxe_clump_chiSqG <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)


#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots ----
#-----------------------------------------------------------------------------#
plot_exposure <- "aspirin"
plot_covariates <- c("age_ref_imp", "sex", "study_gxe", "PC1", "PC2", "PC3")

# Marginal G Results 
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = test)

# GxE results
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)

# 2DF results
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq2df', df = 2)

# 3DF results
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSq3df', df = 3)

# GE, Case, Control
create_qqplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqGE', df = 1)
create_qqplot_ge(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqControl', df = 1)
create_qqplot_ge(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)
create_manhattanplot(gxe, plot_exposure, plot_covariates, stat = 'chiSqCase', df = 1)

# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')




#-----------------------------------------------------------------------------#
# Two-step method w/ LD clumping chiSqGxE statistic ----
#-----------------------------------------------------------------------------#
# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = '_clump_chiSqGxE')

# G|E 2-step Murcray
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE', filename_suffix = "_clump_chiSqGxE")

# EDGE 2-step Gauderman
gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqGxE, 'chiSqEDGE', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE', filename_suffix = "_clump_chiSqGxE")


#-----------------------------------------------------------------------------#
# Two-step method w/ LD clumping chiSqG statistic  ----
#-----------------------------------------------------------------------------#
# Marginal G Results
create_qqplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_clump_chiSqG")
create_manhattanplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_clump_chiSqG")

# GxE results
create_qqplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1)
create_manhattanplot(gxe_clump_chiSqG, plot_exposure, plot_covariates, stat = 'chiSqGxE', df = 1, filename_suffix = "_clump_chiSqG")

# D|G 2-step Kooperberg
gxe_twostep <- format_twostep_data(dat = gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG', filename_suffix = "_clump_chiSqG")



#-----------------------------------------------------------------------------#
# Two-step method w/ bins+adjustment based on expectation under null  ----
#-----------------------------------------------------------------------------#

# in this specific instance, the only signficant result is for step1p_g (no edge, ge)

# perform analysis separately for top hits vs everything else
gxe_top <- filter(gxe, ID %in% fh_annotations$SNP_B) %>% 
  mutate(pGxE = pchisq(chiSqGxE, df = 1, lower.tail = F))
create_manhattanplot(gxe_top, plot_exposure, plot_covariates, stat = 'chiSqG', df = 1, filename_suffix = "_EXPECTATION_GWAShits_chiSqG")
manhattan(gxe_top, chr = "Chromosome", bp = "Location", p = "pGxE", col = c("blue4", "orange3"), suggestiveline = F, genomewideline = -log10(0.05/140), main = "GWAS Hits GxE Significant at P = 0.05/140")

gxe_top_filter <- filter(gxe_top, pGxE < (0.05/140))


# start with ordering by step1p_g (marginal)
gxe_rest <- data.table(gxe)[
  !ID %in% fh_annotations$SNP_B][
  , step1p_g := pchisq(chiSqG, df = 1, lower.tail = F)][
  , step1p_ge := pchisq(chiSqGE, df = 1, lower.tail = F)][
  , step1p_edge := pchisq(chiSqEDGE, df = 2, lower.tail = F)][
  , step2p := pchisq(chiSqGxE, df = 1, lower.tail = F)][
  , logstep2p := -log10(step2p)][
  , MapInfo := Location][order(step1p_g)]


# assume total NSNP = 1,000,000
# base cutoffs and correction on expected distribution of step1 p-values
m = 1000000
nbins = floor(log2(m/5)+1) # 17 bins
sizeBin = c(5 * 2^(0:(nbins-2)), m - 5 * (2^(nbins-1) - 1) ) #bins + the last one
alphaBin <- sizeBin / 1000000
alphaBinCut <- c(-Inf, alphaBin, Inf)

rk.pv <- c(1:m)
grp <- as.numeric(cut(gxe_rest$step1p_g, breaks = alphaBinCut ))      
grp_help <- table(grp)
grp_help2 <- as.vector(as.integer(grp_help))

gxe_rest <- gxe_rest[ ,grp:=grp]
setkey(gxe_rest,grp)
for(i in 1:max(grp))
{
  gxe_rest[J(i),wt:=0.05*2^(-i)/sizeBin[i]] # data.table syntax, create threshold value
}

sig_results <- filter(gxe_rest, step2p < wt)


# Plot

endpointsBin = cumsum(sizeBin) # endpoints of the bins
binsToPlot = 10
min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
last.sig = alphaBin[10]

# create list where each component contains a bin
# log transform bin alpha value
# create 'x', normalized position information for each bin
glist<-list()
for(i in 1:10){
  t <- gxe_rest[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}

# CREATE PLOT
head(glist[[1]]) # for reference


png("figures/TwoStep_WeightedHypothesis_gxe_twostep_aspirin_chiSqG_age_ref_imp_sex_study_gxe_PC1_PC2_PC3_N_72269_EXPECTATION_chiSqGxE.png", width = 1280, height = 720)
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

# plot remaining points
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

axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1)
title(main = "D|G 2-step Procedure Results (expectation based cutoffs/correction)", sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.6, cex.sub = 1.3)
dev.off()



#-----------------------------------------------------------------------------#
# SUGGESTIVE GxE  ----
#-----------------------------------------------------------------------------#

gxe_suggest <- gxe %>% 
  mutate(pGxE = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>% 
  filter(pGxE < 1e-5)



#-----------------------------------------------------------------------------#
# Extract dosage values in significant loci from binary dosage files ----
# output rds files with index position vectors
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





#-----------------------------------------------------------------------------#
# Extract SNP ranges for functional DB lookups ----
#-----------------------------------------------------------------------------#

# start with clumped findings (chr 5 40252294 +- 500kb)

func <- gxe %>% 
  filter(Chromosome == 5 & between(Location, 40252294-500000, 40252294+500000)) %>% 
  mutate(p = pchisq(chiSqGxE, df = 1, lower.tail = F),
         logp = -log10(p),
         start = Location - 1)

func_sig <- filter(func, logp > 4) %>% 
  dplyr::select(Chromosome, start, Location)


func_sig_annovar <- filter(func, logp > 4) %>% 
  mutate(start = Location) %>% 
  dplyr::select(Chromosome, start, Location, Location, Reference, Alternate)

write.table(func_sig_annovar, file = "~/annovar/example/figi_aspirin_wthyp_clump_sig.txt", quote =  F, row.names = F, col.names = F, sep = '\t')


# let's try a more targeted region (essentially all markers in betweeen the significant hits after filtering chiSqGxE -log10 p value > 4)
# chr5:40234224-40286497

func <- gxe %>% 
  filter(Chromosome == 5 & between(Location, 40234224, 40286497)) %>% 
  mutate(p = pchisq(chiSqGxE, df = 1, lower.tail = F),
         logp = -log10(p),
         start = Location - 1)

func_sig_annovar <- func %>% 
  mutate(start = Location) %>% 
  dplyr::select(Chromosome, start, Location, Location, Reference, Alternate)


write.table(func_sig_annovar, file = "~/annovar/example/figi_aspirin_wthyp_clump_sig.txt", quote =  F, row.names = F, col.names = F, sep = '\t')

