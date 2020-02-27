#=============================================================================#
# Recreating some plots for the gecco presentation
#
# perhaps incorporate these changes to the rest of the results reports 
# (makes manhattan and two-step plots much more legible)
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
library(RColorBrewer)
rm(list = ls())

exposure = 'asp_ref'
filename = 'FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan'
covars = c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3')
annotation = "gwas_140_ld_annotation_new.txt"

setwd(paste0("~/Dropbox/FIGI/Results/", exposure))
gxe <- readRDS(paste0("~/data/results/", exposure, "/processed/", filename, "_results.rds"))


# ---- manhattan plots ----
create_manhattanplot <- function(dat, exposure, covars, stat, df, annotation_file, filename_suffix = "") {
  
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])
  total <- cases + controls
  
  dat <- dat %>%
    mutate(P = calculate_pval(dat[,stat], df)) %>%
    filter(!(P > 0.01)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(SNP, CHR, BP, P)
  
  write.table(dat, file = paste0("/media/work/tmp/manhattan_", stat, "_", exposure, "_", paste0(covars, collapse = "_"), filename_suffix), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0("~/Dropbox/FIGI/Results/", exposure, "/plots")
  ecf2 <- "SNP;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- paste0("/media/work/tmp/manhattan_", stat, "_", exposure, "_", paste0(covars, collapse = "_"), filename_suffix)
  ecf_file_name <- paste0("~/Dropbox/FIGI/Results/", exposure, "/files/EasyStrata_", stat, "_", exposure, "_", paste0(covars, collapse = "_"), filename_suffix, ".ecf")
  
  cat(paste0("DEFINE	--pathOut ", ecf1, "
      --acolIn ", ecf2, "
      --acolInClasses ", ecf3, "
      --strMissing NA
      --strSeparator TAB

      EASYIN	--fileIn ", ecf4, "

      START EASYX

      ################
      ## MHplot
      ################

      MHPLOT
      --colMHPlot P
      --colInChr CHR
      --colInPos BP
      --numWidth 1280
      --numHeight 720
      #--astrDefaultColourChr gray51;gray66
      --astrDefaultColourChr gray70;gray80
      --blnYAxisBreak 1
      --numYAxisBreak 22
      #--numPvalOffset 0.01
      # Annotation
      --fileAnnot /home/rak/data/Annotations/", annotation_file, "
      --numAnnotPosLim 1
      # Horizontal lines
      --anumAddPvalLine 5e-6;5e-8
      --anumAddPvalLineLty 6;6
      --astrAddPvalLineCol blue;red
      # Other Graphical Params
      --anumParMar 7;5;7;3
      --numDefaultSymbol 1
      --numDefaultCex 0.6
      --numCexAxis 1.7
      --numCexLab 1.7
      --arcdSymbolCrit  (P>5e-8 & P<5e-6);P<5e-8
      --anumSymbol 20;19
      --arcdColourCrit (P>5e-8 & P<5e-6);P<5e-8
      --astrColour gray55;black
      --arcdCexCrit (P>5e-8 & P<5e-6);P<5e-8
      --anumCex 1.6;1.3

      STOP EASYX"), file = ecf_file_name, append = F)
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}


#create_manhattanplot(gxe_chr8, exposure, covars, 'chiSqG', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_chr8")


# create manhattan plots to replace the ones in the presentation
create_manhattanplot(gxe, exposure, covars, 'chiSqG', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation")
create_manhattanplot(gxe, exposure, covars, 'chiSqGxE', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation")
create_manhattanplot(gxe, exposure, covars, 'chiSq2df', 2, annotation_file = annotation, filename_suffix = "_GECCO_Presentation")
create_manhattanplot(gxe, exposure, covars, 'chiSq3df', 3, annotation_file = annotation, filename_suffix = "_GECCO_Presentation")


snps_to_exclude <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")
gxe_nogwas <- gxe %>%
  mutate(chr_bp = paste0(Chromosome, ":", Location)) %>%
  dplyr::filter(!chr_bp %in% snps_to_exclude$SNP) %>% 
  mutate(chiSqG_pval = pchisq(chiSqG, 1, lower.tail = F),
         chiSq2df_pval = pchisq(chiSq2df, 2, lower.tail = F))

create_manhattanplot(gxe_nogwas, exposure, covars, 'chiSqG', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_NOGWAS")
create_manhattanplot(gxe_nogwas, exposure, covars, 'chiSq2df', 2, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_NOGWAS")
create_manhattanplot(gxe_nogwas, exposure, covars, 'chiSq3df', 3, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_NOGWAS")



# generate locuszoom plots (for the 'novel' 2df peaks)
gxe_chr18 <- filter(gxe, Chromosome == 18) %>% 
  mutate(chiSqG_pval = pchisq(chiSqG, 1, lower.tail = F),
         chiSq2df_pval = pchisq(chiSq2df, 2, lower.tail = F))

loc <- filter(gxe_chr18, chiSqG_pval < 5e-8) %>% arrange(desc(chiSqG))
head(loc)
loc2 <- filter(gxe_nogwas, chiSq2df_pval < 5e-8) %>% arrange(desc(chiSq2df))
head(loc2)



gxe_chr8 <- filter(gxe, Chromosome == 8) %>% 
  mutate(chiSqG_pval = pchisq(chiSqG, 1, lower.tail = F),
         chiSq2df_pval = pchisq(chiSq2df, 2, lower.tail = F))

create_manhattanplot(gxe_chr8, exposure, covars, 'chiSqG', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_chr8")


loc <- filter(gxe_chr8, chiSqG_pval < 5e-8) %>% arrange(desc(chiSqG))
head(loc)
loc2 <- filter(gxe_nogwas, chiSq2df_pval < 5e-8) %>% arrange(desc(chiSq2df))
head(loc2)


gxe_chr1 <- filter(gxe, Chromosome == 1) %>% 
  mutate(chiSqG_pval = pchisq(chiSqG, 1, lower.tail = F),
         chiSq2df_pval = pchisq(chiSq2df, 2, lower.tail = F))

create_manhattanplot(gxe_chr1, exposure, covars, 'chiSqG', 1, annotation_file = annotation, filename_suffix = "_GECCO_Presentation_chr1")

loc <- filter(gxe_chr1, chiSqG_pval < 5e-6) %>% arrange(desc(chiSqG))
head(loc)
loc2 <- filter(gxe_nogwas, chiSq2df_pval < 5e-8) %>% arrange(desc(chiSq2df))
head(loc2)







# two step plots -----

format_twostep_data <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(dat, stats_step1, df=1) {
    data.table(dat)[, step1p := pchisq(dat[, eval(quote(stats_step1))], df = df, lower.tail = F)
                    ][
                      , step2p := pchisq(dat[,'chiSqGxE'],  df = 1, lower.tail = F)
                      ][
                        , y := -log10(step2p)
                        ][
                          order(step1p)
                          ][
                            , MapInfo := Location
                            ]
  }
  
  if(stats_step1 == 'chiSqEDGE') {
    pv <- create_pval_info(dat, stats_step1, df = 2)
  } else {
    pv <- create_pval_info(dat, stats_step1, df = 1)
  }
  
  
  # format output for plotting..
  m = nrow(pv)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  # add grp, wt (p value threshold for each bin), rename p value to 'y'
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
  pv[,grp:=grp] # assigning group to the p values..
  setkey(pv,grp)
  for(i in 1:max(grp))
  {
    pv[J(i),wt:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
  }
  
  # return the data.table
  return(pv)
}


create_twostep_weighted_plot <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "") {
  
  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls
  
  # some plot title and file names
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  # write_twostep_weightedHT_plot_title("chiSqG", exposure, covars, total)
  
  # bin information
  m = nrow(dat)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} # add +1 bin if condition met
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha level for each bbin 1,2, ... N bin tested
  
  # create list and vars for plotting
  min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]
  
  # create list where each component contains a bin
  # log transform bin alpha value
  # create 'x', normalized position information for each bin
  glist<-list()
  for(i in 1:binsToPlot){
    t <- dat[J(i)][order(Chromosome, Location)]
    t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
    t[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    t[, x := 0.8*((t[,mapinfo_tmp]-min(t[,mapinfo_tmp])) / (max(t[,mapinfo_tmp])-min(t[,mapinfo_tmp]))) + 0.1 + i - 1]
    glist[[i]]<-t[order(Chromosome, Location)]
    rm(t)
  }
  
  significant_hits <- dat[step2p <= wt]
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(paste0("~/Dropbox/FIGI/Results/", exposure, "/plots/twostep_wht_", statistic, "_", exposure, "_", paste0(covars, collapse = "_"), filename_suffix, ".png"), height = 720, width = 1280)
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  plot(glist[[1]][,x], glist[[1]][,y],
       col = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], '#E41A1C','#377EB8'),
       pch = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 19, 20),
       cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,min.p),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  # cex = 1.2)
  lines(glist[[1]][,x], glist[[1]][,ref], col = "black", lwd=1)
  
  # the rest of the points
  for(i in 2:binsToPlot){
    points(glist[[i]][,x], glist[[i]][,y],
           col = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, 19, 20),
           cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    # cex = 1.2)
    lines(glist[[i]][,x], glist[[i]][,ref], col = "black",lwd = 1)
  }
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
  
}



wtf <- format_twostep_data(gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(wtf, exposure, covars, 5, 0.05, 10, 'chiSqG', filename_suffix = "_GECCO_Presentation")



tmp <- do.call(rbind, lapply(list.files(paste0("~/data/results/", exposure, "/clump/"), full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>%
  filter(SNP %in% tmp$SNP)

wtf2 <- format_twostep_data(gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot(wtf2, exposure, covars, 5, 0.05, 10, 'chiSqG', filename_suffix = "_GECCO_Presentation_ld_clump_chiSqGxE")





# ---- two-step plots expectation ----

format_twostep_data_expectation <- function(dat, stats_step1, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(dat, stats_step1, df=1) {
    data.table(dat)[, step1p := pchisq(dat[, eval(quote(stats_step1))], df = df, lower.tail = F)
                    ][
                      , step2p := pchisq(dat[,'chiSqGxE'],  df = 1, lower.tail = F)
                      ][
                        , logstep2p := -log10(step2p)
                        ][
                          order(step1p)
                          ][
                            , MapInfo := Location
                            ]
  }
  
  if(stats_step1 == 'chiSqEDGE') {
    pv <- create_pval_info(dat, stats_step1, df = 2)
  } else {
    pv <- create_pval_info(dat, stats_step1, df = 1)
  }
  
  # format output for plotting..
  # m = nrow(pv)
  m = 1000000
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  
  ## sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin_step1 = endpointsBin/1000000
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  
  # expectation based grouping of step 1 pvalues
  alphaBinCut <- c(-Inf, alphaBin_step1, Inf)
  pv[ , grp:=as.numeric(cut(pv[,step1p], breaks = alphaBinCut )) ]
  rep_helper <- c(table(pv[,grp]))
  pv[ , wt:=rep(alphaBin, rep_helper)]
  
  # return the data.table
  return(pv)
}


create_twostep_weighted_plot_expectation <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "") {
  
  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls
  
  # plot title and file name
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }
  
  # get number of bins and bin sizes based on expectation
  m = 1000000
  nbins = floor(log2(m/sizeBin0 + 1))
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), sizeBin0 * (2^(nbins-1)) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin_step1 = endpointsBin/1000000
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested
  
  # create list of bin specific data.tables
  logp_plot_limit = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]
  
  glist<-list()
  setkey(dat,grp)
  for(i in 1:binsToPlot){
    z <- dat[J(i)][order(Chromosome, Location)]
    z[, log_binalpha := -1*log10(min(z[,wt]))] # -log10 of bin specific alpha
    z[, mapinfo_tmp := seq(Chromosome, Chromosome+1, length.out = .N), by = Chromosome]
    z[, mapinfo := 0.9*( (z[,mapinfo_tmp]-min(z[,mapinfo_tmp])) / (max(z[,mapinfo_tmp])-min(z[,mapinfo_tmp])) ) + 0.05 + i - 1]
    glist[[i]] <- z[order(Chromosome, Location)]
    rm(z)
  }
  
  significant_hits <- dat[step2p <= wt]
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(paste0("~/Dropbox/FIGI/Results/", exposure, "/plots/twostep_wht_", statistic, "_", exposure, "_", paste0(covars, collapse = "_"), filename_suffix, ".png"), height = 720, width = 1280)
  
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  # first bin
  plot(glist[[1]][,mapinfo], glist[[1]][,logstep2p],
       col = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], '#E41A1C','#377EB8'),
       pch = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 19, 20),
       cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,logp_plot_limit),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  lines(glist[[1]][,mapinfo], glist[[1]][,log_binalpha], col = "black", lwd=1)
  
  # remaining bins
  for(i in 2:binsToPlot){
    points(glist[[i]][,mapinfo], glist[[i]][,logstep2p],
           col = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, '#E41A1C', color[i]),
           pch = ifelse(glist[[i]][,SNP] %in% significant_hits$SNP, 19, 20),
           cex = ifelse(glist[[1]][,SNP] %in% significant_hits[,SNP], 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    lines(glist[[i]][,mapinfo], glist[[i]][,log_binalpha],
          col = "black",lwd = 1)
  }
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.7)
  
  dev.off()
  
}



wtf3 <- format_twostep_data_expectation(gxe, 'chiSqG', 5, 0.05)
create_twostep_weighted_plot_expectation(wtf3, exposure, covars, 5, 0.05, 10, 'chiSqG', filename_suffix = "_GECCO_Presentation_expectation_based")



# # meta-analysis ----
# # (show main effects in bigger for presentation)
# 
# 
# gxe <- figi_gwas %>%
#   dplyr::filter(vcfid %in% exposure_subset) %>%
#   dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
#                 sex = ifelse(sex == "Female", 0, 1))
# 
# if(params$is_exposure_categorical == T) {
#   gxe[,params$exposure] <- (as.numeric(gxe[, params$exposure]) - 1) # this is only applicable for categorical variables 
# }
# 
# 
# study_info <- readRDS("~/data/Annotations/FIGI_studygxe_info.rds") %>% 
#   mutate(study_gxe = as.character(study_gxe))
# 
# 
# format_data_metacounts <- function(d, exposure, is_e_categorical, min_cell_size = 0) {
#   tmp <- d
#   drops <- data.frame(table(tmp$outcome, tmp$study_gxe)) %>%
#     filter(Freq <= min_cell_size)
#   tmp <- filter(tmp, !study_gxe %in% unique(drops$Var2)) %>%
#     dplyr::mutate(study_gxe = fct_drop(study_gxe))
#   return(tmp)
# }
# 
# figi_gxe <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_HRC_v2.3_GXE_EUR.rds")
# figi_gxe <- format_data_metacounts(figi_gxe, exposure = 'asp_ref', is_e_categorical = T, 0)
# figi_gxe <- filter(figi_gxe, study_gxe != "ColoCare_1")
# 
# 
# 
# exposure_subset <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", "asp_ref", "_basic_covars_glm.rds"))$vcfid
# figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")
# gxe <- figi_gwas %>%
#   dplyr::filter(vcfid %in% exposure_subset) %>%
#   dplyr::mutate(outcome = ifelse(outcome == "Control", 0, 1),
#                 sex = ifelse(sex == "Female", 0, 1))
# 
#   gxe[,'asp_ref'] <- (as.numeric(gxe[, 'asp_ref']) - 1) # this is only applicable for categorical variables 
# 
# 
# 
# format_d_metaanalysis <- function(d, exposure, is_e_categorical, min_cell_size = 0) {
#   
#   # drop studies with low count cells, keep vars_to_keep
#   if (is_e_categorical == T) {
#     
#     drops <- data.frame(table(d$outcome, d[, exposure], d$study_gxe)) %>%
#       filter(Freq <= min_cell_size)
#     
#     tmp <- filter(d, !study_gxe %in% unique(drops$Var3)) %>%
#       dplyr::mutate(study_gxe = fct_drop(study_gxe))}
#   else {
#     drops <- data.frame(table(d$outcome, d$study_gxe)) %>%
#       filter(Freq <= min_cell_size)
#     tmp <- filter(d, !study_gxe %in% unique(drops$Var2)) %>%
#       dplyr::mutate(study_gxe = fct_drop(study_gxe))}
#   
#   return(tmp)
#   
# }
# 
# gxe_no_low_cell_count <- format_d_metaanalysis(gxe, exposure = 'asp_ref', is_e_categorical = T, 5)
# 
# 
# 
# tmp1 <- filter(study_info, study_gxe %in% figi_gxe$study_gxe)
# tmp2 <- get_counts_outcome_by_group(gxe_no_low_cell_count, 'outcome', 'study_gxe')
# gxe_counts <- full_join(tmp1, tmp2, 'study_gxe')
# 
#   gxe_glm <- get_estimates_e_by_group(gxe_no_low_cell_count, 'outcome', 'asp_ref', 'study_gxe', 'age_ref_imp', 'sex')
#   forest_plot_title <- paste0("All subjects\nOutcome ~ ", 'asp_ref', " + age_ref_imp + sex + study_gxe")
# 
# gxe_meta <- dplyr::full_join(gxe_counts, gxe_glm, 'study_gxe')
# 
# 
# 
# results_meta <- meta::metagen(estimate,
#                               std.error,
#                               data=gxe_meta,
#                               studlab=paste(study_gxe),
#                               comb.fixed = FALSE,
#                               comb.random = TRUE,
#                               method.tau = "SJ",
#                               hakn = TRUE,
#                               prediction=TRUE,
#                               sm="OR",
#                               byvar = study_design)
# 
# 
# png(paste0("~/Dropbox/wtf.png"), height = 17, width = 8.5, units = 'in', res = 150)
# meta::forest(results_meta,
#              layout = "JAMA",
#              # text.predict = "95% CI",
#              # col.predict = "black",
#              leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
#              digits.addcols=0,
#              study.results=F,
#              prediction = F,
#              col.random = 'red')
# grid.text(forest_plot_title, 0.5, .98, gp=gpar(cex=1))
# dev.off()





# 
# 
# # manhattan plot, remove gwas hits ----
# snps_to_exclude <- readRDS("~/data/Annotations/gwas_140_snp_ld_and_region_based_to_exclude_from_gxe.rds")
# 
# gxe <- readRDS("~/data/results/asp_ref/processed/FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan_results.rds")
# gxe_nogwas <- gxe %>%
#   mutate(chr_bp = paste0(Chromosome, ":", Location)) %>%
#   dplyr::filter(!chr_bp %in% snps_to_exclude$SNP)
# 
# 
# 
# wtf <- filter(gxe_chr8, chiSqG > 8)
# wtf2 <- filter(gxe_nogwas, Chromosome == 8) %>% arrange(desc(chiSqG))
# head(wtf2)
# 
# 
# gxe_chr18 <- filter(gxe, Chromosome == 18)
# 
# wtf <- filter(gxe_chr18, chiSqG > 8) %>% arrange(desc(chiSqG))
# head(wtf)
# wtf2 <- filter(gxe_nogwas, Chromosome == 18) %>% arrange(desc(chiSqG))
# head(wtf2)
# 
# 
# 
# gxe_chr1 <- filter(gxe, Chromosome == 1)
# 
# wtf <- filter(gxe_chr1, chiSqG > 30) %>% arrange(desc(chiSqG))
# head(wtf)
# wtf2 <- filter(gxe_nogwas, Chromosome == 1) %>% arrange(desc(chiSqG))
# head(wtf2)
# 
# 
# create_manhattanplot(gxe, 'asp_ref', covars, stat = 'chiSqG', annotation_file = annotation,df = 1, filename_suffix = "_GECCO_Presentation")
# create_manhattanplot(gxe, 'asp_ref', covars, stat = 'chiSqGxE', annotation_file = annotation,df = 1, filename_suffix = "_GECCO_Presentation")
# 
# 
# 
# create_manhattanplot(gxe_nogwas, 'asp_ref', covars, stat = 'chiSqG', annotation_file = annotation,df = 1, filename_suffix = "_GECCO_Presentation_NOGWAS")
# 
# 
# 
# 
# # 2DF results ----
# create_qqplot(gxe_nogwas, exposure, covariates, stat = 'chiSq2df', df = 2, filename_suffix = "_no_gwas")
# create_manhattanplot(gxe_nogwas, 'asp_ref', covars, stat = 'chiSq2df', annotation_file = annotation_file,df = 2, filename_suffix = "_no_gwas")
# 
# # 3DF results ----
# create_qqplot(gxe_nogwas, exposure, covariates, stat = 'chiSq3df', df = 3, filename_suffix = "_no_gwas")
# create_manhattanplot(gxe_nogwas, exposure, covariates, stat = 'chiSq3df', annotation_file = annotation_file,df = 3, filename_suffix = "_no_gwas")
# 
