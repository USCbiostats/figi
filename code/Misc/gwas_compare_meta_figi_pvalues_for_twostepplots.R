#=============================================================================#
# 01/07/2020
#
# Compare fredHutch gwas meta-analysis results with gxescan gwas results
# Decide which results to use in gxe two-step plots
#
#=============================================================================#
library(tidyverse)
library(data.table)
library(EasyStrata)


meta <- fread("~/gecco_125k_gwas_eur/MarginalMeta_HRC_EUR_only_Results.tsv.gz") %>% 
  mutate(MarkerName = `#MarkerName`) %>% 
  dplyr::select(MarkerName, Effect, P.value)

figi <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results.rds") %>% 
  mutate(MarkerName = paste0(SNP, "_", Reference, "/", Alternate)) %>% 
  dplyr::select(MarkerName, Chromosome, Location, betaG, chiSqG)

out <- inner_join(meta, figi, 'MarkerName') %>% 
  mutate(pval_meta = P.value, 
         pval_figi = pchisq(chiSqG, df = 1, lower.tail = F),
         # logp_meta = -log10(P.value), 
         # logp_figi = -log10(pval_figi), 
         CHR = Chromosome, 
         BP = Location) %>% 
  dplyr::select(MarkerName, CHR, BP, pval_meta, pval_figi)

# generate miami plots as a start? this way you don't need to filter by p value, just plot all 


wtf <- filter(out, CHR == 1)
names(wtf)
write.table(wtf, file = "/media/junk/test_chr1.txt", quote = F, row.names = F, col.names = T, sep = '\t')

EasyStrata("~/Dropbox/test.ecf")

plot(wtf$logp_meta, wtf$logp_figi)




# ------ Mass produce these things ----- #
# use this function to filter data by chromosome and output miami plot
# you need the 'out' object (inner_join of meta + figi.. )


create_manhattanplot_gwas_check <- function(dat, chr) {
  
  # format data for EasyStrata package, output to file
  tmp <- dat %>% 
    filter(CHR == chr)
  write.table(tmp, file = paste0("/media/junk/EasyStrata_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0("/home/rak/Dropbox/")
  ecf2 <- "MarkerName;CHR;BP;pval_meta;pval_figi"
  ecf3 <- "character;numeric;numeric;numeric;numeric"
  ecf4 <- paste0("/media/junk/EasyStrata_chr", chr, ".txt")
  ecf_file_name <- paste0("/media/junk/EasyStrata_chr", chr, ".ecf")
  
  cat(paste0("DEFINE	--pathOut ", ecf1, "
      --acolIn ", ecf2, "
      --acolInClasses ", ecf3, "
      --strMissing NA
      --strSeparator TAB

      EASYIN	--fileIn ", ecf4, "

      START EASYX

      ################
      ## Miami Plot 
      ################
      
      MIAMIPLOT
      --colInChr CHR
      --colInPos BP
      --colMIAMIPlotUp pval_meta
      --colMIAMIPlotDown pval_figi
      --numPvalOffset 0.05
      --numWidth 1280
      --numHeight 720	
      --blnYAxisBreak 1
      --numYAxisBreak 22
      ## Annotations
      --fileAnnot /home/rak/data/Annotations/gwas_140_ld_annotation.txt
      # i'm defining locus myself, so leave numAnnotPosLim = 1
      --numAnnotPosLim 1
      ## Colors + axis breaks
      #--astrDefaultColourChrUp dodgerblue3;dodgerblue4
      #--astrDefaultColourChrDown dodgerblue3;dodgerblue4
      ## horizontal line
      --anumAddPvalLine 5e-8
      --astrAddPvalLineCol coral3
      --anumAddPvalLineLty 6
      ## others
      --numDefaultSymbol 20
      --numDefaultCex 0.4
      --numCexAxis 1.5
      --numCexLab 2
      --arcdSymbolCritUp pval_meta <5e-8
      --arcdSymbolCritDown pval_figi <5e-8
      --anumSymbolUp 18
      --anumSymbolDown 18
      --arcdColourCritUp pval_meta <5e-8
      --arcdColourCritDown pval_figi <5e-8
      --astrColourUp black
      --astrColourDown black
      --arcdCexCritUp pval_meta <5e-8
      --arcdCexCritDown pval_figi <5e-8
      --anumCexUp 1.2
      --anumCexDown 1.2
      
      STOP EASYX"), file = ecf_file_name, append = F)
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}

create_manhattanplot_gwas_check(out, 1)

for(x in 3:22){
  create_manhattanplot_gwas_check(out, x)
}



# ------ One more thing - are there any significant markers in the meta-analysis that don't overlap with figi gwas results?

library(tidyverse)
library(data.table)
library(EasyStrata)
rm(list = ls())


meta <- fread("~/gecco_125k_gwas_eur/MarginalMeta_HRC_EUR_only_Results.tsv.gz") %>% 
  mutate(MarkerName = `#MarkerName`) %>% 
  dplyr::select(MarkerName, Chr, Position, Effect, P.value)

figi <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results.rds") %>% 
  mutate(MarkerName = paste0(SNP, "_", Reference, "/", Alternate)) %>% 
  dplyr::select(MarkerName, betaG, chiSqG)

out <- anti_join(meta, figi, 'MarkerName') %>% 
  mutate(pval_meta = P.value,
         CHR = Chr, 
         BP = Position) %>% 
  dplyr::select(MarkerName, CHR, BP, pval_meta)

write.table(out, file = "/media/junk/test_antijoin.txt", quote = F, row.names = F, col.names = T, sep = '\t')

EasyStrata("~/Dropbox/test2.ecf")


# take the markers, extract from results and look at the mafs...

wtf <- filter(out, pval_meta < 5e-8) %>% 
  mutate(chromStart = BP - 1) %>% 
  dplyr::select(CHR, chromStart, BP)
write.table(wtf, file = "~/gecco_125k_gwas_eur/meta_analysis_only_sig_snps.bed", quote = F, row.names = F, col.names = F, sep = '\t')  


# after tabix ..
xx <- fread("~/gecco_125k_gwas_eur/meta_analysis_only_sig_snps.out")
















# ---------------- make log/log plot of 'suggestive hits' (in both analyses) -------------------- #

meta <- fread("~/gecco_125k_gwas_eur/MarginalMeta_HRC_EUR_only_Results.tsv.gz") %>% 
  mutate(MarkerName = `#MarkerName`) %>% 
  dplyr::select(MarkerName, Effect, P.value)


figi <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results.rds") %>% 
  mutate(MarkerName = paste0(SNP, "_", Reference, "/", Alternate),
         pval_figi = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
  dplyr::select(MarkerName, Chromosome, Location, betaG, chiSqG, pval_figi) %>% 
  dplyr::filter(pval_figi < 5e-6)


out <- inner_join(meta, figi, 'MarkerName') %>% 
  mutate(pval_meta = P.value,
         logp_meta = -log10(pval_meta),
         logp_figi = -log10(pval_figi),
         CHR = Chromosome, 
         BP = Location) %>% 
  dplyr::select(MarkerName, CHR, BP, pval_meta, pval_figi, logp_meta, logp_figi)

plot(out$logp_meta, out$logp_figi)


# the other way around

meta <- fread("~/gecco_125k_gwas_eur/MarginalMeta_HRC_EUR_only_Results.tsv.gz") %>% 
  mutate(MarkerName = `#MarkerName`) %>% 
  dplyr::select(MarkerName, Effect, P.value) %>% 
  filter(P.value < 5e-6)


figi <- readRDS("~/data/results/gwas/processed/FIGI_v2.3_gwasset_basic_covars_gxescan_results.rds") %>% 
  mutate(MarkerName = paste0(SNP, "_", Reference, "/", Alternate),
         pval_figi = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
  dplyr::select(MarkerName, Chromosome, Location, betaG, chiSqG, pval_figi)


out <- inner_join(meta, figi, 'MarkerName') %>% 
  mutate(pval_meta = P.value,
         logp_meta = -log10(pval_meta),
         logp_figi = -log10(pval_figi),
         CHR = Chromosome, 
         BP = Location) %>% 
  dplyr::select(MarkerName, CHR, BP, pval_meta, pval_figi, logp_meta, logp_figi)

plot(out$logp_meta, out$logp_figi, main = "logp scatter plot - pval_meta < 5e-6")


