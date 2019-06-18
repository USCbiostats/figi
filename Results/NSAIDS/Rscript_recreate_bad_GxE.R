#=======================================================================================#
# just recreating shitty GxE
#=======================================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(Cairo)
rm(list = ls())

# new annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_72145"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

# ge_only results (for 2-step methods, case only, control only)
colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
ge <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref/ge_only/", full.names = T, pattern = "results_GEOnly_asp_ref_sex_age_pc10_studygxe_72145"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID)) %>% dplyr::select(-drop)

# create merged file, only use this one for starting point for all over analyses
names(gxe)
names(ge)
results <- inner_join(gxe, ge[, c("betaGE", "chiSqGE", "betaCase", "ChiSqCase", "betaControl", "chiSqControl", "ID")], by = "ID")

rm(ge);rm(gxe);rm(ukb_filter)


# Functions ------
getlambda <- function(pvals) {
  chisq <- qchisq(1-pvals, 1)
  lambda <- round(median(chisq)/qchisq(0.5,1),4)
  lambda
}
getlambda2df <- function(pvals) {
  chisq <- qchisq(1-pvals, 2)
  lambda <- round(median(chisq)/qchisq(0.5,2),4)
  lambda
}
getlambda3df <- function(pvals) {
  chisq <- qchisq(1-pvals, 3)
  lambda <- round(median(chisq)/qchisq(0.5,3),4)
  lambda
}
getlambda1000 <- function(lambda, cases, controls, cases1000, controls1000) {
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  lambda1000
}

wrapper_qq <- function(df, p, title, filename) {
  cases <- unique(df$Cases)
  controls <- unique(df$Subjects) - unique(df$Cases)
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda <- getlambda(df[,p])
  lambda1000 <- getlambda1000(lambda, cases, controls, cases1000, controls1000)
  
  # QQ Plot
  CairoPNG(height = 720, width = 1280, file = paste0("figures/", filename))
  qqman::qq(df[,p], 
            xlab = "Expected -log10(p)", 
            ylab = "Observed -log10(p)",
            main = title,
            cex.main = 1.8, 
            cex.axis = 1.5, 
            cex.lab = 1.5,
            cex.sub = 1.5,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.5)
  dev.off()
}

wrapper2df_qq <- function(df, p, title, filename) {
  cases <- unique(df$Cases)
  controls <- unique(df$Subjects) - unique(df$Cases)
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda <- getlambda2df(df[,p])
  lambda1000 <- getlambda1000(lambda, cases, controls, cases1000, controls1000)
  
  # QQ Plot
  CairoPNG(height = 720, width = 1280, file = paste0("figures/", filename))
  qqman::qq(df[,p], 
            xlab = "Expected -log10(p)", 
            ylab = "Observed -log10(p)",
            main = title,
            cex.main = 1.8, 
            cex.axis = 1.5, 
            cex.lab = 1.5,
            cex.sub = 1.5,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.5)
  dev.off()
}

#--------------------------------------------------------#
# GxE results ----
#--------------------------------------------------------#
results_GxE <- results %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

cases <- unique(results_GxE$Cases)
controls <- unique(results_GxE$Subjects) - unique(results_GxE$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_GxE$P)
lambda1000 <- getlambda1000(lambda, cases, controls, cases1000, controls1000)

CairoPNG(height = 720, width = 1280, file = "figures/GxEScanR_asp_ref_chiSqGxE_qq_BAD.png")
qqman::qq(results_GxE$P, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "GxNSAIDs QQ Plot\noutcome ~ dosage*asp_ref + age_ref_imp + sex + study_gxe + PC1-10",
          cex.main = 1.8, 
          cex.axis = 1.5, 
          cex.lab = 1.5,
          cex.sub = 1.5,
          col = 'blue4')
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.5)
dev.off()


# Manhattan Plot
results_GxE_EasyStrata <- results_GxE %>% 
  mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
  filter(!(P > 0.05 & annot == 0)) %>% 
  dplyr::rename(CHR = Chromosome,
                BP = Location, 
                A1 = Reference, 
                A2 = Alternate) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, betaGxE, P)

write.table(results_GxE_EasyStrata, file = "~/data/Annotations/EasyStrata_GxE_asp_ref_sex_age_pcs_studygxe_N72145_BAD.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("results_asp_ref_EasyStrata_GxE_LDAnnot_BAD.ecf")


# Further investigate the 10 'fake' significant results
results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

results_GxE_EasyStrata <- results_GxE %>%
  filter(P < 0.05) %>%
  dplyr::rename(CHR = Chromosome,
                BP = Location,
                A1 = Reference,
                A2 = Alternate) %>%
  dplyr::select(ID, CHR, BP, A1, A2, betaGxE, P)

write.table(results_GxE_EasyStrata, file = "~/data/Annotations/EasyStrata_GxE_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_GxE_LDAnnot.ecf")


# extract 10 SNPs from VCF files
results_GxE_filter <- filter(results_GxE, P <= 5E-8)
vcfout <- results_G_filter %>% 
  dplyr::select(CHR, BP) # i checked to make sure the dosage issue is present in original files (they are)


# Extract 10 SNPs from BDose files..
chr <- 22
for(chr in c(2,3,6,7,10,16,17,18,21,22)) {
  snpsToGet <- results_GxE_filter[which(results_GxE_filter$Chromosome == chr), 'SNP']
  saveRDS(snpsToGet, file = paste0("./qc/GxEScanR_GxE_asp_ref_faketop10_extract_from_bdose_chr", chr, ".rds"))
}
