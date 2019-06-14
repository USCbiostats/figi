#========================================================#
# Manhattan and QQ Plots
# (qqman)
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)

rm(list = ls())

# original results
results <- do.call(rbind, lapply(list.files(path = './results/', full.names = T, pattern = "results_n71783_aspirin_sex_age_pcs_studyname_chr"), fread, stringsAsFactors = F))
names(results)

# calculate lambdas ------
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
getlambda1000 <- function(lambda) {
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  lambda1000
}


#--------------------------------------------------------#
# Marginal G Results
#--------------------------------------------------------#
results_G <- results %>%
  #filter(chiSqG > 3.841) %>% 
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaG, P)

cases <- unique(results_G$Cases)
controls <- unique(results_G$Subjects) - unique(results_G$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_G$P)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_aspirin_chiSqG_qq.png", height = 600, width = 800)
qqman::qq(results_G$P, 
          xlab = "Expected chi-squared value", 
          ylab="Observed chi-squared value",
          main = "G Main Effects Results\noutcome ~ Gdosage+age_ref_imp+sex+studyname+PC1-10+aspirin",
          sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4)))) # ~~ adds space lol
dev.off()

results_man <- filter(results_G, P < 0.05)
png("./figures/GxEScanR_aspirin_chiSqG_manhattan.png", width = 800, height = 600)
qqman::manhattan(results_man, main = "G Main Effects Results\noutcome ~ Gdosage+age_ref_imp+sex+studyname+PC1-10+aspirin")
dev.off()


# same manhattan graph but using EasyStrata
# (annotation file with 95 GWAS hits created in another script..)

# two sets of annotations - previously reported in literature vs newly reported in Huyghe 2019
jh_lit <- fread("~/bin/EasyStrata/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no",
         PREVIOUSLY_REPORTED == "yes") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "green", 
         SNP = paste(Chr, Pos, sep = ":"))
write.table(jh_lit, file = "~/bin/EasyStrata/JH_Literature.tsv", quote = F, row.names = F, sep = '\t')

jh_new <- fread("~/bin/EasyStrata/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE == "no",
         PREVIOUSLY_REPORTED == "no") %>% 
  dplyr::rename(Chr = CHROM,
                Pos = POSITION) %>% 
  mutate(Colour = "red", 
         SNP = paste(Chr, Pos, sep = ":"))
write.table(jh_new, file = "~/bin/EasyStrata/JH_New.tsv", quote = F, row.names = F, sep = '\t')

jh <- rbind(jh_lit, jh_new)
write.table(jh, file = "~/bin/EasyStrata/JH_Annotation.tsv", quote = F, row.names = F, sep = '\t')

# want to make sure i have all the significant hits in the results object... but can't. filtered out? 
# so nevermind. 
jh_lit_join <- inner_join(jh_lit, results_G, by = 'SNP') # 53 - missing 2
jh_new_join <- inner_join(jh_new, results_G, by = 'SNP') # 37 - missing 3

results_G_EasyStrata <- results_G %>%
  mutate(A1 = Reference, A2 = Alternate) %>% 
  filter(P < 0.05) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, betaG, P)

jh_lit_join <- inner_join(jh_lit, results_G_EasyStrata, by = 'SNP') # 51 - missing 4 ...
jh_new_join <- inner_join(jh_new, results_G_EasyStrata, by = 'SNP') # 37 - missing 3 ...

write.table(results_G_EasyStrata, file = "~/bin/EasyStrata/results_G_asp_ref_sex_age_pcs_studyname_N72438.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/results_asp_ref_EasyStrata_G.ecf")





#------------------------------------------------------------#
# output list of top hits to extract from VCF and BDose files
results_G_filter <- filter(results_G, P <= 5E-8)
vcfout <- results_G_filter %>% 
  dplyr::select(CHR, BP)
write.table(vcfout, file = "GxEScanR_asp_ref_top1262_vcf_CHR_BP.txt", quote = F, row.names = F, col.names = F, sep = '\t')

for(chr in 1:22) {
  snpsToGet <- results_G_filter[which(results_G_filter$CHR == chr), 'SNP']
  saveRDS(snpsToGet, file = paste0("GxEScanR_asp_ref_top1262_bdose_SNP_chr", chr, ".rds"))
}
#------------------------------------------------------------#


# table for presentation of top hits.. 
jh_results_G <- inner_join(jh, results_G_filter, by = 'SNP') %>% 
  mutate(Author = AUTHOR_FIRST_REPORTED, 
         Year = YEAR_FIRST_REPORTED,
         Beta = BETA, 
         P = P.x,
         Ref_Alt = paste(OTHER_ALLELE, RISK_ALLELE, sep = "/"),
         Ref_Alt_new = paste(Reference, Alternate, sep = "/"),
         Beta_new = betaG, 
         P_new = P.y) %>% 
  dplyr::select(LOCUS, SNP, PREVIOUSLY_REPORTED, RSID, Author, Year, Ref_Alt, Ref_Alt_new, Beta, P, Beta_new, P_new)
saveRDS(jh_results_G, file = "asp_ref_G_tophits_table.rds")



#--------------------------------------------------------#
# GxE results
#--------------------------------------------------------#
results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P)

cases <- unique(results_GxE$Cases)
controls <- unique(results_GxE$Subjects) - unique(results_GxE$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_GxE$P)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_aspirin_chiSqGxE_qq.png", height = 600, width = 800)
qqman::qq(results_GxE$P, 
          xlab = "Expected chi-squared value", 
          ylab = "Observed chi-squared value",
          main = "GxE Results\noutcome ~ Gdosage*aspirin+age_ref_imp+sex+studyname+PC1-10",
          sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4)))) # ~~ adds space lol
dev.off()

results_man <- filter(results_GxE, P < 0.05)
png("./figures/GxEScanR_aspirin_chiSqGxE_manhattan.png", width = 800, height = 600)
qqman::manhattan(results_man, main = "GxE Results\noutcome ~ Gdosage*aspirin+age_ref_imp+sex+studyname+PC1-10")
dev.off()



results_GxE_EasyStrata <- results_GxE %>%
  mutate(A1 = Reference, A2 = Alternate) %>% 
  filter(P < 0.05) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, betaGxE, P)

write.table(results_GxE_EasyStrata, file = "~/bin/EasyStrata/results_GxE_asp_ref_sex_age_pcs_studyname_N72438.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/results_asp_ref_EasyStrata_GxE.ecf")

z <- inner_join(annot, results_GxE, by = 'SNP') # 90
zz <- inner_join(annot, results_GxE_EasyStrata, by = 'SNP') # 3

results_GxE_filter <- filter(results_GxE, P <= 5E-8)
vcfout <- results_G_filter %>% 
  dplyr::select(CHR, BP)



# lets see
results_GxE_filter <- filter(results_GxE, P <= 5E-8) %>% 
  mutate(v1 = "chromosome",
         v2 = 1) %>% 
  dplyr::select(v1, CHR, BP, Reference, Alternate, v2)
write.table(results_GxE_filter, "~/Dropbox/tmp_getannot.txt", quote= F, row.names = F, col.names = F)


# annotation results from: https://snp-nexus.org/test/snpnexus_18324/results.html
snp_nexus <- fread("~/Dropbox/ncsnp_18324.txt") %>% 
  mutate(Chromosome = gsub("chr", "", Chromosome))


#--------------------------------------------------------#
# G+GxE results
#--------------------------------------------------------#
# Marginal G Results
results_GGxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chi2df, df = 2, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, P)

cases <- unique(results_GGxE$Cases)
controls <- unique(results_GGxE$Subjects) - unique(results_GGxE$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda2df(results_GGxE$P)
lambda1000 <- getlambda1000(lambda)

png("GxEScanR_asp_ref_chiSq2df_qq.png", height = 600, width = 800)
qqman::qq(results_GGxE$P, 
          xlab = "Expected chi-squared value (2df)", 
          ylab = "Observed chi-squared value (2df)",
          main = "G/GxE 2DF Results\noutcome ~ GDosage+GDosage*asp_ref+age_ref_imp+sex+studyname+PC1-10",
          sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4)))) # ~~ adds space lol
dev.off()

results_man <- filter(results_GGxE, P < 0.05)
png("GxEScanR_asp_ref_chiSq2df_manhattan.png", width = 800, height = 600)
qqman::manhattan(results_man, main = "G/GxE 2DF Results\noutcome ~ GDosage+GDosage*asp_ref+age_ref_imp+sex+studyname+PC1-10")
dev.off()

results_GGxE_EasyStrata <- results_GGxE %>%
  mutate(A1 = Reference, A2 = Alternate) %>% 
  filter(P < 0.05) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, P)
write.table(results_GGxE_EasyStrata, file = "~/bin/EasyStrata/results_2DF_asp_ref_sex_age_pcs_studyname_N72438.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("~/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/results_asp_ref_EasyStrata_2DF.ecf")

# which one of these is like... the one from GxE
ztmp <- results_GxE_filter %>% 
  mutate(SNP = paste(CHR, BP, sep = ":"), 
         Chr = CHR, Pos = BP, Colour = 'orchid')

jh_edited <- bind_rows(ztmp, jh)
write.table(jh_edited, file = "~/bin/EasyStrata/JH_Annotation_edited.tsv", quote = F, row.names = F, sep = '\t')

z <- inner_join(results_GGxE_EasyStrata, ztmp, by = 'SNP')







