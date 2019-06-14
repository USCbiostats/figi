#=======================================================================================#
# NSAIDS - asp_ref Results
# NO UKBIOBANK!!!!!!!!!!!!!!!!!!!!!!!
# 02/28/2019
# 
# Generate Plots
# Implement 2-step methods
#
# Notes:
# - filtered 'bad' SNPs
#   dosage~UKB+PC @ < 5E-8 (i think we should stick to this for the presentation, instead of comparing this and that)
#=======================================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
rm(list = ls())

# remove markers based on linear model dosage ~ ukb indicator var + PCs
# if you're annotating marginal G hits, there are some GWAS hits that get filtered out through this step. 
# Add those markers back (don't filter from results)

# new annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
  mutate(SNP = paste(Chr, Pos, sep = ":"))

# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_noUKB"), fread, stringsAsFactors = F)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID))

# ge_only results (for 2-step methods, case only, control only)
colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
ge <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/ge_only/", full.names = T, pattern = "results_GEOnly_asp_ref_sex_age_pc10_studygxe_noUKB"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
  mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
  filter(!duplicated(ID)) %>% dplyr::select(-drop)


# you want to make sure NOT to drop annotations from manhattan plots. let's make sure there's good overlap...
# 90 OF THE MARKERS ARE IN MY RESULTS, MISSING 5, 2 of them were previously reported
annotation_original <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_061918.tsv")[, -c('LL_95PCT_CI_BETA', 'UL_95PCT_CI_BETA')] %>% 
  filter(SUGGESTIVE=="no") %>% 
  separate(VARIANT, into = c("SNP", "alleles"), sep = "_")

check <- inner_join(gxe, annotation_original, by = "SNP")
check <- anti_join(annotation_original, gxe, by = "SNP")

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
getlambda3df <- function(pvals) {
  chisq <- qchisq(1-pvals, 3)
  lambda <- round(median(chisq)/qchisq(0.5,3),4)
  lambda
}
getlambda1000 <- function(lambda) {
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  lambda1000
}

#-----------------------------------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------------------------------#
results_G <- gxe %>%
  mutate(P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P)

cases <- unique(results_G$Cases)
controls <- unique(results_G$Subjects) - unique(results_G$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_G$P)
lambda1000 <- getlambda1000(lambda)

# QQ Plot
png("./figures/GxEScanR_asp_ref_noUKB_chiSqG_qq.png", height = 720, width = 1280)
qqman::qq(results_G$P, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "G Main Effects Results (N = 57440)\noutcome ~ dosage+age_ref_imp+sex+studyname+PC1-10+asp_ref",
          #sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3)
dev.off()


# Manhattan Plot
results_G_EasyStrata <- results_G %>% 
  mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
  filter(!(P > 0.05 & annot == 0)) %>% 
  dplyr::rename(CHR = Chromosome,
                BP = Location, 
                A1 = Reference, 
                A2 = Alternate) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, betaG, P)

check <- filter(results_G_EasyStrata, P > 0.05)

write.table(results_G_EasyStrata, file = "~/data/Annotations/EasyStrata_MarginalG_asp_ref_noUKB_sex_age_pcs_studygxe_N57440.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("noUKB_results_asp_ref_EasyStrata_G_LDAnnot.ecf")


# Manhattan Plot
# for every chromosome? 
# for(chr in c(1:22)) {
#   results_G_EasyStrata_ukb_filter <- results_G %>% 
#     filter(P < 0.05,
#            Chromosome == chr ) %>% 
#     dplyr::rename(CHR = Chromosome,
#                   BP = Location, 
#                   A1 = Reference, 
#                   A2 = Alternate) %>% 
#     dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
#   
#   write.table(results_G_EasyStrata_ukb_filter, file = paste0("~/data/Annotations/EasyStrata_MarginalG_asp_ref_sex_age_pcs_studygxe_N72145_ukb_filter_chr", chr, ".txt"), quote = F, row.names = F, sep = '\t')
# }
# 
# EasyStrata("results_asp_ref_EasyStrata_G_LDAnnot_ukb_filter.ecf") # have to change file manually lol.. for now. 





#------------------------------------------------------------#
# output list of top hits to extract from VCF and BDose files
# results_G_filter <- filter(results_G, P <= 5E-8)
# vcfout <- results_G_filter %>% 
#   dplyr::select(CHR, BP)
# write.table(vcfout, file = "GxEScanR_asp_ref_top1262_vcf_CHR_BP.txt", quote = F, row.names = F, col.names = F, sep = '\t')
# 
# for(chr in 1:22) {
#   snpsToGet <- results_G_filter[which(results_G_filter$CHR == chr), 'SNP']
#   saveRDS(snpsToGet, file = paste0("GxEScanR_asp_ref_top1262_bdose_SNP_chr", chr, ".rds"))
# }
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
# GxE results ----
#--------------------------------------------------------#
results_GxE <- gxe %>%
  mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

cases <- unique(results_GxE$Cases)
controls <- unique(results_GxE$Subjects) - unique(results_GxE$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_GxE$P)
lambda1000 <- getlambda1000(lambda)


png("./figures/GxEScanR_asp_ref_noUKB_chiSqGxE_qq.png", height = 720, width = 1280)
qqman::qq(results_GxE$P, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "GxNSAIDs Results (N = 57440)\noutcome ~ Gdosage*asp_ref+age_ref_imp+sex+studyname+PC1-10",
          #sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3)
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

write.table(results_GxE_EasyStrata, file = "~/data/Annotations/EasyStrata_GxE_asp_ref_noUKB_sex_age_pcs_studygxe_N57440.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("noUKB_results_asp_ref_EasyStrata_GxE_LDAnnot.ecf")


#--------------------------------------------------------#
# 2DF (G+GxE) results ----
#--------------------------------------------------------#
results_2df <- gxe %>%
  mutate(P = pchisq(chi2df, df = 2, lower.tail = F)) %>%
  dplyr::select(ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P) %>% 
  filter(!ID %in% ukb_filter$ID)

cases <- unique(results_2df$Cases)
controls <- unique(results_2df$Subjects) - unique(results_2df$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda2df(results_2df$P)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_asp_ref_chiSq2df_qq.png", height = 720, width = 1280)
qqman::qq(results_2df$P, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "G_GxNSAIDs 2DF Results\noutcome ~ GDosage+GDosage*asp_ref+age_ref_imp+sex+studyname+PC1-10",
          #sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3)
dev.off()


results_2df_EasyStrata <- results_2df %>%
  filter(P < 0.05) %>%
  dplyr::rename(CHR = Chromosome,
                BP = Location,
                A1 = Reference,
                A2 = Alternate) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, P)

write.table(results_2df_EasyStrata, file = "~/data/Annotations/EasyStrata_2DF_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("~/Dropbox/FIGI/Results/NSAIDS/results_asp_ref_EasyStrata_2DF.ecf")

# which one of these hits in the 2df test are recapitulating the ones from GxE (at this point should be... NONE)
# ztmp <- results_GxE_filter %>% 
#   mutate(SNP = paste(CHR, BP, sep = ":"), 
#          Chr = CHR, Pos = BP, Colour = 'orchid')
# 
# jh_edited <- bind_rows(ztmp, jh)
# write.table(jh_edited, file = "~/bin/EasyStrata/JH_Annotation_edited.tsv", quote = F, row.names = F, sep = '\t')
# 
# z <- inner_join(results_GGxE_EasyStrata, ztmp, by = 'SNP')


#--------------------------------------------------------#
# Miami Plot comparing G and GxE? ----
#--------------------------------------------------------#
names(results_G)
names(results_2df)

results_G_2df_miami <- inner_join(results_G, results_2df[, c("ID", "P")], by = "ID") %>% 
  rename(pval_G = P.x, 
         pval_2df = P.y,
         CHR = Chromosome, 
         BP = Location, 
         A1 = Reference, 
         A2 = Alternate) %>% 
  filter(pval_G < 0.05 | pval_2df < 0.05) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, pval_G, pval_2df)

names(results_G_2df_miami)

write.table(results_G_2df_miami, file = "~/data/Annotations/EasyStrata_G_2DF_MIAMI_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')


EasyStrata("~/Dropbox/FIGI/Results/NSAIDS/results_asp_ref_EasyStrata_Miami_G_2DF.ecf")


#--------------------------------------------------------#
# GE, Case, Control QQ and/or Manhattan plots ----
#--------------------------------------------------------#
ge_qq <- inner_join(gxe[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE",  "chiSqG", "chiSqGxE")], ge[, c("SNP", "chiSqGE", "betaCase", "ChiSqCase", "betaControl", "chiSqControl")], by = "SNP") %>% 
  mutate(pval.GE = pchisq(chiSqGE, df = 1, lower.tail = F),
         pval.Case = pchisq(ChiSqCase, df = 1, lower.tail = F),
         pval.Control = pchisq(chiSqControl, df = 1, lower.tail = F)) %>% 
  arrange(pval.GE)



# GE statistic qq plot
cases <- unique(ge_qq$Cases)
controls <- unique(ge_qq$Subjects) - unique(ge_qq$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000
lambda <- getlambda(ge_qq$pval.GE)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_asp_ref_chiSqGE_qq.png", height = 720, width = 1280)
qqman::qq(ge_qq$pval.GE, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "GxNSAIDs G|E Results\nGDosage ~ asp_ref+age_ref_imp+sex+studyname+PC1-10",
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3)
dev.off()




# Case-Only statistic qq plot
cases <- unique(ge_qq$Cases)
controls <- unique(ge_qq$Subjects) - unique(ge_qq$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000
lambda <- getlambda(ge_qq$pval.Case)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_asp_ref_ChiSqCase_qq.png", height = 720, width = 1280)
qqman::qq(ge_qq$pval.Case, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "GxNSAIDs Case-Only Results\nGDosage ~ asp_ref+age_ref_imp+sex+studyname+PC1-10",
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4))), cex.sub = 1.3)
dev.off()

ge_qq_Case_EasyStrata <- ge_qq %>%
  rename(P = pval.Case) %>% 
  filter(P < 0.05) %>%
  dplyr::rename(CHR = Chromosome,
                BP = Location,
                A1 = Reference,
                A2 = Alternate) %>% 
  dplyr::select(SNP, CHR, BP, A1, A2, P)
write.table(ge_qq_Case_EasyStrata, file = "~/data/Annotations/EasyStrata_CaseOnly_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("~/Dropbox/FIGI/Results/NSAIDS/results_asp_ref_EasyStrata_G_CASEONLY.ecf")






# Control-Only statistic qq plot
png("./figures/GxEScanR_asp_ref_chiSqControl_qq.png", height = 720, width = 1280)
qqman::qq(ge_qq$pval.Control, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "GxNSAIDs Control-Only Results\nGDosage ~ asp_ref+age_ref_imp+sex+studyname+PC1-10",
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4))), cex.sub = 1.3)
dev.off()



#--------------------------------------------------------#
# D|G 2-step Kooperberg ----
#--------------------------------------------------------#
# no need for joins (filtering on D|G)
koop <- gxe %>% 
  mutate(pval_G = pchisq(chiSqG, df = 1, lower.tail = F),
         pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F), 
         logpval_GxE = -log10(pval_GxE)) %>% 
  arrange(pval_G) %>% 
  mutate(rank = seq(1, nrow(.)))

koop_filter <- koop %>% 
  dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval_G, logpval_GxE)

# koop_filter_sortGxE <- arrange(koop_filter, logpval_GxE)
# head(koop_filter)


#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 

m = 5485386
sizeBin0 = 5
alpha = 0.05
# alpha1 = 0.001 # first bin alpha level.. 
# N=72438; Nbio=100; ncp_GxE=0

nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction 
nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes 
endpointsBin = cumsum(sizeBin) # endpoints of the bins 
alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested 

# not sure what this is for...
# subs = cumsum(sizeBin) <= alpha1*m 
# alphaBin_subset = alpha / sizeBin * subs/ sum(subs==TRUE) 
# alphaBin_linear = alpha / sizeBin / nbins 
pv <- data.table(koop_filter)
head(pv)

## Function 11: Create Bin and significance level for weighted testing approach (needs to be sorted by p-value)
num = sizeBin[1] # initial bin size

rk.pv <- c(1:nrow(pv))
grp <- ceiling(log(rk.pv/num+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
pv[,Bin:=grp] # assigning group to the p values.. 
setkey(pv,Bin)
for(i in 1:max(grp)) 
{
  pv[J(i),Threshold:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
}

# see which are significant...
results <- pv %>% 
  mutate(logThreshold = -log10(Threshold))

results_sig <- results %>% 
  filter(logpval_GxE >= logThreshold)

summary(results$logpval_GxE)
results_suggestive <- filter(results, logpval_GxE >= 5)
table(results_suggestive$Bin)





# for PLOTS
num = nbins
num = 15
min.p = 12 # log scale...
last.sig = -log10(1.394709e-13)

# mapinfo is for scaling x axis plot
setnames(pv, c(8,9,10,3), c('p', 'grp','wt','MapInfo'))
pv[, y:= pv[,p]]
head(pv)

# create list where each component contains a bin
# log transfor bin alpha value
# create 'x', which is ordered by SNP BP sort of
glist<-list()
for(i in 1:num){
  t <- pv[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}
glist[[1]]
glist[[2]]

# trying to understand code above (just arrange BP by bins looks like)
# (Scale mapinfo for each Bin to range between 0.1-0.9 for neatness, and add a unit increase for successive Bin)
# x <- pv[1:5, MapInfo]
# normalized = (x-min(x))/(max(x)-min(x))
# normalized_scaled = 0.8 * normalized + 0.1
# x;normalized;normalized_scaled

# CREATE PLOT
head(glist[[1]]) # for reference

png("./figures/GxEScanR_asp_ref_KooperbergDG_Weighted.png", width = 1280, height = 720)
color <- rep(c("blue","olivedrab4"),100)
plot(glist[[1]][,x], glist[[1]][,y], 
     col = "blue", 
     xlab="Bin # for step1 p-values", 
     ylab="-log10(step2 p-values)", 
     xlim=c(0,num), 
     ylim=c(0,min.p), 
     axes=F, pch=19, cex=0.5)
lines(glist[[1]][,x], glist[[1]][,ref],
      col = "black",lwd=2)

# the rest of the points ...  =|
# (adding to current plot..)
for(i in 2:num){
  points(glist[[i]][,x], glist[[i]][,y], 
         col = color[i], pch = 19, cex = 0.5)      
  lines(glist[[i]][,x], glist[[i]][,ref],
        col = "black",lwd = 2)
}

## the last bin..
## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
## let's only plot the first 15 bins for now, so change code a bit above. 
# points(glist[[num]][,x], glist[[num]][,y], 
#        col= color[num], pch = 19, cex = 0.5)
# lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 0.8)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=0.8)
title(main = "GxNSAIDs D|G Weighted Hypothesis Testing", sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.5, cex.sub = 1.2)

dev.off()




#--------------------------------------------------------#
# G|E 2-step Murcray ----
#--------------------------------------------------------#
# start with weighted hypothesis testing
# need to obtain ranks from the step 1 statistic. Then test hypothesis using bins..
# make sure GxE p value is in the log scale for the plots
murcray <- inner_join(gxe[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE",  "chiSqG", "chiSqGxE")], ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
         pval_GE = pchisq(chiSqGE, df = 1, lower.tail = F),
         logpval_GxE = -log10(pval_GxE)) %>% 
  arrange(pval_GE) %>% 
  dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval.GE, logpval.GxE)

# PRUNE?
prune <- do.call(rbind, lapply(list.files(path = "~/data/PriorityPruner/pval_gxe", full.names = T, pattern = "results"), fread, stringsAsFactors = F)) %>% 
  rename(ID = name) %>% 
  filter(selected == 1)

murcray <- inner_join(gxe[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE",  "chiSqG", "chiSqGxE")], ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(ID = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate),
         pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
         pval_GE = pchisq(chiSqGE, df = 1, lower.tail = F),
         logpval_GxE = -log10(pval_GxE)) %>% 
  arrange(pval_GE) %>% 
  dplyr::select(ID, Chromosome, Location, Reference, Alternate, betaGxE, pval_GE, logpval_GxE) %>% 
  filter(ID %in% prune$ID)

head(murcray)





#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 

m = 5413820
m = 554443
sizeBin0 = 5
alpha = 0.05
# alpha1 = 0.001 # first bin alpha level.. 
# N=72438; Nbio=100; ncp_GxE=0

nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction 
nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes 
endpointsBin = cumsum(sizeBin) # endpoints of the bins 
alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested 

# not sure what this is for...
# subs = cumsum(sizeBin) <= alpha1*m 
# alphaBin_subset = alpha / sizeBin * subs/ sum(subs==TRUE) 
# alphaBin_linear = alpha / sizeBin / nbins 
pv <- data.table(murcray)
head(pv)

## Function 11: Create Bin and significance level for weighted testing approach (needs to be sorted by p-value)
num = sizeBin[1] # initial bin size

rk.pv <- c(1:nrow(pv))
grp <- ceiling(log(rk.pv/num+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
pv[,Bin:=grp] # assigning group to the p values.. 
setkey(pv,Bin)
for(i in 1:max(grp)) 
{
  pv[J(i),Threshold:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
}

# see which are significant...
results <- pv %>% 
  mutate(logThreshold = -log10(Threshold))

results_sig <- results %>% 
  filter(logpval.GxE >= logThreshold)

summary(results$logpval.GxE)
results_suggestive <- filter(results, logpval.GxE >= 5)
table(results_suggestive$Bin)





# for PLOTS
num = nbins
num = 15
min.p = 12 # log scale...
last.sig = -log10(1.394709e-13)

# mapinfo is for scaling x axis plot
setnames(pv, c(8,9,10,3), c('p', 'grp','wt','MapInfo'))
pv[, y:= pv[,p]]
head(pv)

# create list where each component contains a bin
# log transfor bin alpha value
# create 'x', which is ordered by SNP BP sort of
glist<-list()
for(i in 1:num){
  t <- pv[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}
glist[[1]]
glist[[2]]

# trying to understand code above (just arrange BP by bins looks like)
# (Scale mapinfo for each Bin to range between 0.1-0.9 for neatness, and add a unit increase for successive Bin)
# x <- pv[1:5, MapInfo]
# normalized = (x-min(x))/(max(x)-min(x))
# normalized_scaled = 0.8 * normalized + 0.1
# x;normalized;normalized_scaled

# CREATE PLOT
head(glist[[1]]) # for reference

png("./figures/GxEScanR_asp_ref_noUKB_MurcrayGE_Weighted_PRUNED.png", width = 1280, height = 720)
color <- rep(c("blue","olivedrab4"),100)
plot(glist[[1]][,x], glist[[1]][,y], 
     col = "blue", 
     xlab="Bin # for step1 p-values", 
     ylab="-log10(step2 p-values)", 
     xlim=c(0,num), 
     ylim=c(0,min.p), 
     axes=F, pch=19, cex=0.5)
lines(glist[[1]][,x], glist[[1]][,ref],
      col = "black",lwd=2)

# the rest of the points ...  =|
# (adding to current plot..)
for(i in 2:num){
  points(glist[[i]][,x], glist[[i]][,y], 
         col = color[i], pch = 19, cex = 0.5)      
  lines(glist[[i]][,x], glist[[i]][,ref],
        col = "black",lwd = 2)
}

## the last bin..
## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
## let's only plot the first 15 bins for now, so change code a bit above. 
# points(glist[[num]][,x], glist[[num]][,y], 
#        col= color[num], pch = 19, cex = 0.5)
# lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 0.8)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=0.8)
title(main = "GxNSAIDs G|E Weighted Hypothesis Testing", sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.5, cex.sub = 1.2)

dev.off()








#--------------------------------------------------------#
# EDGxE ----
#--------------------------------------------------------#
edgxe <- inner_join(gxe[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE",  "chiSqG", "chiSqGxE")], ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
         chiSqG_chiSqGE = chiSqG + chiSqGE,
         pval_edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
         logpval_GxE = -log10(pval_GxE)) %>% 
  arrange(pval_edge) %>% 
  dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval_edge, logpval_GxE)

head(edgxe)


# PRUNE?
prune <- do.call(rbind, lapply(list.files(path = "~/data/PriorityPruner/pval_gxe", full.names = T, pattern = "results"), fread, stringsAsFactors = F)) %>% 
  rename(ID = name) %>% 
  filter(selected == 1)

edgxe <- inner_join(gxe[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE",  "chiSqG", "chiSqGxE")], ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(ID = paste0(Chromosome, ":", Location, "_", Reference, "_", Alternate), 
         pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
         chiSqG_chiSqGE = chiSqG + chiSqGE,
         pval_edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
         logpval_GxE = -log10(pval_GxE)) %>% 
  arrange(pval_edge) %>% 
  dplyr::select(ID, Chromosome, Location, Reference, Alternate, betaGxE, pval_edge, logpval_GxE) %>% 
  filter(ID %in% prune$ID)

head(edgxe)



#------ weighted hypothesis testing ------#
m = 5485386
m = 554443
sizeBin0 = 5
alpha = 0.05

nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction 
nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes 
endpointsBin = cumsum(sizeBin) # endpoints of the bins 
alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested 

pv <- data.table(edgxe)
head(pv)

## Function 11: Create Bin and significance level for weighted testing approach (needs to be sorted by p-value)
num = sizeBin[1] # initial bin size

rk.pv <- c(1:nrow(pv))
grp <- ceiling(log(rk.pv/num+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
pv[,Bin:=grp] # assigning group to the p values.. 
setkey(pv,Bin)
for(i in 1:max(grp)) 
{
  pv[J(i),Threshold:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
}

# see which SNPs are significant...
# results <- pv %>% 
#   mutate(logThreshold = -log10(Threshold))
# 
# results_sig <- results %>% 
#   filter(logpval.GxE >= logThreshold)
# 
# summary(results$logpval.GxE)
# results_suggestive <- filter(results, logpval.GxE >= 5)
# table(results_suggestive$Bin)

# for PLOTS
num = nbins
num = 15
min.p = 12 # log scale...

# mapinfo is for scaling x axis plot
setnames(pv, c(8,9,10,3), c('p', 'grp','wt','MapInfo'))
pv[, y:= pv[,p]]
head(pv)

# create list where each component contains a bin
# log transfor bin alpha value
# create 'x', which is ordered by SNP BP sort of
glist<-list()
for(i in 1:num){
  t <- pv[J(i)]
  t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
  t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
  glist[[i]]<-t
  rm(t)
}
glist[[1]]
glist[[2]]

# CREATE PLOT
head(glist[[1]]) # for reference

png("./figures/GxEScanR_asp_ref_noUKB_EDGxE_Weighted_PRUNED.png", width = 1280, height = 720)
color <- rep(c("blue","olivedrab4"),100)
plot(glist[[1]][,x], glist[[1]][,y], 
     col = "blue", 
     xlab="Bin # for step1 p-values", 
     ylab="-log10(step2 p-values)", 
     xlim=c(0,num), 
     ylim=c(0,min.p), 
     axes=F, pch=19, cex=0.5, cex.lab = 1.4)
lines(glist[[1]][,x], glist[[1]][,ref],
      col = "black",lwd=2)

# the rest of the points ...  =|
# (adding to current plot..)
for(i in 2:num){
  points(glist[[i]][,x], glist[[i]][,y], 
         col = color[i], pch = 19, cex = 0.5)      
  lines(glist[[i]][,x], glist[[i]][,ref],
        col = "black",lwd = 2)
}

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 1.5)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.5)
title(main = "GxNSAIDs EDGxE Weighted Hypothesis Testing (N = 57440)", sub = "iBin Size = 5, alpha = 0.05", cex.main = 2, cex.sub = 1.5)

dev.off()






table(prune$selected)


#--------------------------------------------------------#
# 3DF Procedure ----
#--------------------------------------------------------#
# start with weighted hypothesis testing
# need to obtain ranks from the step 1 statistic. Then test hypothesis using bins..
# make sure GxE p value is in the log scale for the plots
df3 <- inner_join(gxe, ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(stat_3df = chiSqG + chiSqGxE + chiSqGE, 
         pval_3df = pchisq(stat_3df, df = 3, lower.tail = F)) %>% 
  filter(!SNP %in% ukb_filter$SNP)

head(df3)

cases <- unique(df3$Cases)
controls <- unique(df3$Subjects) - unique(df3$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000
lambda <- getlambda3df(df3$pval_3df)
lambda1000 <- getlambda1000(lambda)

png("./figures/GxEScanR_asp_ref_chiSq3df_qq.png", height = 720, width = 1280)
qqman::qq(df3$pval_3df, 
          xlab = "Expected -log10(p)", 
          ylab = "Observed -log10(p)",
          main = "G_GxNSAIDs 3DF Results",
          #sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
          cex.main = 1.6, 
          cex.axis = 1.3, 
          cex.lab = 1.3,
          cex.sub = 1.3,
          col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3)
dev.off()


results_3df_EasyStrata <- df3 %>%
  rename(P = pval_3df) %>% 
  filter(P < 0.05) %>%
  dplyr::rename(CHR = Chromosome,
                BP = Location,
                A1 = Reference,
                A2 = Alternate) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, P)

write.table(results_3df_EasyStrata, file = "~/data/Annotations/EasyStrata_3DF_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("~/Dropbox/FIGI/Results/NSAIDS/results_asp_ref_EasyStrata_3DF.ecf")





#--------------------------------------------------------#
# 3DF MIAMI (vs G) ----
#--------------------------------------------------------#

names(results_G)
names(df3)

results_G_3df_miami <- inner_join(results_G, df3[, c("ID", "pval_3df")], by = "ID") %>% 
  rename(pval_G = P, 
         CHR = Chromosome, 
         BP = Location, 
         A1 = Reference, 
         A2 = Alternate) %>% 
  filter(pval_G < 0.05 | pval_3df < 0.05) %>% 
  dplyr::select(ID, CHR, BP, A1, A2, pval_G, pval_3df)

names(results_G_3df_miami)

write.table(results_G_3df_miami, file = "~/data/Annotations/EasyStrata_G_3DF_MIAMI_asp_ref_sex_age_pcs_studygxe_N72145.txt", quote = F, row.names = F, sep = '\t')


EasyStrata("~/Dropbox/FIGI/Results/NSAIDS/results_asp_ref_EasyStrata_Miami_G_3DF.ecf")




#--------------------------------------------------------#
# output for Prioritypruner ----
#--------------------------------------------------------#

# Required Columns:
#   name - Name of the SNP (e.g., rs222).
# chr - Name of the chromosome (e.g., 1, chr1). Chromosome X must be denoted as either 'X', 'chrX' or '23'. Chromosome Y and MT are not supported.
# pos - Physical position of the SNP.
# a1 - First allele of the SNP.
# a2 - Second allele of the SNP.
# p - P-value or other prioritization metric between 0 and 1. This is used for prioritizing the selection of SNPs, where lower numbers are prioritized.
# forceSelect - Flag indicating if the SNP should be selected (kept) regardless of its LD with other selected SNPs or other filtering criteria specified, such as MAF or design score (1=true, 0=false).
# designScore - Design score of the SNP (any positive real number). Can be filled in with a constant value (e.g., 1) if unknown.


# let's start with p values based on D|G statistic.. 

# didn't filter by UKB... don't filter here

results_G <- gxe %>%
  rename(name = SNP, 
         chr = Chromosome, 
         pos = Location,
         a1 = Reference,
         a2 = Alternate) %>% 
  mutate(p = pchisq(chiSqG, df = 1, lower.tail = F),
         forceSelect = 0,
         designScore = 1) %>%
  # filter(!ID %in% ukb_filter$ID,
  #        chr == 22) %>% 
  filter(chr == 22) %>% 
  dplyr::select(name, chr, pos, a1, a2, p, forceSelect, designScore)

write.table(results_G, file = "~/test.txt", quote = F, row.names = F)

# UGH need to update sex information on tfam file
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")

fam <- fread("~/final.tfam") %>% 
  inner_join(cov)

# also keep in mind multiallelics...
# simply removed them all. 
# code got lost because of rstudio crash, but easily write again. 




# results

wtfff <- fread("~/wtfff.results")


