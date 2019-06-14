#=======================================================================================#
# NSAIDS (NO UKB samples) - asp_ref results processing
# 02/28/2019
# 
# Generate Plots
# Implement 2-step methods 
#=======================================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
rm(list = ls())

# new annotations
fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>% 
	mutate(SNP = paste(Chr, Pos, sep = ":"))

# original results
gxe <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/gxe/", full.names = T, pattern = "results_GxE_asp_ref_sex_age_pc10_studygxe_"), fread, stringsAsFactors = F)) %>% 
	mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
	filter(!duplicated(ID))

# ge_only results (for 2-step methods, case only, control only)
colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
ge <- do.call(rbind, lapply(list.files(path = "~/data/Results/NSAIDS/asp_ref_noUKB/ge_only/", full.names = T, pattern = "results_GEOnly_asp_ref_sex_age_pc10_studygxe_"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
	mutate(ID = paste(SNP, Reference, Alternate, sep = "_")) %>% 
	filter(!duplicated(ID)) %>% dplyr::select(-drop)

# create merged file, only use this one for starting point for all over analyses
names(gxe)
names(ge)
results <- inner_join(gxe, ge[, c("betaGE", "chiSqGE", "betaCase", "ChiSqCase", "betaControl", "chiSqControl", "ID")], by = "ID")

rm(ge);rm(gxe);


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

#-----------------------------------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------------------------------#
results_G <- results %>%
	mutate(P = pchisq(chiSqG, df = 1, lower.tail = F))

# QQ Plot
wrapper_qq(df = results_G, 
					 p = 'P', 
					 title = 'G Main Effects QQ Plot\noutcome ~ dosage + age_ref_imp + sex + study_gxe + PC1-10 + asp_ref',
					 filename = "qq_plot_asp_ref_ukbfilter_G.png")


# Manhattan Plot
results_G_easystrata <- results_G %>% 
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
	dplyr::rename(CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
write.table(results_G_easystrata, file = "~/data/Annotations/EasyStrata_MarginalG_asp_ref_sex_age_pcs_studygxe_N72145_ukb_filter.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_G_LDAnnot_ukb_filter.ecf")


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


#-----------------------------------------------------------------------------------------------------#
# Marginal G Results - Table of Top Hits ----
#-----------------------------------------------------------------------------------------------------#
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
results_GxE <- results %>%
	mutate(P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
	dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaGxE, P)

cases <- unique(results_GxE$Cases)
controls <- unique(results_GxE$Subjects) - unique(results_GxE$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda(results_GxE$P)
lambda1000 <- getlambda1000(lambda)

CairoPNG(height = 720, width = 1280, file = "figures/GxEScanR_asp_ref_chiSqGxE_qq.png")
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

write.table(results_GxE_EasyStrata, file = "~/data/Annotations/EasyStrata_GxE_asp_ref_noUKB_sex_age_pcs_studygxe_N58109.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("easystrata/asp_ref_noUKB_GxE.ecf")


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


#----------------------------------------------------------#
# Future reference - get functional annotation of top hits:
results_GxE_filter <- filter(results_GxE, P <= 5E-8) %>% 
	mutate(v1 = "chromosome",
				 v2 = 1) %>% 
	dplyr::select(v1, CHR, BP, Reference, Alternate, v2)
write.table(results_GxE_filter, "~/Dropbox/tmp_getannot.txt", quote= F, row.names = F, col.names = F)

# annotation results from: https://snp-nexus.org/test/snpnexus_18324/results.html
snp_nexus <- fread("~/Dropbox/ncsnp_18324.txt") %>% 
	mutate(Chromosome = gsub("chr", "", Chromosome))





#--------------------------------------------------------#
# 2DF (G+GxE) results ----
#--------------------------------------------------------#
results_2df <- results %>%
	mutate(P = pchisq(chi2df, df = 2, lower.tail = F)) %>%
	dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P)

cases <- unique(results_2df$Cases)
controls <- unique(results_2df$Subjects) - unique(results_2df$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000

lambda <- getlambda2df(results_2df$P)
lambda1000 <- getlambda1000(lambda)

CairoPNG(height = 720, width = 1280, file = "figures/GxEScanR_asp_ref_chiSq2df_qq.png")
qqman::qq(results_2df$P, 
					xlab = "Expected -log10(p)", 
					ylab = "Observed -log10(p)",
					main = "G + GxNSAIDs (2DF) QQ Plot\noutcome ~ dosage + dosage*asp_ref + age_ref_imp + sex + study_gxe + PC1-10",
					cex.main = 1.8, 
					cex.axis = 1.5, 
					cex.lab = 1.5,
					cex.sub = 1.5,
					col = 'blue4')
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.5)
dev.off()


results_2df_EasyStrata <- results_2df %>%
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
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
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	rename(pval_G = P.x, 
				 pval_2df = P.y,
				 CHR = Chromosome, 
				 BP = Location, 
				 A1 = Reference, 
				 A2 = Alternate) %>% 
	filter(!(pval_2df > 0.05 & annot == 0)) %>% 
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
dg_gxe <- results %>% 
	mutate(pval_G = pchisq(chiSqG, df = 1, lower.tail = F),
				 pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F), 
				 logpval_GxE = -log10(pval_GxE)) %>% 
	arrange(pval_G) %>% 
	mutate(rank = seq(1, nrow(.))) %>% 
	dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval_G, logpval_GxE)

#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 
nrow(dg_gxe)
m = nrow(dg_gxe)
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
pv <- data.table(dg_gxe)
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

# for PLOTS
num = nbins
num = 15
min.p = 12 # log scale...
# last.sig = -log10(1.394709e-13)

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

## the last bin..
## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
## let's only plot the first 15 bins for now, so change code a bit above. 
# points(glist[[num]][,x], glist[[num]][,y], 
#        col= color[num], pch = 19, cex = 0.5)
# lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 1.4)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.4)
title(main = "GxNSAIDs DG|GxE (2-Step) Weighted Hypothesis Testing", sub = "Initial Bin Size = 5, alpha = 0.05", cex.main = 1.7, cex.sub = 1.4)

dev.off()



#--------------------------------------------------------#
# D|G 2-step Kooperberg subset testing ----
#--------------------------------------------------------#
dg_gxe <- results %>%
	mutate(P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
	dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P)

# QQ plot
wrapper_qq(df = dg_gxe, 
					 p = 'P', 
					 title = 'DG|GxE Step 1 (DG) QQ Plot',
					 filename = "qq_plot_asp_ref_ukbfilter_DG_GxE_Step1.png")

# Manhattan plot
# use alpha = 0.0016
dg_gxe_easystrata <- dg_gxe %>% 
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
	dplyr::rename(CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
write.table(dg_gxe_easystrata, file = "~/data/Annotations/kooperberg_step1_DG.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_G_ukb_filter_Kooperberg.ecf")


# step 2
alpha_2step <- 0.0016

dg_gxe_subset <- results %>% 
	mutate(pval_G = pchisq(chiSqG, df = 1, lower.tail = F),
				 pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F), 
				 logpval_GxE = -log10(pval_GxE)) %>% 
	filter(pval_G <= alpha_2step) %>% 
	arrange(Chromosome, Location)

adj1 <- alpha_2step * nrow(dg_gxe)
adj2 <- nrow(dg_gxe_subset)
adj1; adj2
bonferroni_p_liberal <- 0.05/adj1
bonferroni_p <- 0.05/adj2
bonferroni_p_liberal; bonferroni_p

# qq plot
wrapper_qq(df = dg_gxe_subset, 
					 p = 'pval_GxE', 
					 title = 'DG|GxE Step 2 (GxE) QQ Plot',
					 filename = "qq_plot_asp_ref_ukbfilter_DG_GxE_Step2.png")

# Manhattan
dg_gxe_subset_easystrata <- dg_gxe_subset %>% 
	dplyr::rename(P = pval_GxE, 
								CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
fudge <- c("1:1026708_C_A", 1, 1026708, "C", "A", 0.0405004, 4.5e-8)
dg_gxe_subset_easystrata <- rbind(dg_gxe_subset_easystrata, fudge)

write.table(dg_gxe_subset_easystrata, file = "~/data/Annotations/kooperberg_step2_GxE.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_G_ukb_filter_KooperbergSTEP2.ecf")



#--------------------------------------------------------#
# G|E 2-step Murcray ----
#--------------------------------------------------------#
ge_gxe <- results %>%
	mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
				 pval_GE = pchisq(chiSqGE, df = 1, lower.tail = F),
				 logpval_GxE = -log10(pval_GxE)) %>% 
	arrange(pval_GE) %>% 
	dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval_GE, logpval_GxE)


#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 

m = nrow(ge_gxe)
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
pv <- data.table(ge_gxe)
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

# for PLOTS
num = nbins
num = 15
min.p = 12 # log scale...

# mapinfo is for scaling x axis plot
setnames(pv, c(8,9,10,3), c('p', 'grp','wt','MapInfo'))
pv[, y:= pv[,p]]
head(pv)

# create list where each component contains a bin
# log transfor bin alpha value *** GOOD lol 
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

png("./figures/GxEScanR_asp_ref_MurcrayGE_Weighted.png", width = 1280, height = 720)
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

## the last bin..
## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
## let's only plot the first 15 bins for now, so change code a bit above. 
# points(glist[[num]][,x], glist[[num]][,y], 
#        col= color[num], pch = 19, cex = 0.5)
# lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 1.4)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.4)
title(main = "GxNSAIDs GE|GxE (2-Step) Weighted Hypothesis Testing", sub = "Initial Bin Size = 5, alpha = 0.05", cex.main = 1.7, cex.sub = 1.4)

dev.off()

#--------------------------------------------------------#
# G|E 2-step Murcray subset testing ----
#--------------------------------------------------------#
ge_gxe <- results %>%
	mutate(P = pchisq(chiSqGE, df = 1, lower.tail = F)) %>%
	dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, betaG, P)

## qq plot step 1
wrapper_qq(df = ge_gxe, 
					 p = 'P', 
					 title = 'GE|GxE Step 1 (GE) QQ Plot',
					 filename = "qq_plot_asp_ref_ukbfilter_GE_GxE_Step1.png")

# manhattan step 1
# use 9e-4
ge_gxe_easystrata <- ge_gxe %>% 
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
	dplyr::rename(CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
write.table(ge_gxe_easystrata, file = "~/data/Annotations/murcray_step1_GE.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_GE_ukb_filter_Murcray.ecf")



# step 2
alpha_2step <- 9e-4

ge_gxe_subset <- results %>% 
	mutate(pval_GE = pchisq(chiSqGE, df = 1, lower.tail = F),
				 pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F), 
				 logpval_GxE = -log10(pval_GxE)) %>% 
	filter(pval_GE <= alpha_2step) %>% 
	arrange(Chromosome, Location)

adj1 <- alpha_2step * nrow(dg_gxe)
adj2 <- nrow(dg_gxe_subset)
adj1; adj2
bonferroni_p_liberal <- 0.05/adj1
bonferroni_p <- 0.05/adj2
bonferroni_p_liberal; bonferroni_p

# qq plot
wrapper_qq(df = ge_gxe_subset, 
					 p = 'pval_GxE', 
					 title = 'GE|GxE Step 2 (GxE) QQ Plot',
					 filename = "qq_plot_asp_ref_ukbfilter_GE_GxE_Step2.png")

# Manhattan
ge_gxe_subset_easystrata <- ge_gxe_subset %>% 
	dplyr::rename(P = pval_GxE, 
								CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
head(ge_gxe_subset_easystrata)
fudge <- c("1:2889106_G_A", 1, 2889106, "C", "A", 0.0405004, 4.5e-8)
ge_gxe_subset_easystrata <- rbind(ge_gxe_subset_easystrata, fudge)

write.table(ge_gxe_subset_easystrata, file = "~/data/Annotations/murcray_step2_GxE.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_G_ukb_filter_MurcraySTEP2.ecf")

# check
df <- mutate(ge_gxe_subset_easystrata, 
						 CHR = as.numeric(CHR),
						 BP = as.numeric(BP))
manhattan(ge_gxe_subset_easystrata, p = 'P')


#--------------------------------------------------------#
# EDGxE ----
#--------------------------------------------------------#

edgxe <- results %>% 
	mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
				 chiSqG_chiSqGE = chiSqG + chiSqGE,
				 pval_edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
				 logpval_GxE = -log10(pval_GxE)) %>% 
	arrange(pval_edge) %>% 
	dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, pval_edge, logpval_GxE)

head(edgxe)


#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 

m = nrow(edgxe)
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

png("./figures/GxEScanR_asp_ref_EDGxE_Weighted.png", width = 1280, height = 720)
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

axis(1, at = c(-1.5, seq(0.5, num-0.5, 1)), label = c(0, seq(1, num, 1)), cex.axis = 1.4)
axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1.4)
title(main = "GxNSAIDs EDGxE (2-Step) Weighted Hypothesis Testing", sub = "Initial Bin Size = 5, alpha = 0.05", cex.main = 1.7, cex.sub = 1.4)

dev.off()



#--------------------------------------------------------#
# EDGxE Subset Testing ----
#--------------------------------------------------------#

edgxe <- results %>%
	mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
				 chiSqG_chiSqGE = chiSqG + chiSqGE,
				 pval_edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
				 logpval_GxE = -log10(pval_GxE)) %>% 
	arrange(pval_edge)

## qq plot step 1
wrapper2df_qq(df = edgxe, 
							p = 'pval_edge', 
							title = 'EDGxE Step 1 (G_GE) QQ Plot',
							filename = "qq_plot_asp_ref_ukbfilter_EDGxE_Step1.png")

# manhattan step 1
# use 6.3e-5
edgxe_easystrata <- edgxe %>% 
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	dplyr::rename(P = pval_edge, 
								CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
write.table(edgxe_easystrata, file = "~/data/Annotations/edgxe_step1_GE.txt", quote = F, row.names = F, sep = '\t')
EasyStrata("results_asp_ref_EasyStrata_GE_ukb_filter_edgxe.ecf")


## Step 2


# step 2
alpha_2step <- 6.3e-5
edgxe_subset <- results %>%
	mutate(pval_GxE = pchisq(chiSqGxE, df = 1, lower.tail = F),
				 chiSqG_chiSqGE = chiSqG + chiSqGE,
				 pval_edge = pchisq(chiSqG_chiSqGE, df = 2, lower.tail = F),
				 logpval_GxE = -log10(pval_GxE)) %>% 
	filter(pval_edge <= alpha_2step) %>% 
	arrange(Chromosome, Location)

adj1 <- alpha_2step * nrow(edgxe)
adj2 <- nrow(edgxe_subset)
adj1; adj2
bonferroni_p_liberal <- 0.05/adj1
bonferroni_p <- 0.05/adj2
bonferroni_p_liberal; bonferroni_p

# qq plot
wrapper_qq(df = edgxe_subset, 
					 p = 'pval_GxE', 
					 title = 'EDGxE Step 2 (GxE) QQ Plot',
					 filename = "qq_plot_asp_ref_ukbfilter_EDGxE_Step2.png")

# Manhattan
edgxe_subset_easystrata <- edgxe_subset %>% 
	dplyr::rename(P = pval_GxE, 
								CHR = Chromosome,
								BP = Location, 
								A1 = Reference, 
								A2 = Alternate) %>% 
	dplyr::select(ID, CHR, BP, A1, A2, betaG, P)
head(edgxe_subset_easystrata)
fudge <- c("1:38397370_C_T", 1, 38397370, "C", "A", 0.0405004, 4.5e-8)
edgxe_subset_easystrata <- rbind(edgxe_subset_easystrata, fudge)

write.table(edgxe_subset_easystrata, file = "~/data/Annotations/edgxe_step2_GxE.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("results_asp_ref_EasyStrata_GE_ukb_filter_edgxeSTEP2.ecf")

# check
df <- mutate(edgxe_subset_easystrata, 
						 CHR = as.numeric(CHR),
						 BP = as.numeric(BP),
						 P = as.numeric(P))
manhattan(df, p = 'P')











# filter markers at edge stat significant @ 6.3x10-5
adj1 <- 6.3e-5*5413820 # 341

# N = 5375
edgxe_subset_testing <- edgxe_filter %>%
	filter(pval_edge <= 6.3e-5) %>%
	arrange(Chromosome, Location)

adj2 <- nrow(edgxe_subset_testing) # 5375

edgxe_bonferroni_p <- 0.05/adj1
edgxe_bonferroni_p_true <- 0.05/adj2

table(edgxe_subset_testing$Chromosome)

# step2 plots
cases <- unique(edgxe_subset_testing$Cases)
controls <- unique(edgxe_subset_testing$Subjects) - unique(edgxe_subset_testing$Cases)
total <- cases + controls
cases1000 <- (cases/total) * 1000
controls1000 <- (controls/total) * 1000
lambda <- getlambda(edgxe_subset_testing$pval_GxE)


png("./figures/GxEScanR_asp_ref_EDGxE_Step2_qq.png", height = 720, width = 1280)
qqman::qq(edgxe_subset_testing$pval_GxE, 
					xlab = "Expected -log10(p)", 
					ylab = "Observed -log10(p)",
					main = "EDGxE Statistic - Step 2 GxE Pval (N=5375) QQ Plot",
					#sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))),
					cex.main = 1.6, 
					cex.axis = 1.3, 
					cex.lab = 1.3,
					cex.sub = 1.3,
					col = 'blue4') # ~~ adds space lol
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4))), cex.sub = 1.3)
dev.off()

edgxe_subset_testing_easystrata <- edgxe_subset_testing %>%
	# filter(pval_GxE < 0.05) %>% 
	dplyr::rename(CHR = Chromosome,
								BP = Location,
								A1 = Reference,
								A2 = Alternate) %>%
	dplyr::select(SNP, CHR, BP, A1, A2, betaGxE, pval_GxE)

write.table(edgxe_subset_testing_easystrata, file = "~/data/Annotations/EasyStrata_GxE_asp_ref_sex_age_pcs_studygxe_N72145_EDGxE.txt", quote = F, row.names = F, sep = '\t')

EasyStrata("results_asp_ref_EasyStrata_EDGxE_SubsetTesting.ecf")



# quick and dirty version
# manhattan(edgxe_subset_testing, chr = "Chromosome", bp = "Location", p = "pval.GxE", genomewideline = -log10(edgxe_bonferroni_p), suggestiveline = F)

edgxe_sig <- edgxe_subset_testing %>% 
	filter(pval_GxE < edgxe_bonferroni_p)

x <- fread("./rds/GxE_asp_ref_EDGxE_subsetTesting_sig_snplist_63286.txt") %>% arrange(Position) %>% inner_join(edgxe_sig[, c("Location", "pval_edge", "betaGxE", "pval_GxE")], by = c("Position" = "Location"))

write.table(x, file = "./rds/GxE_asp_ref_EDGxE_subsetTesting_sig_snplist_final.txt", quote = F, row.names = F)


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


CairoPNG(height = 720, width = 1280, file = "figures/GxEScanR_asp_ref_chiSq3df_qq.png")
qqman::qq(df3$pval_3df, 
					xlab = "Expected -log10(p)", 
					ylab = "Observed -log10(p)",
					main = "G + GxNSAIDs + G|E (3DF) QQ Plot",
					cex.main = 1.8, 
					cex.axis = 1.5, 
					cex.lab = 1.5,
					cex.sub = 1.5,
					col = 'blue4')
par(adj = 1)
title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.5)
dev.off()


mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
	
	
	results_3df_EasyStrata <- df3 %>%
	rename(P = pval_3df) %>% 
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	filter(!(P > 0.05 & annot == 0)) %>% 
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
	mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
	rename(pval_G = P, 
				 CHR = Chromosome, 
				 BP = Location, 
				 A1 = Reference, 
				 A2 = Alternate) %>% 
	filter(!(pval_3df > 0.05 & annot == 0)) %>% 
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



#--------------------------------------------------------#
# GWAS95 Replicates Table (main) ----
#--------------------------------------------------------#
results_G <- results %>%
	mutate(P = pchisq(chiSqG, df = 1, lower.tail = F)) %>% 
	filter(P < 5e-8)

annot <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata.tsv")
names(annot)
table(annot$PREVIOUSLY_REPORTED)

z <- inner_join(results_G, annot, by = "ID") %>% 
	arrange(Chromosome, Location) %>% 
	dplyr::select(ID, LOCUS, PREVIOUSLY_REPORTED, RSID, AUTHOR_FIRST_REPORTED, YEAR_FIRST_REPORTED) %>% 
	arrange(PREVIOUSLY_REPORTED)

saveRDS(z, file = "rds/marginal_g_previoushits_replicates.rds")

library(kableExtra)

kable(z, format = 'html') %>%
	kable_styling(bootstrap_options = c("condensed"), full_width = F, position = "left") %>%
	as_image(file = "figures/test.png")



#--------------------------------------------------------#
# p value plots (2df vs G, 3df vs G) ----
#--------------------------------------------------------#

results_2df_pval <- results %>%
	mutate(pval_G = pchisq(chiSqG, df = 1, lower.tail = F),
				 pval_2df = pchisq(chi2df, df = 2, lower.tail = F)) %>%
	dplyr::select(SNP, ID, Chromosome, Location, Subjects, Cases, Reference, Alternate, pval_G, pval_2df) %>% 
	filter(pval_G < 5e-8)

huh <- results_2df_pval %>% 
	filter(pval_2df < pval_G)


png("./figures/pval_comparison_2df_G_sig_only.png", width = 1280, height = 720)
plot(-log10(results_2df_pval$pval_G), -log10(results_2df_pval$pval_2df), cex = 0.5, ylab = "2DF log P values", xlab = "G Main Effects log P values", main = "2DF vs G Main Effects log P value (pval_G < 5e-8)", cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
abline(a = 0, b = 1, col = 'red')
dev.off()



results_3df_pval <- results %>%
	mutate(stat_3df = chiSqG + chiSqGxE + chiSqGE, 
				 pval_3df = pchisq(stat_3df, df = 3, lower.tail = F),
				 pval_G = pchisq(chiSqG, df = 1, lower.tail = F),
				 stat_test = chi2df + chiSqGE,
				 pval_test = pchisq(stat_test, df = 3, lower.tail = F)) %>%
	filter(pval_G < 5e-8)

huh <- results_3df_pval %>% 
	filter(pval_3df < pval_G)

png("./figures/pval_comparison_3df_G_sig_only.png", width = 1280, height = 720)
plot(-log10(results_3df_pval$pval_G), -log10(results_3df_pval$pval_3df), cex = 0.5, ylab = "3DF log P values", xlab = "G Main Effects log P values", main = "3DF vs G Main Effects log P value (pval_G < 5e-8)", cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
abline(a = 0, b = 1, col = 'red')
dev.off()
