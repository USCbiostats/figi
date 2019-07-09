#=============================================================================#
# Estimating Rsq
#
# 1) var(dosage) / var(expected HWE) --- 'Medland'
# 2) var(imputedAltAlleleProb) / var(expected) --- 'Minimac'
#=============================================================================#
library(BinaryDosage)

# ------ Input Data ------
# mecc imputation batch, chr22
mecc_chr22_index <- GetBDoseInfo("/home/rak/data/mecc/mecc_chr22.bdose", index = T)
mecc_chr22_info <- fread("/home/rak/data/mecc/chr22.info")

# ------ Estimate Rsq using ratio of variances (medland) ------

# These values match those in the bdose object data.frame
GetAlternateAlleleFrequencies(mecc_chr22_index, 1:10)

# first, ensure you're getting the right variance values..
# (they match, very well minus rounding)
# mecc_chr22_dose <- GetSNPValues(mecc_chr22_index, 1:10000, geneProb = F)
# tmp1 <- apply(mecc_chr22_dose, 2, var)
# tmp2 <- GetAlternateAlleleVariance(mecc_chr22_index, 1:10000)
# plot(tmp1, tmp2)

# actual calculation
rsq_medland <- GetRsqEstimateMinimac(mecc_chr22_index, 1:mecc_chr22_index$NumSNPs)





# Start with common SNPs MAF > 0.05
snps_common <- which(mecc_chr22_info$MAF > 0.05)

a <- GetAlternateAlleleFrequencies(mecc_chr22_index, snps_common)
b <- GetAlternateAlleleVariance(mecc_chr22_index, snps_common)

e <- 2*a*(1-a)
Rsq_est <- b/e
Rsq_info <- as.numeric(mecc_chr22_info[snps_common, ]$Rsq)

png("/home/rak/Dropbox/FIGI/Code/Rsq_Estimation/mecc_chr22_Rsq_comparison_medland_common.png", height = 720, width = 1280)
plot(Rsq_info, Rsq_est, main = "Comparison of Imputation Quality Metric - MAF > 5%", xlab = "Original HRC/Minimac", ylab = "Dosage Variance / Expected Variance (HWE)")
abline(0,1, col = 'red')
dev.off()



# Rare SNPs MAF < 0.01
snps_rare <- which(mecc_chr22_info$MAF < 0.01)

a <- GetAlternateAlleleFrequencies(mecc_chr22_index, snps_rare)
b <- GetAlternateAlleleVariance(mecc_chr22_index, snps_rare)

e <- 2*a*(1-a)
Rsq_est <- b/e
Rsq_info <- as.numeric(mecc_chr22_info[snps_rare, ]$Rsq)

png("/home/rak/Dropbox/FIGI/Code/Rsq_Estimation/mecc_chr22_Rsq_comparison_medland_rare.png", height = 720, width = 1280)
plot(Rsq_info, Rsq_est, main = "Comparison of Imputation Quality Metric - MAF < 0.1%", xlab = "Original HRC/Minimac", ylab = "Dosage Variance / Expected Variance (HWE)")
abline(0,1, col = 'red')
dev.off()



# ------ Estimate Rsq using ratio of variances, alt allele probability (minimac) ------

# Start with common SNPs MAF > 0.05
snps_common <- which(mecc_chr22_info$MAF > 0.05)

Rsq_est <- GetRsqEstimateMinimac(mecc_chr22_index, SNPs = snps_common)
Rsq_info <- as.numeric(mecc_chr22_info[snps_common, ]$Rsq)

png("/home/rak/Dropbox/FIGI/Code/Rsq_Estimation/mecc_chr22_Rsq_comparison_minimac_common.png", height = 720, width = 1280)
plot(Rsq_info, Rsq_est, main = "Comparison of Imputation Quality Metric - MAF > 5%", xlab = "Original HRC/Minimac", ylab = "AA Probability Variance / Expected Variance",
     xlim = c(0,1), ylim = c(0,1))
abline(0,1, col = 'red')
dev.off()



# Rare SNPs MAF < 0.01
snps_rare <- which(mecc_chr22_info$MAF < 0.01)
Rsq_est <- GetRsqEstimateMinimac(mecc_chr22_index, SNPs = snps_rare)
Rsq_info <- as.numeric(mecc_chr22_info[snps_rare, ]$Rsq)


png("/home/rak/Dropbox/FIGI/Code/Rsq_Estimation/mecc_chr22_Rsq_comparison_minimac_rare.png", height = 720, width = 1280)
plot(Rsq_info, Rsq_est, main = "Comparison of Imputation Quality Metric - MAF < 0.1%", xlab = "Original HRC/Minimac", ylab = "AA Probability Variance / Expected Variance",
     xlim = c(0,1), ylim = c(0,1))
abline(0,1, col = 'red')
dev.off()



# mostly MAF = 0 etc. none good quality markers.
bad <- data.frame(mecc_chr22_info[snps_rare, ], Rsq_est) %>% filter(is.na(Rsq_est))
bad_snps <- which(mecc_chr22_info$SNP %in% bad$SNP)








# ------ Explore drops ------

# 22:27796846

# explore bad markers
snps_common <- which(mecc_chr22_info$MAF > 0.05)
Rsq_est <- GetRsqEstimateMinimac(mecc_chr22_index, SNPs = snps_common)
Rsq_info <- as.numeric(mecc_chr22_info[snps_common, ]$Rsq)
plot(Rsq_info, Rsq_est, main = "Comparison of Imputation Quality Metric - MAF > 5%", xlab = "Original HRC/Minimac", ylab = "AA Probability Variance / Expected Variance",
     xlim = c(0,1), ylim = c(0,1))

wtf <- test(mecc_chr22_index, SNPs = c('22:27796846'))


bad <- data.frame(mecc_chr22_info[snps_common, ], Rsq_est) %>% filter(is.na(Rsq_est))
bad_snps <- which(mecc_chr22_info$SNP %in% bad$SNP)


wtf <- test(mecc_chr22_index, SNPs = c("22:27796846"))
x <- sum(GetSNPValues(mecc_chr22_index, SNPs = c("22:27796846"), geneProb = F)[, 1])



x <- GetSNPValues(mecc_chr22_index, bad_snps, geneProb = F)[, 1:100]
apply(x, 2, sum)
# choose one, so calculation by hand


badddd <- GetRsqEstimateMinimac(mecc_chr22_index, SNPs = c("22:27796846"))
