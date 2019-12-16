set.seed(123456)
# A test on a correlation coefficient.
n<- 20
nrep<- 50000
price <- c(5.0, 4.8, 4.7, 4.0, 5.3, 4.1, 5.5, 4.7, 3.3, 4.0, 4.0, 4.6, 5.3, 3.0, 3.5, 3.9, 4.7, 5.0, 5.2, 4.6)
x<- quantity <- c(60, 59, 58, 47, 65, 48, 67, 70, 55, 63, 62, 65, 71, 56, 59, 60, 74, 77, 78, 62)
r_obs <- cor(price, quantity)
r_obs
r_perm <- vector()

for (i in 1:nrep) {
  
  y <- sample(price, n, replace = FALSE)
  r_perm <- c(r_perm,cor(x,y))
  
}

pval <- length(r_perm[r_perm >= r_obs])/nrep
pval

r_obs <- round(r_obs, digits = 3)
pval<-round(pval, digits=4)
hist(r_perm, breaks = 50, main = "Distribution of Sample Correlations",col="blue",freq=FALSE, xlab = expression(r[pi]))
legend(0.3,1.3,paste("p-value = ", pval))
legend(.5, 0.65, r_obs)
arrows(0.65,0.5,0.62, 0.05, col = "red", lwd=2, lty=1, length=0.1)



#=============================================================================#
# let's see if we can adapt this to gxe
# it's slightly different - you're not creating a distribution of a value like correlation coef
# instead you're counting number of significant hits @ some pre-determined cutoff value
#
# let's start with aspirin, since the significant finding is in bin #2
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(lmtest)
# ---- firstly, which bin is the significant finding in...
# it's in the fifth bin..

gxe <- readRDS('~/data/Results/aspirin/processed/FIGI_GxESet_aspirin_age_sex_pc3_studygxe_72269_results.rds')

tmp <- do.call(rbind, lapply(list.files("~/data/Results/aspirin/clump_controls/", full.names = T, pattern = "*.clumped"), fread, stringsAsFactors = F))
gxe_chiSqGxE_ld <- gxe %>% 
  filter(ID %in% tmp$SNP)
rm(tmp)

gxe_twostep <- format_twostep_data(dat = gxe_chiSqGxE_ld, 'chiSqG', 5, 0.05) %>%
  filter(grp == 2)
saveRDS(gxe_twostep, file = "~/permutation_test.rds", version = 2)



# ---- read in the dosage values you got from BinaryDosage GetSNPValues function
# filter vcfid included in the aspirin analysis

aspirin <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_72269_GLM.rds")

dosages <- readRDS("~/out.rds")

dosages_clean <- data.frame(do.call(cbind, lapply(dosages, function(x) x[which(rownames(x) %in% aspirin$vcfid), ]))) %>% 
  rownames_to_column(var = 'vcfid')
colnames(dosages_clean) <- gsub("_.*$", "", colnames(dosages_clean))


perm_data <- inner_join(aspirin, dosages_clean, by = 'vcfid') %>% 
  dplyr::select(outcome, age_ref_imp, sex, study_gxe, PC1, PC2, PC3, aspirin, colnames(dosages_clean))

saveRDS(perm_data, file = "~/permutation/perm_data.rds", version = 2)

#-----------------------------------------------------------------------------#
# actual permutation?
#-----------------------------------------------------------------------------#

# each permutation consists of fitting models for all 10 SNPs, then counting the number of SNPs that are 'significant'. 

# pseudo code

# for each permutation - 
#   apply columns glm with same covariates
#   count if snps p value < sig number 1
#   append number to vector 1
#   count if snps p value < sig number 2
#   append number to vector 2
#   
# do this for each snp in each iteration (testing if its smaller than sig levels)
# so.. internal loop counter



  
# lrtest for a single chromosome



glm_func <- function(y) glm(outcome ~ y * aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# run GLM + lrtest
basemodel <- map(perm_data[ , 10:19], function(x) glm_func_base(x))
intmodel <- map(perm_data[ , 10:19], function(x) glm_func(x))
results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)

results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
  tibble::rownames_to_column('SNP') %>%
  filter(!is.na(Chisq)) %>%
  mutate(SNP = gsub('.{2}$', '', SNP))



#-----------------------------------------------------------------------------#
# gather up results, create counts
#-----------------------------------------------------------------------------#
filelist <- list.files("~/permutation/", pattern = "permutation_", full.names = T)

get_counts_binp <- function(x) {
  sum(readRDS(x)[, c("Pr(>Chisq)")] < 0.00125)
}


filelist <- list.files("~/permutation/", pattern = "permutation_", full.names = T)
out <- do.call(c, lapply(filelist, get_counts_binp))
sum(out)
sum(out)/length(filelist)
hist(out)


get_counts_binp <- function(x) {
  sum(readRDS(x)[, c("Pr(>Chisq)")] < 0.0001197766)
}
out <- do.call(c, lapply(filelist, get_counts_binp))

hist(out)



#-----------------------------------------------------------------------------#
# gather up results for updated permutations
#-----------------------------------------------------------------------------#

x1 <- readRDS("~/permutation/permutation_ver2_p1_200.rds")
x2 <- readRDS("~/permutation/permutation_ver2_p201_400.rds")
x3 <- readRDS("~/permutation/permutation_ver2_p401_1000.rds")

x <- c(x1, x2, x3)

y <- lapply(x, function(x) arrange(x, `Pr(>Chisq)`))

y <- lapply(x, function(x) arrange(x, `Pr(>Chisq)`) %>% .[1, "Pr(>Chisq)"])

z <- -log10(do.call(c, y))

hist(z, xlim = c(0, 5), breaks = 40, main = "Distribution of -log10 p values")
abline(v = 3.921628, col = 'red')

sum(z < 3.921628) 
