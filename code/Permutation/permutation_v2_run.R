#-----------------------------------------------------------------------------#
# actual permutation?
#-----------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
library(figifs)
library(lmtest)
library(foreach)
library(doParallel)
library(doMC)
library(doSNOW)

dat <- readRDS("~/data/permutation_aspirin/perm_data_clumped.rds")
registerDoMC(cores=2)

results_gxe_out <- foreach(i=1:1000) %dopar% {
  
  N <- length(dat$aspirin)
  
  sample_perm <- sample(N, N, replace = FALSE)
  
  glm_func <- function(y) glm(dat[,'outcome'] ~ y * dat[sample_perm, 'aspirin'] + dat[, 'nsaids'] + dat[, 'age'] + dat[,'sex'] + dat[, 'PC1'] + dat[,'PC2'] + dat[, 'PC3'] + dat[, 'study'], family = binomial)
  glm_func_base <- function(y) glm(dat[,'outcome'] ~ y + dat[sample_perm, 'aspirin'] + dat[, 'nsaids'] + dat[, 'age'] + dat[,'sex'] + dat[, 'PC1'] + dat[,'PC2'] + dat[, 'PC3'] + dat[, 'study'], family = binomial)
  run_lrtest <- function(x,y) lrtest(x, y)
  
  basemodel <- map(dat[ , 11:30], function(x) glm_func_base(x))
  intmodel <- map(dat[ , 11:30], function(x) glm_func(x))
  results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)
  
  do.call(rbind, results_gxe_glm) %>%
    tibble::rownames_to_column('SNP') %>%
    filter(!is.na(Chisq)) %>%
    mutate(SNP = gsub('.{2}$', '', SNP))
  
  saveRDS(results_gxe_out[[1]], file = paste0("~/permutation/permutation_perm_data_clump_number_", i, ".rds"), version = 2)
}

# saveRDS(results_gxe_out, file = paste0("~/permutation/permutation_ver2_p401_1000.rds"), version = 2)
