#-----------------------------------------------------------------------------#
# actual permutation?
#-----------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
library(figifs)
library(lmtest)
set.seed(2019)

# convenience functions
glm_func <- function(y) glm(outcome ~ y * aspirin_p + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
glm_func_base <- function(y) glm(outcome ~ y + aspirin_p + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = perm_data, family = binomial)
run_lrtest <- function(x,y) lrtest(x, y)

# data set
perm_data <- readRDS("~/permutation/perm_data.rds")

# loop
for (i in 1:10) {
  
  perm_data$aspirin_p <- sample(perm_data$aspirin, 72269, replace = FALSE)
  
  basemodel <- map(perm_data[ , 10:19], function(x) glm_func_base(x))
  intmodel <- map(perm_data[ , 10:19], function(x) glm_func(x))
  results_gxe_glm <- mapply(run_lrtest, basemodel, intmodel, SIMPLIFY = F)
  
  results_gxe_out <- do.call(rbind, results_gxe_glm) %>%
    tibble::rownames_to_column('SNP') %>%
    filter(!is.na(Chisq)) %>%
    mutate(SNP = gsub('.{2}$', '', SNP))
  
  saveRDS(results_gxe_out, file = "~/permutation/permutation_", i, ".rds", version = 2)
}


