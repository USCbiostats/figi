#=============================================================================#
# GxEScan post-hoc analysis
# 06/01/2019
# 
# Commonly used functions to generate plots etc
# (I know you tried this before, but this time it might be worth it)
# (make sure john's output absolutely never changes again)
# BAD PRACTICE - RELIES ON GLOBAL ENVIRONMENTS: 
# - E: exposure (character)
# - covs: covariates (vector)
# - N: samplesize
#
# keep in mind that for now, you're only creating a set number of plots:
# - G, GxE, 2DF, 3DF, GE, Case, Control. don't go nuts with overcomplication
#=============================================================================#
library(tidyverse)
library(data.table)
library(EasyStrata)
library(qqman)


# ------ Create Bar Plots ------
# no choice but to specify factors levels etc
# for now creating barplots by study_gxe
createBarPlot <- function(data, outcome, exposure, flevels, flabels, fcolors) {
  
  outcome_quo <- enquo(outcome)
  exposure_quo <- enquo(exposure)
  
  data <- data %>%
    dplyr::mutate(!! exposure_quo := factor(!! exposure_quo, exclude = NULL, levels = flevels, labels = flabels)) # %>%
  # dplyr::select(!! outcome_quo, study_gxe, !! exposure_quo) %>%
  # group_by(study_gxe)
  
  ggplot(data=data, aes(x=outcome)) +
    geom_bar(aes(fill = !! exposure_quo), position = 'fill') +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    scale_fill_manual(values = fcolors) +
    facet_wrap(vars(study_gxe), ncol = 5)
}


# ------ Create Tables (Table 1) ------

# function very limited but works for simple cases (like by outcome, with p values)
# make sure 'covarCat' is a vector of variables names (strings) to include in the table
# (covariates has all vars you want to include in table)
# only filter will be !(is.na(exposure))
createTable1 <- function(data, outcome, exposure, covariates, covarCat) {
  
  outcome_quo <- enquo(outcome)
  exposure_quo <- enquo(exposure)
  
  # temporary data
  # create 3 level outcome factor variable
  # make sure categorical variables are FACTORS
  data <- data %>%
    dplyr::filter(!is.na(!! exposure_quo)) %>%
    dplyr::mutate(table1_outcome = factor(!! outcome_quo, levels = c(0,1,2), labels=c("Control", "Case", "P-value"))) %>%
    dplyr::mutate_at(., setdiff(covariates, covarCat), as.numeric) %>% 
    dplyr::mutate_at(., covarCat, as.factor)
  
  my.render.cont <- function(x) {
    with(stats.apply.rounding(stats.default(x), digits=2),
         c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
  }
  
  my.render.cat <- function(x) {
    c("", sapply(stats.default(x),
                 function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
  }
  
  # p values
  # hacky - requires an outcome factor with third level = p value. 
  # Also need to specify values explicitly in the 'rndr' function
  rndr <- function(x, name, ...) {
    if (length(x) == 0) {
      y <- data[[name]]
      s <- rep("", length(render.default(x=y, name=name, ...)))
      if (is.numeric(y)) {
        p <- t.test(y ~ data$table1_outcome)$p.value
      } else {
        p <- chisq.test(table(y, droplevels(data$table1_outcome)))$p.value
      }
      s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
      s
    } else {
      render.default(x=x, name=name, ...)
    }
  }
  
  rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
  }
  
  # (paste0("~ ", get_expr(exposure_quo), "+", paste(covariates, collapse = "+"), " | table1_outcome"))
  
  table1(as.formula(paste0("~ ", get_expr(exposure_quo), "+", paste(covariates, collapse = "+"), "| table1_outcome")),
         data=data,
         render.continuous=my.render.cont,
         render.categorical=my.render.cat,
         render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)
}

# x <- createTable1(cov, outcome, aspirin, c("sex", "age_ref_imp"), c("sex"))





# ------ Run meta-analysis ------
# conducted by study/platform (study_gxe)

# quick function to format data for meta analysis
# might have to handle situations where cell counts are way too small (N < 10?)
create_data_4_meta <- function(data, outcome, exposure) {
  
  outcome_quo <- enquo(outcome)
  exposure_quo <- enquo(exposure)
  
  # filter cases missing E
  # filter case/control only studies
  data <- data %>%
    dplyr::filter(!is.na(!! exposure_quo))
  
  drops <- data.frame(table(data$study_gxe, dplyr::pull(data, !! outcome_quo))) %>%
    dplyr::filter(Freq == 0)
  exclude_studies <- as.vector(unique(drops$Var1))
  
  data <- filter(data, !(study_gxe %in% exclude_studies))
  
  # check for cell sizes, drop if study_gxe vs outcome cell N < 10
  drops2 <- data.frame(table(data[, 'study_gxe'], data[, quo_name(outcome_quo)])) %>%
    filter(Freq < 10)
  exclude_studies2 <- as.vector(unique(drops2$Var1))
  data <- filter(data, !(study_gxe %in% exclude_studies2))
  
}

# x <- create_data_4_meta(cov, outcome, aspirin)
# x <- create_data_4_meta(tmp, outcome, aspirin)




# get outcome counts by study_gxe
# outputs a dataframe (study, control, case, N total)
getCounts_byOutcome <- function(data, outcome, exposure) {
  
  # need to read up NSE very soon
  data <- create_data_4_meta(data, !! enquo(outcome), !! enquo(exposure))
  
  # get counts for rmeta function
  tmp_counts <- data %>%
    dplyr::group_by(!! enquo(outcome), study_gxe) %>%
    dplyr::summarize(count = n()) %>%
    tidyr::spread(key = !! enquo(outcome), value = count) %>%
    dplyr::rename(Control = `0`, Case = `1`) %>%
    dplyr::mutate(N = Control + Case)
}

# x <- getCounts_byOutcome(cov, outcome, aspirin)


# run GLM by study_gxe
# just use string formula_txt as input
# make sure only variables that you need indicators for are FACTORS (otherwise just continuous)
getGLM_byGroup <- function(data, outcome, exposure, covars) {
  
  # need to read up NSE very soon
  data <- create_data_4_meta(data, !! enquo(outcome), !! enquo(exposure)) %>% 
    group_by(study_gxe)
  
  # quo_name(exposure_quo)
  # rlang::sym(rlang::quo_name(exposure_quo))
  
  results_beta <- do(data, tidy(glm(as.formula(paste0(quo_name(enquo(outcome)), " ~ ", quo_name(enquo(exposure)), "+", paste(covars, collapse = "+"))), data = . , family = 'binomial')))
  results_ci   <- do(data, confint_tidy(glm(as.formula(paste0(quo_name(enquo(outcome)), " ~ ", quo_name(enquo(exposure)), "+", paste(covars, collapse = "+"))), data = . , family = 'binomial')))
  
  results <- bind_cols(results_beta, results_ci) %>%
    ungroup %>%
    dplyr::filter(grepl(quo_name(enquo(exposure)), term)) %>%
    dplyr::select(-study_gxe1)
  results
  
}
# x <- getGLM_byGroup(cov, outcome, aspirin, c("age_ref_imp", "sex"))

# x <- getCounts_byOutcome(tmp, outcome, aspirin)
# x <- getGLM_byGroup(tmp, outcome, aspirin, c("age_ref_imp", "sex"))
# table(x$study_gxe, x$outc)



# quick experiment
# 
# f <- function(a,b,c) {
#   # enquo(c) # returns quosure
#   # c # evaluates it
#   # get_expr(c) # evaluates it
#   # get_expr(enquo(c)) # class 'name'
#   # substitute(c) # same as above
#   # deparse(quote(c)) # returns 'c'
#   # deparse(substitute(c)) # "aspirin"
#   # quo_name(enquo(c)) # "aspirin"
# }
# 
# f(cov, outcome, aspirin)
# class(f(cov, outcome, aspirin))



# take outputs from above to perform meta analysis and plot forest plots
# just insert title yourself, make sure to include model
run_meta_analysis_create_forestplot <- function(counts_df, results_df, title) {
  
  # run meta analysis (summary OR + ConfInt)
  tmp_meta <- meta.summaries(results_df$estimate, results_df$std.error, method = 'random') # returns list object
  tmp_meta_b_ci <- c(tmp_meta[[3]], tmp_meta[[3]]-(tmp_meta[[4]]*1.96), tmp_meta[[3]]+(tmp_meta[[4]]*1.96))
  
  # format forestplot text (counts, study specific ORs)
  tmp_fp <- dplyr::select(results_df, estimate, conf.low, conf.high)
  tmp_fp <- rbind(rep(NA, 3), tmp_fp, rep(NA, 3), tmp_meta_b_ci) # add meta summaries, need conf.int for forestplot
  tmp_label_summary <- c("Summary", sum(counts_df$N), sum(counts_df$Case), round(tmp_meta[[3]], 2), round(tmp_meta[[4]], 3), round(exp(tmp_meta[[3]]), 2), round(tmp_meta[[5]][2], 4))
  
  tmp_label <- inner_join(counts_df, results_df, by = 'study_gxe') %>%
    mutate(Study = as.character(study_gxe),
           OR = round(exp(estimate), 2),
           Pvalue = round(p.value, 4),
           SE = round(std.error, 3),
           beta = round(estimate, 2)) %>%
    dplyr::select(Study, N, Case, beta, SE, OR, Pvalue)
  
  
  tmp_label2 <- rbind(colnames(tmp_label), tmp_label, rep(NA, 7), tmp_label_summary)
  rmeta::forestplot(label=tmp_label2,
                    as.numeric(tmp_fp$estimate), as.numeric(tmp_fp$conf.low), as.numeric(tmp_fp$conf.high),
                    is.summary = c(T, rep(F, nrow(tmp_label2)-2), T),
                    col=meta.colors(box="royalblue",line="darkblue",zero='gray0', summary="royalblue"),
                    clip=c(-1.5,1.5),
                    xlog=T,
                    xlab='logOR',
                    zero=0, title(main = title, line = 0))
  mtext(paste("het.pval=", formatC(tmp_meta$het[3], format='e', digits=2)), side = 1, line = 1)
}
# 
# x <- getCounts_byOutcome(cov, outcome, aspirin)
# y <- getGLM_byGroup(cov, outcome, aspirin, c("age_ref_imp", "sex"))
# run_meta_analysis_create_forestplot(x, y, "testing")







# getMeta <- function(data, exposure, covars) {
#   # create a tmp df that only includes covariates, complete case, and remove samples if within study strata becomes case only or control only
#   metadf <- data %>% 
#     dplyr::select(outcome, exposure, studyname, covars) %>% 
#     dplyr::filter(complete.cases(.))
#   cc <- data.frame(table(metadf$studyname, metadf$outcome)) %>% 
#     filter(Freq == 0)
#   metadf <- filter(metadf, !studyname %in% unique(cc$Var1)) %>% 
#     mutate(studyname = factor(studyname))
#   
#   formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"), "+", exposure)
#   
#   r_counts <- getCounts_outc(data = metadf, group = 'studyname')
#   saveRDS(r_counts, file = paste0('./meta/r_counts_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
#   r_glm <- glm_by_group_df(metadf, 'studyname', exposure, formula_txt = formula_txt)
#   saveRDS(r_glm, file = paste0('./meta/r_glm_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
#   #meta_fp(r_counts, r_glm, title = title) # find this function in the rmarkdown file NSAID_PostHarmonization.rmd
# }
# 
# # covars <- c("age_ref_imp", "sex", paste0(rep("PC", 10), seq(1,10)), "smoke")
# # covars <- c("age_ref_imp", "sex")
# # covars <- c("age_ref_imp", "sex", "smoke")
# # covars <- c("age_ref_imp", "sex")
# 
# covars <- c("age_ref_imp", "sex")
# getMeta(df, 'asp_ref', covars = covars)
# getMeta(df, 'aspirin', covars = covars)
# getMeta(df, 'nsaids', covars = covars)
# 
# 
# covars <- c("age_ref_imp", "sex", "smoke")
# getMeta(df, 'asp_ref', covars = covars)





#-----------------------------------------------------------------------------#
# Extract dosages from BinaryDosage files ------
#-----------------------------------------------------------------------------#

# write simple wrapper that takes SNP IDs vector as input, spits out index location
# do it by chromosome to avoid reading in RDS files over and over

write_binarydosage_vcfid_filename <- function() {
  return(paste0("files/GetSNPValues_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, "_vcfid"))
}

write_binarydosage_index_filename <- function(chr) {
  return(paste0("files/GetSNPValues_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, "_chr", chr))
}


get_binarydosage_index <- function(snplist, chr) {
  tmp_chr <- readRDS(paste0("/home/rak/data/BinaryDosage_InfoFile/FIGI_chr", chr, ".rds"))$SNPs %>% 
    dplyr::mutate(ID = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
  out <- as.integer(which(tmp_chr$ID %in% snplist))
  saveRDS(out, file = paste0(write_binarydosage_index_filename(chr), ".rds"), version = 2)
}



#-----------------------------------------------------------------------------#
# Plot allele frequencies for top hits ------
#-----------------------------------------------------------------------------#

# you should do it by study_gxe, so need to merge covariate file with G dosages
# generalize later, let's get it working first

# covariate_file <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_aspirin_sex_age_pc3_studygxe_66485_GLM.rds")
# dosages <- data.frame(readRDS("files/GetSNPValues_aspirin_age_ref_imp_sex_study_gxe_PC1-3_N_66485_chr5_out.rds")) %>% 
#   rownames_to_column(var = 'vcfid')
# posthoc_df <- inner_join(covariate_file, dosages, by = 'vcfid')
# 
# posthoc_df_maf <- posthoc_df %>% 
#   group_by(study_gxe) %>% 
#   summarise_at(vars(X5.40252294, X5.40273441), function(x) 0.5 - abs( (sum(x) / (2*nrow(.))) - 0.5))
# 
# ggplot(posthoc_df_maf) +
#   geom_point(aes(y = study_gxe, x = X5.40252294))
# 
# 
# ggplot(posthoc_df_maf) +
#   geom_point(aes(y = study_gxe, x = X5.40252294)) + 
#   theme_bw() + 
#   xlim(0,0.5)
#   # theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# aaf_X5.40252294 <- sum(posthoc_df$X5.40252294) / (2*nrow(posthoc_df))
# 
# 
# reate_plotMAF_study_gxe <- function()





#-----------------------------------------------------------------------------#
# Meta analysis of GxE ------
#-----------------------------------------------------------------------------#
# run GLM by study_gxe
# just use string formula_txt as input
# make sure only variables that you need indicators for are FACTORS (otherwise just continuous)
getGLM_byGroup_gxe <- function(data, outcome, exposure, genos, covars) {
  
  # data prep
  tmp <- create_data_4_meta(data, !! enquo(outcome), !! enquo(exposure)) %>% 
    group_by(study_gxe)
  
  results_beta <- do(tmp, tidy(glm(as.formula(paste0(quo_name(enquo(outcome)), " ~ ", quo_name(enquo(exposure)), "*", quo_name(enquo(genos)), "+", paste(covars, collapse = "+"))), data = . , family = 'binomial')))
  results_ci   <- do(tmp, confint_tidy(glm(as.formula(paste0(quo_name(enquo(outcome)), " ~ ", quo_name(enquo(exposure)), "*", quo_name(enquo(genos)), "+", paste(covars, collapse = "+"))), data = . , family = 'binomial')))
  
  results <- bind_cols(results_beta, results_ci) %>%
    ungroup %>%
    dplyr::filter(grepl(paste0(quo_name(enquo(exposure)), ":"), term)) %>%
    dplyr::select(-study_gxe1)
  results
  
}


# take outputs from above to perform meta analysis and plot forest plots
# just insert title yourself, make sure to include model
run_meta_analysis_create_forestplot_gxe <- function(counts_df, results_df, title) {
  
  # run meta analysis (summary OR + ConfInt)
  tmp_meta <- meta.summaries(results_df$estimate, results_df$std.error, method = 'random') # returns list object
  tmp_meta_b_ci <- c(tmp_meta[[3]], tmp_meta[[3]]-(tmp_meta[[4]]*1.96), tmp_meta[[3]]+(tmp_meta[[4]]*1.96))
  
  # format forestplot text (counts, study specific ORs)
  tmp_fp <- dplyr::select(results_df, estimate, conf.low, conf.high)
  tmp_fp <- rbind(rep(NA, 3), tmp_fp, rep(NA, 3), tmp_meta_b_ci) # add meta summaries, need conf.int for forestplot
  tmp_label_summary <- c("Summary", sum(counts_df$N), sum(counts_df$Case), round(tmp_meta[[3]], 2), round(tmp_meta[[4]], 3), round(tmp_meta[[5]][2], 4))
  
  tmp_label <- inner_join(counts_df, results_df, by = 'study_gxe') %>%
    mutate(Study = as.character(study_gxe),
           Pvalue = round(p.value, 4),
           SE = round(std.error, 3),
           beta = round(estimate, 2)) %>%
    dplyr::select(Study, N, Case, beta, SE, Pvalue)
  
  
  tmp_label2 <- rbind(colnames(tmp_label), tmp_label, rep(NA, 7), tmp_label_summary)
  rmeta::forestplot(label=tmp_label2,
                    as.numeric(tmp_fp$estimate), as.numeric(tmp_fp$conf.low), as.numeric(tmp_fp$conf.high),
                    is.summary = c(T, rep(F, nrow(tmp_label2)-2), T),
                    col=meta.colors(box="royalblue",line="darkblue",zero='gray0', summary="royalblue"),
                    clip=c(-1.5,1.5),
                    xlog=F,
                    xlab='GxE Beta (CI)',
                    zero=0, title(main = title, line = 0))
  mtext(paste("het.pval=", formatC(tmp_meta$het[3], format='e', digits=2)), side = 1, line = 1)
}

# x <- getCounts_byOutcome(posthoc_df, outcome, aspirin)
# y <- getGLM_byGroup_gxe(posthoc_df, outcome, aspirin, X5.40252294, c("age_ref_imp", "sex", "PC1", "PC2", "PC3"))
# z <- run_meta_analysis_create_forestplot_gxe(x, y, "testing")




#-----------------------------------------------------------------------------#
# calculate lambdas ------
# quick function to calculate lambdas and lambda1000
# make sure you use the appropriate DF for the result
# variables cases/controls/cases1000/controls1000 
#-----------------------------------------------------------------------------#
# getlambda <- function(pvals) {
#   chisq <- qchisq(1-pvals, 1)
#   lambda <- round(median(chisq)/qchisq(0.5,1),4)
# }
# 
# getlambda2df <- function(pvals) {
#   chisq <- qchisq(1-pvals, 2)
#   lambda <- round(median(chisq)/qchisq(0.5,2),4)
# }
# 
# getlambda3df <- function(pvals) {
#   chisq <- qchisq(1-pvals, 3)
#   lambda <- round(median(chisq)/qchisq(0.5,3),4)
# }
# 
# # takes gxescan results data + lambda calculated above
# getlambda1000 <- function(data, lambda) {
#   cases <- unique(data[, 'Cases'])
#   controls <- unique(data[, 'Subjects']) - unique(data[, 'Cases'])
#   total <- cases + controls
#   cases1000 <- (cases/total) * 1000
#   controls1000 <- (controls/total) * 1000
#   lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
# }




#-----------------------------------------------------------------------------#
# Results QQ and Manhattan Plots ------
# uses lambda functions above
# needs script called "~/Dropbox/FIGI/code/Functions/Rscript_create_ecf.R"
#-----------------------------------------------------------------------------#
# just add a variable called 'P'
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}

calculate_pval_miami <- function(data, statistic1, statistic2, df1, df2) {
  data$P1 <- pchisq(data[,statistic1], df = df1, lower.tail = F)
  data$P2 <- pchisq(data[,statistic2], df = df2, lower.tail = F)
  data
}



# dumb little helper functions for plot titles and file names
# NOTE SAVE LOCATION
# all arguments should be quote (Std evaluation...)
write_plot_filename <- function(data, statistic) {
  return(paste0("figures/QQ_Plot_", data, "_", statistic, "_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, ".png"))
}

write_easystrata_filename <- function(data, statistic) {
  return(paste0(paste("/media/work/tmp/EasyStrata", data, statistic, global_E, paste0(global_covs, collapse = "_"), "N", global_N, sep = "_"), ".txt"))
}

write_easystrata_filename_ecf <- function(data, statistic) {
  return(paste0(paste("files/EasyStrata", data, statistic, global_E, paste0(global_covs, collapse = "_"), "N", global_N, sep = "_"), ".ecf"))
}

write_plot_title <- function(statistic) {
  gxescan_tests <- c(paste0("G Main Effects Results (N = ", global_N, ")\noutc ~ G+", paste0(global_covs, collapse = "+"),"+", global_E), 
                     paste0("GxE Results (N = ", global_N, ")\noutc ~ G*", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("2DF Results (N = ", global_N, ")\noutc ~ G+G*", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("G|E Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("Case-Only Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("Control-Only Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("3DF Results (N = ", global_N, ")\nchiSqG+chiSqGxE+chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGxE", "chiSq2df", "chiSqGE", "chiSqCase", "chiSqControl", "chiSq3df")
  return(gxescan_tests[statistic])
}

write_plot_filename("gxe", "chiSqGxE")
write_easystrata_filename("gxe", "chiSqGxE")
write_easystrata_filename_ecf("gxe", "chiSqGxE")
# write_plot_title("chiSqG")
# write_plot_title("chiSqGxE")
# write_plot_title("chiSq2df")
# write_plot_title("chiSqGE")
# write_plot_title("chiSqCase")
# write_plot_title("chiSqControl")


# qq plot function
create_qqplot <- function(data, statistic, df) {
  
  # testing
  # statistic = 'chiSqGxE'
  # data = gxe
  # df = 1
  
  # calculate p value
  # calculate lambda
  tmpdata <- calculate_pval(data, statistic, df)
  lambda <- round( (median(qchisq(1-tmpdata[,'P'], df)) / qchisq(0.5, df)), 4)
  
  # calculate lambda1000
  cases <- unique(tmpdata[, 'Cases'])
  controls <- unique(tmpdata[, 'Subjects']) - unique(tmpdata[, 'Cases'])
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  
  # plotting function
  png(write_plot_filename(deparse(substitute(data)), statistic), height = 720, width = 1280)
  qqman::qq(tmpdata[, 'P'], 
            xlab = "Expected -log10(p)", 
            ylab = "Observed -log10(p)",
            main = write_plot_title(statistic),
            cex.main = 1.6, 
            cex.axis = 1.3, 
            cex.lab = 1.3,
            cex.sub = 1.3,
            col = 'blue4') 
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3) # FYI ~~ adds spaces when using signif
  dev.off()
}


# Manhattan Plot
# work pending - modify annotion to match chr:pos:ref:alt instead of just chr:pos
create_manhattanplot <- function(data, statistic, df) {
  
  # testing
  # statistic = 'chiSqG'
  # data = gxe
  # df = 1
  
  # format data for easystrata (just be consistent)
  # remember 'calculate_pval' creates variable 'P'
  tmpdata <- calculate_pval(data, statistic, df)
  
  # need to write results to table. using same location for now (temporary directory in /media/work)
  data_easystrata <- tmpdata %>%
    mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>%
    filter(!(P > 0.05 & annot == 0)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(ID, CHR, BP, P)
  write.table(data_easystrata, file = write_easystrata_filename(deparse(substitute(data)), statistic), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(getwd(), "/figures")
  ecf2 <- "ID;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- write_easystrata_filename(deparse(substitute(data)), statistic)
  ecf_file_name <- write_easystrata_filename_ecf(deparse(substitute(data)), statistic)
  source("/home/rak/Dropbox/FIGI/Code/Functions/Rscript_create_ecf.R", local = T) # edit script as necessary
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}




# MIAMI Plot
# work pending - modify annotion to match chr:pos:ref:alt instead of just chr:pos
create_miamiplot <- function(data, statistic1, statistic2, df1, df2) {
  
  # format data for easystrata (just be consistent)
  # remember 'calculate_pval' creates variable 'P'
  tmpdata <- calculate_pval_miami(data, statistic1, statistic2, df1, df2)
  
  # need to write results to table. using same location for now (temporary directory in /media/work)
  data_easystrata <- tmpdata %>% 
    mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
    filter(!((P1 > 0.05 | P2 > 0.05) & annot == 0)) %>% 
    dplyr::rename(CHR = Chromosome, 
                  BP = Location) %>% 
    dplyr::select(ID, CHR, BP, P1, P2)
  write.table(data_easystrata, file = write_easystrata_filename(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_")), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(getwd(), "/figures")
  ecf2 <- "ID;CHR;BP;P1;P2"
  ecf3 <- "character;numeric;numeric;numeric;numeric"
  ecf4 <- write_easystrata_filename(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_"))
  ecf_file_name <- write_easystrata_filename_ecf(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_"))
  source("/home/rak/Dropbox/FIGI/Code/Functions/Rscript_create_ecf_miami.R", local = T) # edit script as necessary
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}




#-----------------------------------------------------------------------------#
# Weighted Hypothesis Testing for 2-step methods ------
# kooperberg, murcray, edge
#-----------------------------------------------------------------------------#

# variables required
# step1 statistic
# size of initial bin 
# number of SNPs
# overall alpha level = 0.05


# function to format step 1 (arrange by p value) --- assume data is the gxe object from GxEScanR
# output: data.table (because original code used data.tables) with bin number, bin logp threshold, normalized X axis info for plotting
format_2step_data <- function(data, step1_statistic, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(data, statistic, df=1) {
    data.table(data)[
      , step1p := pchisq(data[,statistic], df = df, lower.tail = F)
      ][
        , step2p := pchisq(data[,'chiSqGxE'],  df = 1, lower.tail = F)
        ][
          , y := -log10(step2p)
          ][
            order(step1p)
            ][
              , MapInfo := Location
              ]
  }
  
  if(step1_statistic == 'chiSqEDGE') {
    pv <- create_pval_info(data, step1_statistic, df = 2)
  } else {
    pv <- create_pval_info(data, step1_statistic, df = 1)
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




# file name =D
write_weighted_test_plot_title <- function(statistic) {
  gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", global_N, ")\noutc ~ G+", paste0(global_covs, collapse = "+"),"+", global_E), 
                     paste0("G|E 2-step Procedure Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("EDGE 2-step Procedure Results (N = ", global_N, ")\nchiSqG + chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
  return(gxescan_tests[statistic])
}

write_weighted_plot_filename <- function(data, statistic) {
  return(paste0("figures/TwoStep_WeightedHypothesis_", data, "_", statistic, "_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, ".png"))
}



# create weighted hypothesis plot
# first step is creating a list object based on data.table created by the 'format_2step_data' function
# second step is actually plotting. remember that it only plots first 15 bins for legibility 
create_2step_weighted_plot <- function(data, sizeBin0, alpha, binsToPlot, statistic) {
  
  # m = nrow(koop_filter)
  # sizeBin0 = 5
  # alpha = 0.05
  # binsToPlot = 15
  # data = test
  
  # bin information
  m = nrow(data)
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
    t <- data[J(i)]
    t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
    t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
    glist[[i]]<-t
    rm(t)
  }
  
  # trying to understand code above (just arrange BP by bins looks like)
  # (Scale mapinfo for each Bin to range between 0.1-0.9 for neatness, and add a unit increase for successive Bin)
  # x <- pv[1:5, MapInfo]
  # normalized = (x-min(x))/(max(x)-min(x))
  # normalized_scaled = 0.8 * normalized + 0.1
  # x;normalized;normalized_scaled
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(write_weighted_plot_filename(deparse(substitute(data)), statistic), width = 1280, height = 720)
  color <- rep(c("blue","olivedrab4"),100)
  plot(glist[[1]][,x], glist[[1]][,y], 
       col = "blue", 
       xlab="Bin # for step1 p-values", 
       ylab="-log10(step2 p-values)", 
       xlim=c(0,binsToPlot), 
       ylim=c(0,min.p), 
       axes=F, pch=19, cex=0.5)
  lines(glist[[1]][,x], glist[[1]][,ref],
        col = "black",lwd=2)
  
  # the rest of the points ...  =|
  # (adding to current plot..)
  for(i in 2:binsToPlot){
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
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 0.8)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=0.8)
  title(main = write_weighted_test_plot_title(statistic), sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.5, cex.sub = 1.2)
  
  dev.off()
  
}








