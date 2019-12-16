#=============================================================================#
# FIGI postharmonization analysis
# 06/17/2019
#=============================================================================#

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
  

