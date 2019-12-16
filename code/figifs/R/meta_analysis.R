#=============================================================================#
# Meta Analyses
#
# Functions to organize data, conduct meta-analyses, and create plots
#=============================================================================#


#' remove_low_count_cells
#'
#' With data subsets, cell counts of outcome, exposure, and study_gxe can become low and cause regression errors. If exposure is categorical, 3-way tabulate exposure/outcome/study_gxe, remove studies with cell counts <= 10. If numeric, 2-way tabulate outcome/study_gxe, remove studies with cell counts == 0.
#'
#' @param dat Input data
#' @param is_categorical Logical, categorical TRUE or FALSE
#'
#' @return Data
#' @export
#'
#' @examples remove_low_count_cells(data, is_categorical = TRUE)
remove_low_count_cells <- function(dat, is_categorical) {
  if(is_categorical == T) {
    drops <- data.frame(table(dat$outcome, dat[,params$exposure], dat$study_gxe)) %>%
      filter(Freq <= 10)
    dat <- filter(dat, !study_gxe %in% unique(drops$Var3))
    return(dat)
  } else {
    drops <- data.frame(table(dat$outcome, dat$study_gxe)) %>%
      filter(Freq == 0)
    dat <- filter(dat, !study_gxe %in% unique(drops$Var3))
    return(dat)
  }
}


#' get_counts_outcome_by_group
#'
#' Tabulate outcome by study/platform to inclusion in forest plots. Input data should not include studies with zero counts when tabulating exposure/outcome/study.
#'
#' @param dat Input data
#' @param outcome Outcome variable
#' @param group Group variable e.g. study_gxe
#'
#' @return Tibble with case/control counts by group.
#' @export
#'
#' @examples get_counts_by_outcome(gxe, outcome, aspirin, study_gxe)
get_counts_outcome_by_group <- function(dat, outcome, group) {

  # get counts for rmeta function
  dat <- dat %>%
    dplyr::group_by( .data[[outcome]], .data[[group]] ) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    tidyr::spread(key = {{ outcome }}, value = count) %>%
    dplyr::rename(Control = `0`, Case = `1`) %>%
    dplyr::mutate(N = Control + Case)
  return(dat)

}


#' get_estimates_e_by_group
#'
#' Fit regression for exposure main effect by group. Input data should not include studies with zero counts when tabulating exposure/outcome/study.
#'
#' @param dat Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable e.g. study_gxe
#' @param ... Adjustment covariates
#'
#' @return Tibble with exposure GLM estimates/se/stats/pval/95%CI by group
#' @export
#'
#' @examples run_glm_e_by_group(gxe, outcome, asp_ref, study_gxe, age_ref_imp, sex, PC1, PC2, PC3)
get_estimates_e_by_group <- function(dat, outcome, exposure, group, ...) {

  # model covariates as vector
  covariates <- c(...)

  # group data, perform grouped glm
  dat <- dat %>%
    group_by( .data[[group]] )

  glm_formula <- reformulate(termlabels = c( {{ exposure }}, covariates) , response = {{ outcome }} )
  results_beta <- dplyr::do(dat, broom::tidy(glm(glm_formula, data = . , family = 'binomial')))
  results_ci   <- dplyr::do(dat, broom::confint_tidy(glm(glm_formula, data = . , family = 'binomial')))

  results <- dplyr::bind_cols(results_beta, results_ci) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term == exposure)
  return(results)

}



#' meta_analysis_wrapper
#'
#' Wrapper function to run meta-analysis, create forest and funnel plots. Uses package 'meta'.
#'
#' @param dat Data frame containing counts and regression estimates BY GROUP
#' @param forest_plot_title Plot title
#' @param filename_suffix Filename suffix following exposure (usually covariates)
#' @param forest_height png height
#' @param forest_width png width
#' @param funnel_height png height
#' @param funnel_width png width
#'
#' @return Outputs png files for forest and funnel plots
#' @export
#'
#' @examples meta_analysis_wrapper(gxe_meta, "title", "age_ref_imp_sex_study_gxe", 13, 8.5, 8, 8.5)
meta_analysis_wrapper <- function(dat, forest_plot_title, filename_suffix, forest_height, forest_width, funnel_height, funnel_width) {

  results_meta <- meta::metagen(estimate,
                                std.error,
                                data=dat,
                                studlab=paste(study_gxe),
                                comb.fixed = FALSE,
                                comb.random = TRUE,
                                method.tau = "SJ",
                                hakn = TRUE,
                                prediction=TRUE,
                                sm="OR",
                                byvar = study_design)

  png(paste0("/media/work/tmp_images/meta_analysis_", params$exposure,  "_", filename_suffix, ".png"), height = forest_height, width = forest_width, units = 'in', res = 150)
  meta::forest(results_meta,
               layout = "JAMA",
               # text.predict = "95% CI",
               # col.predict = "black",
               leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
               digits.addcols=0,
               study.results=T,
               prediction = F,
               col.random = 'red')
  grid.text(forest_plot_title, 0.5, .98, gp=gpar(cex=2))
  dev.off()

  png(paste0("/media/work/tmp_images/funnel_plot_", params$exposure,  "_", filename_suffix, ".png"), height = funnel_height, width = funnel_width, units = 'in', res = 150)
  meta::funnel(results_meta, sm="OR", studlab = T, pos = 4, col.random = 'red')
  dev.off()

}








#' get_estimates_gxe_by_group
#'
#' Function to get GLM estimates by group (e.g. study_gxe). Make sure data is properly formatted yeah (maybe add info on what that means exactly..). This is an alternate version that returns the gxe term in the summary table
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable (typically study_g)
#' @param dosage SNP for interaction with E
#' @param ... Adjustment Covariates
#'
#' @return Tibble with GLM estimates/se/stats/pval/95%CI by group
#' @export
#'
#' @examples run_glm_gxe_by_group(gxe, outcome, asp_ref, study_gxe, X5.40252294, age_ref_imp, sex, PC1, PC2, PC3)
get_estimates_gxe_by_group <- function(df, outcome, exposure, group, dosage, ...) {

  outcome <- enquo(outcome)
  exposure <- enquo(exposure)
  group <- enquo(group)
  dosage <- enquo(dosage)
  covariates <- enquos(...)

  # data prep
  df <- format_data_meta_analysis(df, !! outcome, !! exposure, !! group) %>%
    dplyr::group_by(!! group)

  glm_formula <- reformulate(termlabels = c(paste0(quo_name(exposure),"*", quo_name(dosage)),
                                            as.vector(unlist(lapply(covariates, quo_name)))),
                             response = quo_name(outcome))

  results_beta <- dplyr::do(df, broom::tidy(glm(glm_formula, data = . , family = 'binomial')))
  results_ci   <- dplyr::do(df, broom::confint_tidy(glm(glm_formula, data = . , family = 'binomial')))

  results <- dplyr::bind_cols(results_beta, results_ci) %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(paste0(quo_name(exposure), ":"), term))
  return(results)

}




