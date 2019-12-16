# old versions of functions 

#' format_data_meta_analysis
#'
#' Quick convenience function to format FIGI data for meta-analysis. Not exported.
#'
#' @section Steps:
#' 1) Remove samples where exposure is not available
#'
#' 2) remove case only and control only studies
#'
#' 3) check cell sizes, drop study if case or control N < 10 to avoid convergence issues
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable (typically study_g)
#'
#' @return data.frame for meta-analysis downstream steps
#' @export
#'
#' @examples format_data_meta_analysis(gxe, outcome, asp_ref, study_gxe)
format_data_meta_analysis <- function(df, outcome, exposure, group) {
  outcome <- enquo(outcome)
  exposure <- enquo(exposure)
  group <- enquo(group)
  
  # filter missing exposure
  df <- dplyr::filter(df, !is.na(!! exposure))
  
  tab_group_outcome <- data.frame(table(df[, quo_name(group)], df[, quo_name(outcome)]))
  
  # filter case only and control only studies
  drops_caseonly <- as.vector(unique(dplyr::filter(tab_group_outcome, Freq == 0)[ , 1]))
  df <- dplyr::filter(df, !(!! group %in% drops_caseonly))
  
  # check for cell sizes, drop if study_gxe vs outcome cell N < 10
  drops_cellsize <- as.vector(unique(dplyr::filter(tab_group_outcome, Freq < 10)[ , 1]))
  df <- dplyr::filter(df, !(!! group %in% drops_cellsize))
}




#' get_counts_outcome_by_group
#'
#' Function to create
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable (typically study_g)
#'
#' @return Tibble with group counts (case/control)
#' @export
#'
#' @examples get_counts_by_outcome(gxe, outcome, aspirin, study_gxe)
get_counts_outcome_by_group <- function(df, outcome, exposure, group) {
  
  outcome <- enquo(outcome)
  exposure <- enquo(exposure)
  group <- enquo(group)
  
  df <- format_data_meta_analysis(df, !! outcome, !! exposure, !! group)
  
  # get counts for rmeta function
  df <- df %>%
    dplyr::group_by(!! outcome, !! group) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    tidyr::spread(key = !! outcome, value = count) %>%
    dplyr::rename(Control = `0`, Case = `1`) %>%
    dplyr::mutate(N = Control + Case)
  return(df)
}








#' get_estimates_e_by_group
#'
#' Function to get GLM estimates by group (e.g. study_gxe). Make sure data is properly formatted yeah (maybe add info on what that means exactly..)
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param group Group variable (typically study_g)
#' @param ... Adjustment Covariates
#'
#' @return Tibble with GLM estimates/se/stats/pval/95%CI by group
#' @export
#'
#' @examples run_glm_e_by_group(gxe, outcome, asp_ref, study_gxe, age_ref_imp, sex, PC1, PC2, PC3)
get_estimates_e_by_group <- function(df, outcome, exposure, group, ...) {
  
  outcome <- enquo(outcome)
  exposure <- enquo(exposure)
  group <- enquo(group)
  covariates <- enquos(...)
  
  # calling the function with group name = study_gxe
  # complete case analysis
  df <- format_data_meta_analysis(df, !! outcome, !! exposure, !! group) %>%
    dplyr::select(!! outcome, !! exposure, !! group, !!! covariates) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::group_by(!! group)
  
  # as.vector(unlist(lapply(covariates, quo_name)))
  glm_formula <- reformulate(termlabels = c(quo_name(exposure), as.vector(unlist(lapply(covariates, quo_name)))),
                             response = quo_name(outcome))
  
  results_beta <- dplyr::do(df, broom::tidy(glm(glm_formula, data = . , family = 'binomial')))
  results_ci   <- dplyr::do(df, broom::confint_tidy(glm(glm_formula, data = . , family = 'binomial')))
  
  results <- dplyr::bind_cols(results_beta, results_ci) %>%
    dplyr::ungroup() %>%
    dplyr::filter(grepl(quo_name(exposure), term))
  # dplyr::select(-dplyr::contains(quo_name(group)))
  return(results)
  
}