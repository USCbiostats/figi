# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#-----------------------------------------------------------------------------#
# Descriptive Tables ------
#-----------------------------------------------------------------------------#

#' create_descriptive_table
#'
#' Convenience wrapper function around table1 package. Function is limited but works for simple cases like distribution of covariates by outcome. Code includes workaround for adding p values that the package author described in vignette.
#'
#' @param df Input data
#' @param outcome Outcome variable
#' @param exposure Exposure variable
#' @param covarAll String vector of covariates to include in table
#' @param covarCat String vector of categorical variables (factors)
#'
#' @return table1 table (use in Rmarkdown)
#' @export
#'
#' @examples create_descriptive_table(gxe, outcome, aspirin, c('sex', 'age'), c('sex'))
create_descriptive_table <- function(dat, outcome, exposure, covariates, covarCat) {

  outcome <- rlang::enquo(outcome)
  exposure <- rlang::enquo(exposure)

  # create 3 level outcome factor variable for p value
  # categorical variables should be factors, continuous numeric
  dat <- dat %>%
    dplyr::filter(!is.na(!! exposure)) %>%
    dplyr::mutate(table1_outcome = factor(!! outcome, levels = c(0,1,2), labels=c("Control", "Case", "P-value"))) %>%
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

  # insert p values. hacky, requires an outcome factor with third level for p value.
  # Also need to specify values by name (table1_outcome) in the 'rndr' function
  rndr <- function(x, name, ...) {
    if (length(x) == 0) {
      y <- dat[[name]]
      s <- rep("", length(render.default(x=y, name=name, ...)))
      if (is.numeric(y)) {
        p <- t.test(y ~ dat$table1_outcome)$p.value
      } else {
        p <- chisq.test(table(y, droplevels(dat$table1_outcome)))$p.value
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

  table1(as.formula(paste0("~ ", get_expr(exposure), "+", paste(covariates, collapse = "+"), "| table1_outcome")),
         data=dat,
         render.continuous=my.render.cont,
         render.categorical=my.render.cat,
         render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)
}






