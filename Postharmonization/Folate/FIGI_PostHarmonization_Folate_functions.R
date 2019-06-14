# helpful functions

#-------------------------------------------------------#
# Functions
#-------------------------------------------------------#

## Create table of case/control counts
## needs outc variable (define in data step). Group is typically studyname
counts_outc <- function(data, group) {
  data %>% 
    dplyr::group_by_('outc', group) %>% 
    dplyr::summarize(count = n()) %>% 
    tidyr::spread(key = outc, value = count) %>%
    dplyr::mutate(N = Control + Case)
}

## Conduct GLM by studyname
## outputs a table with components necessary for meta analyses
glm_by_group_df <- function(data, model_formula) {
  dd <- data %>%
    group_by(studyname)
  tmp_results_beta <- do(dd, tidy(glm(model_formula, data = . , family = 'binomial')))
  tmp_results_ci   <- do(dd, confint_tidy(glm(model_formula, data = . , family = 'binomial'))) %>%
    filter(!is.na(conf.low))
  tmp_results <- bind_cols(tmp_results_beta, tmp_results_ci) %>% ungroup
  tmp_results
}

# glm_by_group_df <- function(data, group, exposure, formula_txt) {
#   formula.as.text <- paste0(formula_txt, "+", exposure)
#   tmp_group <- data %>%
#     group_by(studyname)
#   tmp_results_beta <- do(tmp_group, tidy(glm(as.formula(formula.as.text), data = . , family = 'binomial')))
#   tmp_results_ci   <- do(tmp_group, confint_tidy(glm(as.formula(formula.as.text), data = . , family = 'binomial'))) %>%
#     filter(!is.na(conf.low))
#   # i'm also creating variables you'd want to explore heterogeneity by - in this case, study design
#   tmp_results <- bind_cols(tmp_results_beta, tmp_results_ci) %>%
#     ungroup %>% 
#     dplyr::filter(term == exposure) 
#   tmp_results
# }


## round with zeros
round_with_zeros <- function(.x, digits = 2) {
  format(round(.x, digits = digits), nsmall = digits)
}


## Wrapper to run malcolm meta-analysis (tidymeta)
## BY STUDYNAME + 1 STUDY ATTRIBUTE (e.g. study design)
exposure <- "folate_tot"
covars <- c("age_ref", "sex.n")
group <- "studydesign"
data <- df

ma_wrapper <- function(exposure, covars, group, data) {
  
  model_form = as.formula(paste0("outcome ~ ", exposure, "+", paste0(covars, collapse = "+")))
  
  # complete cases
  dd <- data %>% 
    dplyr::select(c('outcome', exposure, covars, group, 'studyname', 'outc')) %>% 
    filter(complete.cases(.))
  
  # counts tables (for forest plots...)
  r_counts <- counts_outc(data = dd, group = 'studyname')
  
  r_counts_all <- r_counts %>% 
    summarise_at(c('Control', 'Case', 'N'), funs(sum)) %>% 
    mutate(studyname = "Overall")
 
  r_counts_design <- inner_join(r_counts, unique(dd[, c('studyname', group)]), by = 'studyname') %>% 
    group_by_(group) %>% 
    summarise_at(c('Control', 'Case', 'N'), funs(sum)) %>% ungroup() %>% 
    mutate(studyname = paste0("Subgroup: ", !!sym(group))) # have to call it 'studyname' to match other stuff for Forest Plot.. 
  
  r_counts_use <- bind_rows(r_counts, r_counts_all, r_counts_design) %>%  dplyr::select(-group)
  
  ## Run GLM by studyname + conduct random effects meta analysis
  r_glm_meta <- glm_by_group_df(data = dd, model_formula = model_form) %>% 
    filter(term == exposure) %>% 
    inner_join(., unique(dd[, c('studyname', group)]), by = 'studyname') %>% 
    group_by_(group) %>% 
    meta_analysis(yi = estimate, sei = std.error, slab = studyname, exponentiate = TRUE) %>% ungroup() %>% 
    mutate(studyname = factor(study, levels = rev(study)), # Fyi the 'study' var is created in meta_analysis
           orci = paste0(round_with_zeros(estimate), " (",
                         round_with_zeros(conf.low), ", ",
                         round_with_zeros(conf.high), ")")) %>% 
    full_join(., r_counts_use, by = c("study" = "studyname")) %>% 
    mutate(counts = paste0(Case, " / ", Control))
}

x <- ma_wrapper(exposure = "folate_tot_400", covars <- c("age_ref", "sex.n"), group <- "studydesign", data = df)

## Wrapper to run malcolm meta-analysis + FP
fp_wrapper <- function(ma_wrapper_output, exposure, covars, formula_txt, plot_title, plot_subtitle, aa=2.5, bb=5) {
  
  # clean this up when you get a chance. 
  breaks_subset <- r_glm_meta[which(r_glm_meta$meta != "NULL"), 'study', drop = T]
  breaks_I2 <- filter(r_glm_meta, meta != "NULL") %>% 
    dplyr::select(meta) %>% 
    apply(., 1, function(x) x[[1]]$I2)
  breaks_QEp <- filter(r_glm_meta, meta != "NULL") %>% 
    dplyr::select(meta) %>% 
    apply(., 1, function(x) x[[1]]$QEp)
  
  breaks_h <- as.factor(paste0(breaks_subset, "\n", "(I2=", round(breaks_I2, 2), "%, p=", round(breaks_QEp, 3), ")"))
  breaks <- factor(breaks_h, levels = rev(levels(breaks_h)))
  faces <- as.expression(breaks)
  vline_or_data <- data.frame(z = 1, variable = "estimate")
  
  # dirty (sorry)
  breaks_df <- data.frame(breaks, design = c("Cohort", "CaseControl", "Summary"))
  r_glm_meta_h <- inner_join(r_glm_meta, breaks_df, by = 'design') %>% 
    dplyr::select(-design) %>% mutate(design = breaks)
  
  ## need separate components to plot forest + text with facet grids!
  dat_ggplot_main <- r_glm_meta_h %>% 
    dplyr::select(study, design, type, weight, conf.low, conf.high, estimate, orci, counts) %>% 
    tidyr::gather(facetvar, value, estimate:counts) %>% 
    mutate(facetvar = factor(facetvar, levels = c("estimate", "orci", "counts")))
  
  dat_ggplot_orci <- dat_ggplot_main %>% 
    filter(facetvar == "estimate") %>% 
    mutate(estimate = as.numeric(value), 
           studyname = factor(study, levels = rev(r_glm_meta_h$study)))
  
  dat_ggplot_text <- dat_ggplot_main %>% 
    filter(facetvar != "estimate") %>% 
    mutate(xcoord = ifelse(facetvar == "orci", aa, bb))
  
  labels <- c(estimate = "Forest Plot", orci = "OR (95% CI)", counts = "Case/Control")
  breaks_studyname <- levels(factor(r_glm_meta_h$study, levels = rev(r_glm_meta_h$study)))
  faces <- as.expression(breaks_studyname)
  
  
  ggplot(dat_ggplot_main, y = studyname) +
    theme_minimal() +
    facet_grid(design ~ ., scales = 'free', space = 'free_y', labeller = labeller(variable = labels)) + 
    geom_vline(data = vline_or_data, aes(xintercept = 1), linetype = 'dashed') +
    geom_point(data = dat_ggplot_orci, aes(x = estimate, y = studyname, col = design, size = weight, shape = type), alpha = .75) +
    geom_errorbarh(data = dat_ggplot_orci, aes(x = estimate, y = study, xmin = conf.low, xmax = conf.high, col = design), height = 0, size = .75, alpha = .75) +
    scale_x_continuous(trans = "log",
                       breaks = c(.5, 1, 2.0)) +
    coord_cartesian(xlim=c(0.4,7)) + 
    geom_text(data = dat_ggplot_text, aes(x = xcoord, y = study, label = value), size = 3.5) + 
    scale_y_discrete(label = faces, breaks = breaks_studyname) +
    scale_shape_manual(values = c(15, 18)) + 
    theme(strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color = "grey85", fill = "grey85"), # grey
          #strip.background = element_rect(color = "black", fill = "white", size = .75), # white with black borders
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_rect(color = "black", fill = NA, size = 1), # grid around panels
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 11),
          legend.position = "none") +
    scale_size_continuous(range = c(2, 8)) + 
    ggtitle(plot_title, subtitle = plot_subtitle)
  
  ggsave(paste0("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_", exposure, "_tidymeta.png"), width = 8, height = 7, scale = 1.2)
  
}
