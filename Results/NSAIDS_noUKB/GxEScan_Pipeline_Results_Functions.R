#=============================================================================#
# Useful Functions for analysis pipeline / results
#
# mostly uses NSE (tidyverse)
#=============================================================================#


#-----------------------------------------------------------------------------#
# output data frame of
# NA/No/Yes counts (%) by outcome (case/control) and study/platform
# should format exposure properly
#-----------------------------------------------------------------------------#

getCounts <- function(data, group, exposure, outcome) {
  
  # quote variables
  group_quo <- enquo(group)
  exposure_quo <- enquo(exposure)
  outcome_quo <- enquo(outcome)
  
  # formula for dcast function
  dcastFormula = as.formula(paste0(quo_name(group_quo), "~", quo_name(outcome_quo), "+", quo_name(exposure_quo)))
  
  # get totals first
  countsTotal <- data %>% 
    dplyr::group_by(!! outcome_quo) %>% 
    dplyr::select(!! exposure_quo, !! outcome_quo) %>%
    dplyr::count(!! exposure_quo) %>%
    dplyr::mutate(tot = sum(n),
                  pct = round(100*(n/tot), 1),
                  n_pct = paste0(n, " (", pct, ")"),
                  !! group_quo := "Total") %>% dplyr::select(-n, -pct) %>% ungroup() %>%
    reshape2::dcast(dcastFormula, fill = "", value.var = 'n_pct')
  
  # get counts by group
  countsGroup <- data %>%
    dplyr::group_by(!! group_quo, !! outcome_quo) %>%
    dplyr::select(!! group_quo, !! exposure_quo, !! outcome_quo) %>%
    dplyr::count(!! exposure_quo) %>%
    dplyr::mutate(tot = sum(n),
                  pct = round(100*(n/tot), 1),
                  n_pct = paste0(n, " (", pct, ")")) %>% dplyr::select(-n, -pct) %>% ungroup() %>%
    reshape2::dcast(dcastFormula, fill = "", value.var = 'n_pct')
  
  results <- rbind(countsGroup, countsTotal)
}

# test <- getCounts(gxe_set, group = study_gxe, exposure = folate_totqc2, outcome = outc)



#-----------------------------------------------------------------------------#
# Barplot for exposure variable (ggplot2)
# facet_grid by 'group' variable (usually study_gxe)
# number of columns is set to 9.. 
#-----------------------------------------------------------------------------#
getPlots <- function(data, group, exposure, outcome) {
  
  # quote variables
  group_quo <- enquo(group)
  exposure_quo <- enquo(exposure)
  outcome_quo <- enquo(outcome)
  
  # create a temporary grouped dataset
  cov <- data %>%
    dplyr::select(!! outcome_quo, !! group_quo, !! exposure_quo) %>%
    group_by(!! group_quo)
  
  # color code appropriately
  # quo_text(exposure_quo)
  fac_length <- levels(cov[, quo_text(exposure_quo)])
  
  # fac_length <- length(levels(cov[, quo_text(exposure_quo)]))
  # 
  # if (fac_length == 3) {
  #   colScale <- scale_fill_manual(values = c("black", "cyan", "red"))
  # } else if (fac_length == 5) {
  #   colScale <- scale_fill_manual(values = c("black", "#730BFE", "#00FFFF", "#70FF00", "#FF0006"))
  # } else {
  #   colScale <- scale_fill_manual()
  # }
  # 
  # 
  # p <- ggplot(data = cov, aes(x = !! outcome_quo)) +
  #   geom_bar(aes(fill = !! exposure_quo), position = 'fill') +
  #   theme_bw() +
  #   theme(
  #     panel.grid.major.x = element_blank(),
  #     panel.grid.minor.x = element_blank()) +
  #     colScale
    # scale_fill_manual(values = c( "black", "cyan", "red"))
    # facet_wrap(group_quo, ncol = 9)
    # facet_wrap(group_quo)
}
# 
# test <- getPlots(data = gxe_set, exposure = exposure_fac, group = study_gxe, outcome = outc_fac)
# test






#-----------------------------------------------------------------------------#
# Functions for meta analyses + forestplot using package 'rmeta'
# again, make sure outc is coded "Control", "Case"...

getCounts_outc <- function(data, group, exposure, outcome, ...) {
	
	group_quo <- enquo(group)
	exposure_quo <- enquo(exposure)
	outcome_quo <- enquo(outcome)
	covars_quo <- enexprs(...)
	
	tmp_counts <- data %>% 
		dplyr::select(!! outcome_quo, !! exposure_quo, !!group_quo, !!! covars_quo) %>% 
		filter(complete.cases(.)) %>% # probably doing it in another step, but why not do it again
		group_by(!! group_quo)
	
	tmp_counts <- tmp_counts %>% 
		dplyr::group_by(!! outcome_quo, !! group_quo) %>% 
		dplyr::summarize(count = n()) %>% 
		tidyr::spread(key = !! outcome_quo, value = count) %>% 
		dplyr::mutate(N = Control + Case)
	tmp_counts
}

# x <- getCounts_outc(df, study_gxe, outc)


# ... means covariates (just list them in the function call)
glm_by_group_df <- function(data, group, exposure, outcome, ...) {
	
	group_quo <- enquo(group)
	exposure_quo <- enquo(exposure)
	outcome_quo <- enquo(outcome)
	covars_quo <- enexprs(...) # it's like enquos, but it doesn't capture the environment...
	
	glm_formula = as.formula(paste(quo_name(outcome_quo), "~", quo_name(exposure_quo), "+", paste(covars_quo, collapse = " + ")))

	tmp_group <- data %>% 
		dplyr::select(!! outcome_quo, !! exposure_quo, !!group_quo, !!! covars_quo) %>% 
		filter(complete.cases(.)) %>% # probably doing it in another step, but why not do it again
		group_by(!! group_quo)
	
	tmp_results_beta <- do(tmp_group, tidy(glm(glm_formula, data = . , family = 'binomial')))
	tmp_results_ci   <- do(tmp_group, confint_tidy(glm(glm_formula, data = . , family = 'binomial')))
	tmp_results <- bind_cols(tmp_results_beta, tmp_results_ci) %>% 
		ungroup %>%
		dplyr::filter(term == quo_name(exposure_quo)) %>% 
		dplyr::select(-c(paste0(quo_name(group_quo), 1)))
	# tmp_results
}

# y <- glm_by_group_df(data = df, group = study_gxe, exposure = asp_ref, outcome = outc_01, age_ref_imp, sex)

meta_fp <- function(counts_df, results_df, group, title) {
	
	group_quo <- enexpr(group)
	
	tmp_meta <- meta.summaries(results_df$estimate, results_df$std.error, method = 'random')
	tmp_meta_b_ci <- c(tmp_meta[[3]], tmp_meta[[3]]-(tmp_meta[[4]]*1.96), tmp_meta[[3]]+(tmp_meta[[4]]*1.96))
	tmp_fp <- dplyr::select(results_df, estimate, conf.low, conf.high)
	tmp_fp <- rbind(rep(NA, 3), tmp_fp, rep(NA, 3), tmp_meta_b_ci)
	tmp_label_summary <- c("Summary", sum(counts_df$N), sum(counts_df$Case), round(tmp_meta[[3]], 2), round(tmp_meta[[4]], 3), round(exp(tmp_meta[[3]]), 2), round(tmp_meta[[5]][2], 4))
	tmp_label <- inner_join(counts_df, results_df, by = quo_name(group_quo)) %>% 
		mutate(Study = as.character(!! group_quo),
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
	mtext(paste(results_df$term[1], '\n', "het.pval=", formatC(tmp_meta$het[3], format='e', digits=2)), side = 1, line = 1)
}

# z <- meta_fp(x, y, group = study_gxe, title = "wtf")



# getMeta <- function(data, group, exposure, outcome, covars) {
# 	
# 
# 	# (this bit I used to do to ensure that there are no case only control only studies
# 	# just make sure prior to running function)
# 	# metadf <- data %>% 
# 	# 	dplyr::select(outcome, exposure, studyname, covars) %>% 
# 	# 	dplyr::filter(complete.cases(.)) # is this necessary since it's based on GLM... 
# 	# 
# 	# cc <- data.frame(table(metadf[, group], metadf[, outcome])) %>%
# 	# 	filter(Freq == 0)
# 	# 
# 	# metadf <- filter(metadf, !studyname %in% unique(cc$Var1)) %>%
# 	# 	mutate(studyname = factor(studyname))
# 	
# 	formula_txt <- paste0(outcome, "~", paste(covars, collapse = "+"), "+", exposure)
# 	
# 	# r_counts <- getCounts_outc(data = metadf, group = 'studyname')
# 	# saveRDS(r_counts, file = paste0('./meta/r_counts_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
# 	# r_glm <- glm_by_group_df(metadf, 'studyname', exposure, formula_txt = formula_txt)
# 	# saveRDS(r_glm, file = paste0('./meta/r_glm_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
# 	# meta_fp(r_counts, r_glm, title = title) # find this function in the rmarkdown file NSAID_PostHarmonization.rmd
# }
