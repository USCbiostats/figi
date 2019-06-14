#======================================================#
# Generate meta analysis, forest and funnel plots
#
# using malcolm approach (ggplot2 from his paper)
#======================================================#
library(tidymeta) # malcolm's package
library(ggplot2)
library(dplyr)
library(broom)
library(metafor)
library(rmeta)
library(reshape)
rm(list = ls())

# Functions
## Create table of case/control counts
## needs 'outcome' variable (define in data step). Group is typically studyname
counts <- function(data, group) {
  tmp_counts <- data %>% 
    dplyr::group_by_('outcome', group) %>% 
    dplyr::summarize(count = n()) %>% 
    tidyr::spread(key = outcome, value = count) %>% 
    dplyr::rename(Control = `0`, Case = `1`) %>% 
    dplyr::mutate(N = Control + Case)
}

## Conduct GLM by group (typically studyname for meta)
glm_by_group_df <- function(data, group, exposure, formula_txt) {
  formula.as.text <- paste0(formula_txt, "+", exposure)
  tmp_group <- data %>%
    group_by(studyname)
  tmp_results_beta <- do(tmp_group, tidy(glm(as.formula(formula.as.text), data = . , family = 'binomial')))
  tmp_results_ci   <- do(tmp_group, confint_tidy(glm(as.formula(formula.as.text), data = . , family = 'binomial'))) %>%
    filter(!is.na(conf.low))
  # i'm also creating variables you'd want to explore heterogeneity by - in this case, study design
  tmp_results <- bind_cols(tmp_results_beta, tmp_results_ci) %>%
    ungroup %>% 
    dplyr::filter(term == exposure) 
  tmp_results
}

## round with zeros
round_with_zeros <- function(.x, digits = 2) {
  format(round(.x, digits = digits), nsmall = digits)
}

Epi <- readRDS("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120_EpiOnly.rds")

# Data (please remove case-only studies)
gxe_set <- Epi %>% 
  filter(drop == 0 & gxe == 1) %>% 
  mutate(age_ref = as.numeric(age_ref),
         sex = as.factor(sex),
         energytot = as.numeric(energytot),
         folate_tot = as.numeric(folate_tot),
         folate_tot_400 = folate_tot / 400,
         folate_totqc2 = as.numeric(folate_totqc2),
         folate_diet = as.numeric(folate_diet),
         folate_diet_400 = folate_diet / 400,
         folate_dietqc2 = as.numeric(folate_dietqc2),
         folate_sup = as.numeric(folate_sup),
         folate_sup_400 = folate_sup / 400,
         studyname = ifelse(studyname %in% c("HPFS_1", "HPFS_2"), "HPFS_1_2", studyname),
         studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname),
         outcome = ifelse(outc == "Case", 1, 0),
         sex = ifelse(sex == "Female", 0, 1)) %>% 
  dplyr::select(vcfid, studyname, outc, outcome, age_ref, sex, energytot, starts_with('folate'))
table(gxe_set$studyname, gxe_set$outc)

## create studyname groups
exclude_studies <- c("GALEON", "MOFFITT", "NFCCR_1", "NGCCS", "ColoCare_1")
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DALS_1", "DALS_2", "Kentucky", "LCCS", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "NFCCR_2", "REACH_AD", "SELECT", "SMC_COSM", "SMS_AD")
cohort <- c("ATBC", "CLUEII", "CPSII_1", "CPSII_2", "EPIC", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD","MCCS_1", "MCCS_2","MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD","PPS3", "PPS4",  "UKB_1", "VITAL",  "WHI_1", "WHI_2", "WHI_3")


## create components for meta-analysis
## subsets change by what covariates you include e.g. ENERGYTOT
exposure = 'folate_tot_400'
covars =  c('age_ref', 'sex', 'energytot')
# covars =  c('age_ref', 'sex')
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))

dat <- gxe_set %>% 
  dplyr::select(outcome, exposure, covars, studyname) %>% 
  filter(!studyname %in% exclude_studies,
         complete.cases(.))

## case/control counts
r_counts <- counts(data = dat, group = 'studyname')
## casecontrol vs cohort
r_glm_meta <- glm_by_group_df(data = dat, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(design = ifelse(studyname %in% caco, 'CaseControl', 'Cohort')) %>% 
  group_by(design) %>% 
  meta_analysis(yi = estimate, sei = std.error, slab = studyname, exponentiate = TRUE) %>% ungroup() %>% 
  mutate(studyname = factor(study, levels = rev(study)),
         orci = paste0(round_with_zeros(estimate), " (",
                       round_with_zeros(conf.low), ", ",
                       round_with_zeros(conf.high), ")")) %>% 
  full_join(., r_counts, by = c("study" = "studyname")) %>% 
  mutate(counts = paste0(Case, " / ", Control, " (ctrl)"))


# test (make sure you have the right components in there)
#test <- data.frame(r_glm_meta[which(r_glm_meta$meta != "NULL"), 'meta'])[[2,1]] # just to check

# clean this up when you get a chance. 
breaks_subset <- r_glm_meta[which(r_glm_meta$meta != "NULL"), 'study', drop = T]
breaks_I2 <- filter(r_glm_meta, meta != "NULL") %>% 
  dplyr::select(meta) %>% 
  apply(., 1, function(x) x[[1]]$I2)
breaks_QEp <- filter(r_glm_meta, meta != "NULL") %>% 
  dplyr::select(meta) %>% 
  apply(., 1, function(x) x[[1]]$QEp)

breaks_h <- as.factor(paste0(breaks_subset, "\n", "(I-squared = ", round(breaks_I2, 2), "%, p = ", round(breaks_QEp, 3), ")"))
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
  mutate(xcoord = ifelse(facetvar == "orci", 5, 10),
         #xcoord = 10,
         value = ifelse(value == "NA/NA", "", value))

labels <- c(estimate = "Forest Plot", orci = "OR (95% CI)", counts = "Case/Control")
breaks_studyname <- levels(factor(r_glm_meta_h$study, levels = rev(r_glm_meta_h$study)))
faces <- as.expression(breaks_studyname)


ggplot(dat_ggplot_main, y = studyname) +
  theme_minimal() +
  geom_vline(data = vline_or_data, aes(xintercept = 1), linetype = 'dashed') +
  geom_point(data = dat_ggplot_orci, aes(x = estimate, y = studyname, col = design, size = weight, shape = type), alpha = .75) +
  geom_errorbarh(data = dat_ggplot_orci, aes(x = estimate, y = study, xmin = conf.low, xmax = conf.high, col = design), height = 0, size = .75, alpha = .75) +
  scale_x_continuous(trans = "log",
                     breaks = c(.5, 1, 2.0)) +
  geom_text(data = dat_ggplot_text, aes(x = xcoord, y = study, label = value), size = 3.5) + 
  scale_y_discrete(label = faces, breaks = breaks_studyname) +
  scale_shape_manual(values = c(15, 18)) + 
  # i gave up on faceting by columns too, can't center text on the facets, forget it, move on with life.
  facet_grid(design ~ ., scales = 'free', space = 'free', switch = "y", labeller = labeller(variable = labels)) + 
  theme(strip.text.y = element_text(angle = 180),
        strip.background = element_rect(color = "grey85", fill = "grey85"), # grey
        #strip.background = element_rect(color = "black", fill = "white", size = .75), # white with black borders
        #  remove grid lines
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_rect(color = "black", fill = NA, size = 1), # grid around panels
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 11),
        legend.position = "none") +
  scale_size_continuous(range = c(2, 8)) + 
  ggtitle("Total Folate (400mcg/day) - Study Design")

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_400_tidymeta.png", width = 8, height = 6, scale = 1.5)

