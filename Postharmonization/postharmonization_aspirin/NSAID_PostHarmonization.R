#======================================================#
# FIGI NSAID GxE - Post Harmonization
# 02/07/2019
# 
# 1) meta analyses
# 2) forest plots
# 3) funnel plots
#======================================================#
# library(tidymeta) # malcolm's package
library(tidyverse)
library(data.table)
library(broom)
library(metafor)
library(rmeta)
library(reshape)
library(rlang)

rm(list = ls())
# source("~/Dropbox/code/Functions/PostHarmonization.R")
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")


# sanity
# any(duplicated(Epi$vcfid))
# any(duplicated(Epi$pooledcompassid))
# table(Epi$studyname)
# table(Epi$V2)

#----------------------------------------------------------#
# create covariate file useful for post harmonization
# GxE Set - asp_ref (N = 102,792 --> 72,438)
#
# some questions on how to handle certain issues, see list
#----------------------------------------------------------#
pc30k <- fread("~/Dropbox/code/FIGI_Results/PrincipalComponents/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))

# case only studies (haven't removed yet)
table(Epi$studyname, Epi$outc)
exclude_studies <- c("GALEON", "MOFFITT", "NGCCS", "ColoCare_1") #  ColoCare_1 = 'case-series'?

df <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         smoke = ifelse(smoke == "", NA, smoke),
         smoke = ifelse(smoke == "Never smoker", 0, 
                        ifelse(smoke == "Former smoker", 1, 2)),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname), # combined NFCCR1 and NFCCR2
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname), # remove special character
         studyname = ifelse(studyname %in% c("HPFS_1", "HPFS_2"), "HPFS_1_2", studyname),
         studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA)),
         asp_ref2 = ifelse(asp_ref2 == "Yes", 1, ifelse(asp_ref2 == "No", 0, NA)),
         aspirin = ifelse(aspirin == "Yes", 1, ifelse(aspirin == "No", 0, NA)),
         nsaids = ifelse(nsaids == "Yes", 1, ifelse(nsaids == "No", 0, NA)),
         aspirinyrs_b = ifelse(aspirinyrs > 0, 1, 0),
         nsaidsyrs_b = ifelse(nsaidsyrs > 0, 1, 0),
         aspyrs_ref_b = ifelse(aspyrs_ref > 0, 1, 0))

table(df$studyname, df$outc)
# saveRDS(df, file = '~/df.rds')
# ignore variable 'aspirin_nsaids_ever' - too many missing
#table(aspirinyrs, aspirin_nsaids_ever, useNA = 'ifany')
#table(nsaidsyrs, aspirin_nsaids_ever, useNA = 'ifany')
#table(aspirin_nsaids_ever)
#table(aspirin, aspirin_nsaids_ever)
#table(nsaids, aspirin_nsaids_ever)


#------ Counts Table ------
# Make one for each aspirin/nsaid variable for the Rmd document
# by studyname, case/control. Include a 'totals' row
# include all studies for the summary
# please resist the temptation to mess with it (NSE and all)

getCounts <- function(data, group, exposure) {
  tmp <- enquo(group)
  
  castFormula = as.formula(paste0(group, "~", "outc", "+", exposure))
  
  countsTotal <- df %>% 
    dplyr::group_by_('outc') %>% 
    dplyr::select_(exposure, 'outc') %>% 
    dplyr::count_(exposure) %>% 
    dplyr::mutate(tot = sum(n),
                  pct = round(100*(n/tot), 1),
                  n_pct = paste0(n, " (", pct, ")"),
                  !! tmp := "Total") %>% dplyr::select(-n, -pct) %>% ungroup() %>% 
    reshape2::dcast(castFormula, fill = "", value.var = 'n_pct')
  
  countsGroup <- data %>%
    dplyr::group_by_(group, 'outc') %>%
    dplyr::select_(group, exposure, 'outc') %>%
    dplyr::count_(exposure) %>%
    dplyr::mutate(tot = sum(n),
                  pct = round(100*(n/tot), 1),
                  n_pct = paste0(n, " (", pct, ")")) %>% dplyr::select(-n, -pct) %>% ungroup() %>%
    reshape2::dcast(castFormula, fill = "", value.var = 'n_pct')
  
  results <- rbind(countsGroup, countsTotal)
  saveRDS(results, file = paste0('./tables/counts_', exposure, '.rds'))
}

# generate for Rmd document
getCounts(df, 'studyname', 'asp_ref')
getCounts(df, 'studyname', 'asp_ref2')
getCounts(df, 'studyname', 'aspirin')
getCounts(df, 'studyname', 'nsaids')


#------ Bar Plots ------
# simplify instead of getting 'cute'
getPlots <- function(data, group, exposure) {
  
  group_quo <- enquo(group)
  exposure_quo <- enquo(exposure)

  cov <- data %>%
    dplyr::select(outc, !! group_quo, !! exposure_quo) %>%
    dplyr::mutate(!! exposure_quo := factor(!! exposure_quo, exclude = NULL, labels = c("No", "Yes", "NA")),
                  outc = factor(outc, labels = c("Case", "Control"))) %>%
    group_by(!! group_quo)

  p <- ggplot(data=cov, aes(x=outc)) +
    geom_bar(aes(fill = !! exposure_quo), position = 'fill') +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    scale_fill_manual(values = c( "cyan", "red", "black")) +
    facet_wrap(vars(!! group_quo), ncol = 5)
  
  saveRDS(p, file = paste0("./figures/barplots_", deparse(substitute(exposure)), ".rds"))
}

getPlots(df, studyname, asp_ref)
getPlots(df, studyname, asp_ref2)
getPlots(df, studyname, aspirin)
getPlots(df, studyname, nsaids)


#------ Meta Analyses ------
# outcome is ALWAYS outcome
# group is ALWAYS studyname
#
# here.. to avoid saving data in DropBox, save the counts and GLM as rds objects

getCounts_outc <- function(data, group) {
  tmp_counts <- data %>% 
    dplyr::group_by_('outcome', group) %>% 
    dplyr::summarize(count = n()) %>% 
    tidyr::spread(key = outcome, value = count) %>% 
    dplyr::rename(Control = `0`, Case = `1`) %>% 
    dplyr::mutate(N = Control + Case)
}

glm_by_group_df <- function(data, group, exposure, formula_txt) {
  tmp_group <- data %>%
    group_by(studyname)
  tmp_results_beta <- do(tmp_group, tidy(glm(as.formula(formula_txt), data = . , family = 'binomial')))
  tmp_results_ci   <- do(tmp_group, confint_tidy(glm(as.formula(formula_txt), data = . , family = 'binomial')))
  tmp_results <- bind_cols(tmp_results_beta, tmp_results_ci) %>%
    ungroup %>% 
    dplyr::filter(term == exposure) %>% 
    dplyr::select(-studyname1)
  tmp_results
}

meta_fp <- function(counts_df, results_df, title) {
  tmp_meta <- meta.summaries(results_df$estimate, results_df$std.error, method = 'random')
  tmp_meta_b_ci <- c(tmp_meta[[3]], tmp_meta[[3]]-(tmp_meta[[4]]*1.96), tmp_meta[[3]]+(tmp_meta[[4]]*1.96))
  tmp_fp <- dplyr::select(results_df, estimate, conf.low, conf.high)
  tmp_fp <- rbind(rep(NA, 3), tmp_fp, rep(NA, 3), tmp_meta_b_ci)
  tmp_label_summary <- c("Summary", sum(counts_df$N), sum(counts_df$Case), round(tmp_meta[[3]], 2), round(tmp_meta[[4]], 3), round(exp(tmp_meta[[3]]), 2), round(tmp_meta[[5]][2], 4))
  tmp_label <- inner_join(counts_df, results_df, by = 'studyname') %>% 
    mutate(Study = as.character(studyname), 
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
  mtext(paste(exposure, '\n', "het.pval=", formatC(tmp_meta$het[3], format='e', digits=2)), side = 1, line = 1)
}

getMeta <- function(data, exposure, covars) {
  # create a tmp df that only includes covariates, complete case, and remove samples if within study strata becomes case only or control only
  metadf <- data %>% 
    dplyr::select(outcome, exposure, studyname, covars) %>% 
    dplyr::filter(complete.cases(.))
  cc <- data.frame(table(metadf$studyname, metadf$outcome)) %>% 
    filter(Freq == 0)
  metadf <- filter(metadf, !studyname %in% unique(cc$Var1)) %>% 
    mutate(studyname = factor(studyname))
  
  formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"), "+", exposure)
  
  r_counts <- getCounts_outc(data = metadf, group = 'studyname')
  saveRDS(r_counts, file = paste0('./meta/r_counts_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
  r_glm <- glm_by_group_df(metadf, 'studyname', exposure, formula_txt = formula_txt)
  saveRDS(r_glm, file = paste0('./meta/r_glm_', exposure, '_', paste(covars, collapse = "_"), '.rds'))
  #meta_fp(r_counts, r_glm, title = title) # find this function in the rmarkdown file NSAID_PostHarmonization.rmd
}

# covars <- c("age_ref_imp", "sex", paste0(rep("PC", 10), seq(1,10)), "smoke")
# covars <- c("age_ref_imp", "sex")
# covars <- c("age_ref_imp", "sex", "smoke")
# covars <- c("age_ref_imp", "sex")

covars <- c("age_ref_imp", "sex")
getMeta(df, 'asp_ref', covars = covars)
getMeta(df, 'aspirin', covars = covars)
getMeta(df, 'nsaids', covars = covars)


covars <- c("age_ref_imp", "sex", "smoke")
getMeta(df, 'asp_ref', covars = covars)


#------ MEGA Analyses ------
# maybe different adjustment covariates.. 









#------------------------------------------------------#
## create studyname groups
caco <- c("ASTERISK", "CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CORSA_1", "CORSA_2", "CRCGEN", "CzechCCS", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "EDRN", "EPICOLON", "ESTHER_VERDI","HawaiiCCS_AD","Kentucky", "LCCS" ,"MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "NFCCR_2", "PPS3", "PPS4", "REACH_AD", "SEARCH" ,  "SELECT" , "SMS_AD", "USC_HRT_CRC")

cohort <- c("ATBC", "CLUEII", "CPSII_1", "CPSII_2", "EPIC", "HPFS_1_2", "HPFS_3_AD","HPFS_4",  "HPFS_5_AD",  "MCCS_1" , "MCCS_2", "MEC_1","MEC_2", "NHS_1_2" ,"NHS_3_AD" ,  "NHS_4","NHS_5_AD", "PHS","PLCO_1"  ,"PLCO_2",  "PLCO_3", "PLCO_4_AD" , "SMC_COSM", "UKB_1" , "VITAL"   ,"WHI_1","WHI_2","WHI_3")

# American vs Non-American
american <- c("CCFR_1", "CCFR_3", "CCFR_4", "CLUEII", "Colo2&3", "CPSII_1", 
              "CPSII_2", "DALS_1", "DALS_2", "EDRN", "ESTHER_VERDI", "HawaiiCCS_AD", 
              "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "Kentucky", "MEC_1", 
              "MEC_2", "NCCCSI", "NCCCSII", "NFCCR_2", "NHS_1_2", "NHS_3_AD", 
              "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", 
              "PPS3", "PPS4", "REACH_AD", "SELECT", "SMS_AD", "USC_HRT_CRC", 
              "VITAL", "WHI_1", "WHI_2", "WHI_3")

nonamerican <- c("ASTERISK", "ATBC", "CORSA_1", "CORSA_2", "CRCGEN", "CzechCCS", 
                 "DACHS_1", "DACHS_2", "DACHS_3", "EPIC", "EPICOLON", "LCCS", "MCCS_1", "MCCS_2", "MECC_1", "MECC_2", 
                 "MECC_3", "SEARCH", "SMC_COSM", "UKB_1")

# Folate fortification period... pre, post, overlap
# STUDY VARIABLE
# prefort <- c("CLUEII", "DALS", "MEC", "PHS")
# postfort <- c("Kentucky", "NCCCSII", "NFCCR", "REACH", "SMS", "VITAL")
# overlap <- c("CCFR", "Colo23", "CPSII", "HPFS", "NCCCSI", "NHS", "PLCO", "WHI")


#------ Recreate asp_ref variable - check ------
# Just to double check
df <- dplyr::select(gxe_set, aspirin, nsaids, asp_ref, asp_ref2) %>% 
  mutate(asp_ref_test = ifelse(aspirin == "Yes" | nsaids == "Yes", "Yes", 
                        ifelse(aspirin == "No" & (nsaids == "" | nsaids == "No"), "No", 
                        ifelse(nsaids == "No" & (aspirin == "" | aspirin == "No"), "No", ""))))

table(df$asp_ref_test)
table(df$asp_ref)

# OK OK OK # 

df <- filter(gxe_set, asp_ref != "") # Number checks out
df <- filter(gxe_set, asp_ref2 != "") # Number checks out
73208 - 71429

#------------------------------------------------------#
# Create working dataset (subset of studies with available ANY FOLATE)
# make variables suitable for analysis
df <- gxe_set %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref = as.numeric(age_ref),
         sex.f = as.factor(sex),
         sex.n = ifelse(sex == "Female", 0, 1),
         energytot = as.numeric(energytot),
         folate_tot = as.numeric(folate_tot),
         folate_tot_400 = folate_tot / 400,
         folate_totqc2 = as.numeric(folate_totqc2),
         folate_diet = as.numeric(folate_diet),
         folate_diet_400 = folate_diet / 400,
         folate_dietqc2 = as.numeric(folate_dietqc2),
         folate_sup = as.numeric(folate_sup),
         folate_sup_400 = folate_sup / 400,
         folate_diet400qcm = as.numeric(folate_diet400qcm), 
         folate_tot_400nd = 1000*(folate_tot_400 / energytot), 
         folate_diet_400nd = 1000*(folate_diet_400 / energytot),
         folate_sup_400nd = 1000*(folate_sup_400 / energytot), 
         studyname = ifelse(studyname %in% c("HPFS_1", "HPFS_2"), "HPFS_1_2", studyname),
         studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname)) %>% 
  mutate(studydesign = ifelse(studyname %in% caco, 'CaseControl', 
                              ifelse(studyname %in% cohort, 'Cohort', '')),
         studycountry = ifelse(studyname %in% american, 'American', 
                               ifelse(studyname %in% nonamerican, 'NonAmerican', ''))) %>% 
  dplyr::select(vcfid, studyname, outc, outcome, age_ref, sex, sex.n, sex.f, energytot, starts_with('folate'), studydesign, studycountry)

table(df$studyname, df$outc)



#------ Meta-Analyses (Malcolm) ------

exposure = 'folate_tot'
covars =  c('age_ref', 'sex', 'energytot')
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))
plot_title = "Dietary Folate (Q4 median) by study design"
plot_subtitle = paste0("outcome~age_ref+sex+energytot+E")


dat <- gxe_set %>% 
  dplyr::select(outcome, exposure, covars, studyname) %>% 
  filter(!studyname %in% exclude_studies,
         complete.cases(.))




## case/control counts
# modify to add totals by subgroup
r_counts <- counts_outc(data = dat, group = 'studyname')

r_counts_all <- r_counts %>% 
  summarise_at(c('Control', 'Case', 'N'), funs(sum)) %>% 
  mutate(studyname = "Overall")

r_counts_design <- r_counts %>% 
  mutate(design = ifelse(studyname %in% caco, 'CaseControl', 'Cohort')) %>% 
  group_by(design) %>% 
  summarise_at(c('Control', 'Case', 'N'), funs(sum)) %>% ungroup() %>% 
  mutate(studyname = paste0("Subgroup: ", design))

r_counts_use <- bind_rows(r_counts, r_counts_all, r_counts_design) %>%  dplyr::select(-design)

## casecontrol vs cohort
r_glm_meta <- glm_by_group_df(data = dat, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(design = ifelse(studyname %in% caco, 'CaseControl', 'Cohort')) %>% 
  group_by(design) %>% 
  meta_analysis(yi = estimate, sei = std.error, slab = studyname, exponentiate = TRUE) %>% ungroup() %>% 
  mutate(studyname = factor(study, levels = rev(study)),
         orci = paste0(round_with_zeros(estimate), " (",
                       round_with_zeros(conf.low), ", ",
                       round_with_zeros(conf.high), ")")) %>% 
  full_join(., r_counts_use, by = c("study" = "studyname")) %>% 
  mutate(counts = paste0(Case, " / ", Control))


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
  mutate(xcoord = ifelse(facetvar == "orci", 20, 500))

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

ggsave(paste0("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_", exposure, "_tidymeta.png"), width = 10, height = 6, scale = 1.5)


# outcomes counts
#------ OUTCOME TABLE ------
dat <- filter(gxe_set, !studyname %in% exclude_studies)
counts_df <- counts_control_case_adv(data = gxe_set, group = 'studyname')
saveRDS(counts_df, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_counts.rds")


# year ref?
x <- filter(Epi, folate_tot != "")

xx <- data.frame(table(x$study, x$year_ref)) %>% 
  filter(Freq != 0) %>% 
  arrange(Var1)


# 


# Heterogeneity

#------ Heterogeneity ------
# inverse variance weighted p value.. 
exposure = 'folate_tot_500'
hetero <- metadat %>% 
  filter(studyname %in% caco_cohort)

covars = c('age_ref', 'sex', 'energytot')
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))


# by study design 
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_tot_500_heterogeneity_studydesign.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ study_design, data = r_glm, weights = w)
summary(x)


# by study country
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(country = ifelse(studyname %in% other, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_tot_500_heterogeneity_usa.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ country, data = r_glm, weights = w)
summary(x)




# inverse variance weighted p value.. 
exposure = 'folate_sup_500'
hetero <- metadat %>% 
  filter(studyname %in% caco_cohort,
         !is.na(folate_sup_500))

covars = c('age_ref', 'sex', 'energytot')
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))


# by study design 
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_sup_500_heterogeneity_studydesign.rds")


# by study country
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(country = ifelse(studyname %in% other, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_sup_500_heterogeneity_usa.rds")



# inverse variance weighted p value.. 
exposure = 'folate_diet_500'
hetero <- metadat %>% 
  filter(studyname %in% caco_cohort,
         !is.na(folate_diet_500))

covars = c('age_ref', 'sex', 'energytot')
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))


# by study design 
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_diet_500_heterogeneity_studydesign.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ study_design, data = r_glm, weights = w)
summary(x)


# by study country
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(country = ifelse(studyname %in% other, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_diet_500_heterogeneity_usa.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ country, data = r_glm, weights = w)
summary(x)













############################

# NCCCSII very high aspirin intake

ncccsii <- filter(gxe_set, studyname == "NCCCSII")
table(ncccsii$aspirin, ncccsii$outc)
table(ncccsii$aspirinyrs_b, ncccsii$aspirin, useNA = 'ifany')

hist(as.numeric(ncccsii[which(ncccsii$asp_ref == "Yes"),]$aspirinyrs))



ccfr1 <- filter(gxe_set, studyname == "CCFR_1")
table(ccfr1$aspirin, ccfr1$outc)
table(ccfr1$aspirinyrs_b, ccfr1$aspirin, useNA = 'ifany')

hist(as.numeric(ccfr1[which(ccfr1$asp_ref == "Yes"),]$aspirinyrs))
