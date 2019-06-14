#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI POST HARMONIZATION 11/26/2018
#
# Summary of E Variable - Folate
# Meta Analyses
# Study Heterogeneity
#
# use outputs from here to generate rmarkdown reports
#
# For summary, let's use epi gxe set (gxe == 1)
# change if needed
#
################################################################
# Riki suggestions:
# - investigate distribution of E across studies
#   - check this for each study and study set (when study is included multiple sets, check if distribution of E is similar across sets?)
#   - are there outliers first check the study design, for instance is it a study in smokers or a trial for a specific E and the way the E
#     was assessed by looking at what question were used by checking the mapping document and/or used questionnaire (available on the portal)
#   - Are the overall values within normal range
#   - what iwll you do with outlier? reassign a logical maximum, exclude?
#   - if there is no explanation check with Yi
# - Discuss and evaluate the best approach to code the existing variables
#   - explore which coding schema of the E variable provides the least heterogeneity (I2 and p for heterog), 
#     strongest effect size, most significant p-value
# - please note that we cannot go back to stuies for harmonization of additional variables given available resources and the very 
#   large number of studies involved
# - investigate distribution by know/ well established effect modifiers (sex is usual, could be others like alcohol and folate)
# - check for heterogeneity of effect
#   - evaluate if any factors, such as study design, sex, year of enrollment, country/continent, can explain some or all of the heterogeneity
#
# Andre
# (look, no one cares if your graphs and tables are pretty, they care if they're informative)
# - assess interactions with alcohol by study?
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(broom)
library(data.table)
library(rmeta)

rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData") # newest Epi data subset
source("~/Dropbox/code/Functions_AK/FIGI_PostHarmonization_Functions.R")

#------ epidat ------
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
         studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname)) 

summary(gxe_set$folate_tot)
summary(gxe_set$folate_diet)
summary(gxe_set$folate_sup)



ca <- c("ColoCare_1", "NFCCR_1")
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DALS_1", "DALS_2", "Kentucky", "LCCS", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "NFCCR_2", "SELECT", "SMC_COSM")
cohort <- c("ATBC", "CLUEII", "CPSII_1", "CPSII_2", "EPIC", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD","MCCS_1", "MCCS_2","MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD","PPS3", "PPS4",  "UKB_1", "VITAL",  "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)
ca_caco_cohort <- c(ca, caco, cohort)

# American vs others
usa <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "DALS_1", "DALS_2", "Kentucky", "NCCCSI", "NCCCSII",  "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD",  "MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "NFCCR_2", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "VITAL", "WHI_1", "WHI_2", "WHI_3")
other <- c("CRCGEN", "EPIC", "LCCS", "MECC_1", "MECC_2", "MECC_3", "ATBC", "MCCS_1", "MCCS_2", "SMC_COSM", "UKB_1")


#------ META ANALYSIS ------
# model outc to match Yi shinyapp
metadat <- gxe_set %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0), 
         sex = ifelse(sex == "Female", 0, 1)) %>% 
  dplyr::select(vcfid, studyname, outc, outcome, age_ref, sex, energytot, starts_with('folate'))


#------ Helper Function ------
meta_fp_helper <- function(data, exposure, covars, study_subset, png_name, width = 700, height = 800, title) {
  formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))
  dat <- metadat %>% 
    filter(studyname %in% study_subset) %>% 
    dplyr::select(outcome, exposure, covars, studyname) %>% 
    filter(complete.cases(.))
    r_counts <- counts(data = dat, group = 'studyname')
  r_glm <- glm_by_group_df(data = dat, group = 'studyname', exposure = exposure, formula_txt = formula_txt)
 
  png(png_name, width = width, height = height)
  meta_fp(r_counts, r_glm, title = title)
  dev.off()
}

# parameters in common
covars <- c('age_ref', 'sex', 'energytot')

# CaseControl + Cohort ------
# Total Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_tot_400',
               study_subset = caco_cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_caco_cohort.png",
               width = 600, 
               height = 850, 
               title = 'Total Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_tot_400\nCaseControl + Cohort studies')

# Dietary Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_diet_400',
               study_subset = caco_cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_fp_caco_cohort.png",
               width = 600, 
               height = 850, 
               title = 'Dietary Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_diet_400\nCaseControl + Cohort studies')

# Supplemental Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_sup_400',
               study_subset = caco_cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_sup_fp_caco_cohort.png",
               width = 600, 
               height = 750, 
               title = 'Supplemental Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_sup_400\nCaseControl + Cohort studies')

# CaseControl ------
# Total Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_tot_400',
               study_subset = caco, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_caco.png",
               width = 600, 
               height = 500, 
               title = 'Total Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_tot_400\nCaseControl studies')

# Dietary Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_diet_400',
               study_subset = caco, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_fp_caco.png",
               width = 600, 
               height = 500, 
               title = 'Dietary Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_diet_400\nCaseControl studies')

# Supplemental Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_sup_400',
               study_subset = caco, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_sup_fp_caco.png",
               width = 600, 
               height = 500, 
               title = 'Supplemental Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_sup_400\nCaseControl studies')

# Cohort ------
# Total Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_tot_400',
               study_subset = cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_cohort.png",
               width = 600, 
               height = 600, 
               title = 'Total Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_tot_400\nCohort studies')

# Dietary Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_diet_400',
               study_subset = cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_fp_cohort.png",
               width = 600, 
               height = 600, 
               title = 'Dietary Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_diet_400\nCohort studies')

# Supplemental Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_sup_400',
               study_subset = cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_sup_fp_cohort.png",
               width = 600, 
               height = 600, 
               title = 'Supplemental Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_sup_400\nCohort studies')



# USA ------
# Total Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_tot_400',
               study_subset = usa, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_usa.png",
               width = 600, 
               height = 700, 
               title = 'Total Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_tot_400\nNorthAmerican studies')

# Dietary Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_diet_400',
               study_subset = usa, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_fp_usa.png",
               width = 600, 
               height = 700, 
               title = 'Dietary Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_diet_400\nNorthAmerican studies')

# Supplemental Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_sup_400',
               study_subset = usa, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_sup_fp_usa.png",
               width = 600, 
               height = 700, 
               title = 'Supplemental Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_sup_400\nNorthAmerican studies')



# Other ------
# Total Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_tot_400',
               study_subset = other, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_other.png",
               width = 600, 
               height = 500, 
               title = 'Total Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_tot_400\nNon-American studies')

# Dietary Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_diet_400',
               study_subset = other, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_fp_other.png",
               width = 600, 
               height = 500, 
               title = 'Dietary Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_diet_400\nNon-American studies')

# Supplemental Folate
meta_fp_helper(data = metadat,
               covars = covars,
               exposure = 'folate_sup_400',
               study_subset = other, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_sup_fp_other.png",
               width = 600, 
               height = 500, 
               title = 'Supplemental Folate (400 mcg/day)\nModel: outc~age_ref+sex+energytot+folate_sup_400\nNon-American studies')



#------ OUTCOME TABLE ------
dat <- filter(metadat, studyname %in% c(caco, cohort))
counts_df <- counts_control_case_adv(data = dat, group = 'studyname')
saveRDS(counts_df, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_counts.rds")


#------ OUTCOME PLOT ------
# plot percentage case/control/adv
# factor level order matters - sorting by # of controls
# tmp <- sort(unique(epidat$studyname), decreasing = T)

ggdat <- metadat %>% 
  mutate(outcome = ifelse(outc == "Control", 0, 
                          ifelse(adenoma_adv == "Case", 2, 1))) %>% 
	filter(studyname %in% caco_cohort) %>% 
	group_by(studyname, outcome) %>% 
	summarise(count = n()) %>% 
	mutate(perc=count/sum(count)) %>% 
	ungroup() %>% 
	mutate(outcome_f = factor(outcome, levels = c(2,1,0), labels = c("AdvAd", "CRC", "Control"))) %>% 
	arrange(outcome, perc)

sort_h <- rev(unique(ggdat$studyname)) # create vector to sort factor var
ggdat <- ggdat %>% 
	mutate(studyname_f = factor(studyname, levels = sort_h))

p <- ggplot(aes(y = perc, x = studyname_f), data = ggdat) + 
	geom_bar(aes(fill = outcome_f), stat='identity', alpha = 0.6, width = 0.7) +
	coord_flip() + 
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank()) +
	scale_fill_manual(values = c("grey", "red", "blue"))

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_outcome_barplot.png", plot = p, device = "png", width = 5, height = 5, scale = 1.2)



#------ GGPLOT ----
ggdat <- metadat %>% 
	filter(studyname %in% caco_cohort,
	       !is.na(folate_tot)) %>% 
	mutate(sex = factor(sex, labels = c("Female", "Male")),
				 #legend label for outcome
				 outcomep = ifelse(outcome == 0, 'Control', 'Case'))
ggdat2 <- ggdat %>% 
	group_by(studyname) %>% 
	summarise(folate_med = median(folate_tot)) %>% 
	ungroup() %>% arrange(folate_med)
sort_h <- c(rev(unique(ggdat2$studyname)))
test <- ggdat %>% 
	mutate(studyname_f = factor(studyname, levels = sort_h))


# for folate, lots of outlying values? 
# let's try removing top 99th percentile in all data
# high_folate <- quantile(tmp$folate_tot, c(.99))
# tmp <- filter(tmp, folate_tot < high_folate)
p <- ggplot(aes(y = folate_tot, x = outcomep, fill = outcomep), data = test) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_boxplot.png", plot = p, device = "png", width = 5, height = 8, scale = 1.2)

# case-control only
test2 <- filter(test, studyname %in% caco)
p <- ggplot(aes(y = folate_tot, x = outcomep, fill = outcomep), data = test2) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
  ggtitle("Case-Control") + 
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_boxplot_caco.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)

# Cohort
test2 <- filter(test, studyname %in% cohort)
p <- ggplot(aes(y = folate_tot, x = outcomep, fill = outcomep), data = test2) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
  ggtitle("Cohort") +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_boxplot_cohort.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)

# USA
test2 <- filter(test, studyname %in% usa)
p <- ggplot(aes(y = folate_tot, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("USA") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_boxplot_usa.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)

# Other
test2 <- filter(test, studyname %in% other)
p <- ggplot(aes(y = folate_tot, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("Other") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_boxplot_other.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)


# # GGPLOT Folate - controls only, visualize study_design ----
# 
# # sort table by median folate_tot levels (grouped by study)
# ggdat2 <- ggdat %>% 
# 	filter(outcome == 0) %>% 
# 	group_by(studyname) %>% 
# 	summarise(folate_med = median(folate_tot)) %>% 
# 	arrange(folate_med) %>% 
# 	ungroup() %>% arrange(folate_med)
# 
# sort_h <- c(rev(unique(ggdat2$studyname)))
# 
# test <- ggdat %>% 
# 	mutate(studyname_f = factor(studyname, levels = sort_h), 
# 				 study_design = factor(ifelse(studyname %in% caco, 'case-control', 'cohort')))
# 
# p <- ggplot(aes(y = as.numeric(folate_tot_500), x = study_design, fill = study_design), data = test) + 
# 	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
# 	coord_flip() +
# 	theme_bw() + 
# 	theme(
# 		panel.grid.major.x = element_blank(),
# 		panel.grid.minor.x = element_blank(),
# 		legend.position = 'bottom') + 
# 	scale_fill_manual(values = c("green", "orange"))
# 
# ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot_control_studydesign.png", plot = p, device = "png", width = 5, height = 6, scale = 1.2)
# 


# Heterogeneity ----
# inverse variance weighted p value.. 
exposure = 'folate_tot_500'
hetero <- metadat %>% 
	filter(studyname %in% caco_cohort) %>% 
	dplyr::select(folate_tot_500, c('outcome', 'age_ref', 'sex', 'energytot', 'studyname')) %>% 
	filter(complete.cases(.))

# by study design
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = "outcome ~ age_ref + sex + energytot + ") %>% 
	mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/heterogeneity_studydesign.rds")

w <- (r_glm$std.error)^2
x <- lm(estimate ~ study_design, data = r_glm, weights = w)
summary(x)


# by study country
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = "outcome ~ age_ref + sex + energytot + ") %>% 
	mutate(country = ifelse(studyname %in% other, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/heterogeneity_usa.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ country, data = r_glm, weights = w)
summary(x)




#------ Violin ------
# remember test is defined above!
# all, caco, cohort (x2 - total + dietary)


# folat_tot - caco vs cohort
test2 <- filter(test, studyname %in% caco_cohort) 
p1 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
	geom_violin() +
	geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("CaseControl + Cohort") +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'none')
	
test2 <- filter(test, studyname %in% caco) 
p2 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
	geom_violin() +
	geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("CaseControl") +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'none')

test2 <- filter(test, studyname %in% cohort)
p3 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
	geom_violin() +
	geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("Cohort") +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_violin_caco_cohort_filter.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)

# folat_tot - usa vs other
test2 <- filter(test, studyname %in% caco_cohort) 
p1 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("CaseControl + Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% usa) 
p2 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("USA") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% other)
p3 <- ggplot(aes(y = folate_tot, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2500) + 
  ggtitle("Other") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_violin_usa_other_filter.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)




# folate_diet - caco vs cohort
test2 <- filter(test, studyname %in% caco_cohort) 
p1 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("CaseControl + Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% caco) 
p2 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("CaseControl") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% cohort)
p3 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_violin_caco_cohort_filter.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)



# folate_diet - usa other
test2 <- filter(test, studyname %in% caco_cohort) 
p1 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("CaseControl + Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% usa) 
p2 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("USA") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% other)
p3 <- ggplot(aes(y = folate_diet, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(0, 2000) + 
  ggtitle("Other") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_violin_usa_other_filter.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)




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





#------ MEGA analysis... ------#
dat <- metadat %>% 
  filter(studyname %in% caco_cohort) %>% 
	mutate(studyname = as.factor(studyname))

saveRDS(dat, file = "~/Dropbox/code/FIGI_PostHarmonization/working/folate_mega.rds")


x <- glm(outcome ~ sex + age_ref + energytot + folate_tot_500 + studyname, data = dat, family = 'binomial')
tidy(x) %>% 
  mutate(or = exp(estimate))
glance(x)
coef(summary(x))
exp(coef(summary(x)))

x <- glm(outcome ~ sex + age_ref + energytot + folate_diet_500 + studyname, data = dat, family = 'binomial')
exp(coef(summary(x)))

x <- glm(outcome ~ sex + age_ref + energytot + folate_sup_500 + studyname, data = dat, family = 'binomial')
exp(coef(summary(x)))
