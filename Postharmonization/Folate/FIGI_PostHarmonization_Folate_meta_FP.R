#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI POST HARMONIZATION 10/04/2018
#
# Summary of E Variable - Folate
# Meta Analyses
# Study Heterogeneity
#
# use outputs from here to generate rmarkdown reports
#
# For summary, let's use epi gxe set (gxe == 1)
# change if needed
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(tidyverse)
library(broom)
library(data.table)
library(rmeta)

rm(list = ls())
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_11052018.RData") # newest Epi data subset
source("~/Dropbox/code/Functions_AK/FIGI_PostHarmonization_Helper_continuousE.R")

# STUDY SUBSETS ----
# In my set, I'm missing the following studies: BWHS, DACHS_4
# Exclude the following: ASTERISK, CCFR_2, CGN, ColoCare_1, Colo2&3, CORSA_1, CORSA_2, Czech_CCS, EPIC, EPICOLON, ESTHER_VERDI, FIRE3, GALEON, HawaiiCCS_AD,  HCES-CRC, HispanicCCS, HPFS_1, HPFS_2, LCCS, MAVERICC, MOFFITT, MSKCC, NFCCR_1, NGCCS, NHSII, PMH-CCFR, PURIFICAR, SLRCCS, SMC_COSM, SMHS, SWHS, Taiwan, TRIBE
# also exclude NFCCR_2
# See Riki email to Yi -- some clinical trials should be treated as COHORTS, unless it's the exposure they're trialing
# PPS 3/4 a bit confusing because one of them used aspirin/folate as the trial intervention
# For now, can include as cohort why not 
#cat(unique(sort(Epi$studyname)), sep = '", "')
#clintrials <- c("ATBC", "PPS3", "PPS4", "SELECT")
#cohort <- c("CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")

exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPIC", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMC_COSM", "SMHS", "SWHS", "Taiwan", "TRIBE")
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "Kentucky", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "REACH_AD", "SMS_AD", "USC_HRT_CRC")
cohort <- c("ATBC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)


# testing
exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMHS", "SWHS", "Taiwan", "TRIBE")
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "Kentucky", "LCCS", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "REACH_AD", "SMC_COSM", "SMS_AD", "USC_HRT_CRC")
cohort <- c("ATBC", "EPIC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)




# verified against newest shinyapp
exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPIC", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMC_COSM", "SMHS", "SWHS", "Taiwan", "TRIBE")

caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DALS_1", "DALS_2", "Kentucky", "LCCS", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "NFCCR_2", "REACH_AD", "SMC_COSM", "SMS_AD")
cohort <- c("ATBC", "EPIC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)



nrow(Epi[which(Epi$gxe == 1),])  # number of samples in gxe set, N = 102,799

epidat <- Epi %>% 
	filter(drop == 0,
				 gxe == 1,
				 !studyname %in% c("NFCCR_1", "ColoCare_1")) %>% 
	mutate(studyname = ifelse(studyname %in% c("HPFS_1", "HPFS_2"), "HPFS_1_2", studyname),
				 studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname)) %>% 
	mutate(outcome = ifelse(outc == "Control", 0, 
									 ifelse(adenoma_adv == "Case", 2, 1))) # creating one variable for crc + adv adenoma outcomes


#------ OUTCOME TABLE ------
dat <- filter(epidat, studyname %in% caco_cohort)
counts_df <- counts_control_case_adv(data = dat, group = 'studyname')
saveRDS(counts_df, file = "~/Dropbox/code/FIGI_PostHarmonization/working/counts.rds")


#------ OUTCOME PLOT ------
# plot percentage case/control/adv
# factor level order matters - sorting by # of controls
# tmp <- sort(unique(epidat$studyname), decreasing = T)
case_only <- c("GALEON", "MOFFITT", "NFCCR_1", "NGCCS")

ggdat <- epidat %>% 
	group_by(studyname, outcome) %>% 
	summarise(count = n()) %>% 
	mutate(perc=count/sum(count)) %>% 
	ungroup() %>% 
	mutate(outcome_f = factor(outcome, levels = c(2,1,0), labels = c("AdvAd", "CRC", "Control"))) %>% 
	arrange(outcome, perc)

ggdat_cc <- filter(ggdat, !studyname %in% case_only)
sort_h <- c(rev(unique(ggdat_cc$studyname)), case_only) # create vector to sort factor var

ggdat <- ggdat %>% 
	mutate(studyname_f = factor(studyname, levels = sort_h))

p <- ggplot(aes(y = perc, x = studyname_f), data = ggdat) + 
	geom_bar(aes(fill = outcome_f), stat='identity', alpha = 0.6) +
	coord_flip() + 
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank()) +
	scale_fill_manual(values = c("grey", "red", "blue"))

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_outcome_barplot.png", plot = p, device = "png", width = 9, height = 11, scale = 1.2)


#------ META ANALYSIS ------
# "other" does NOT necessarily mean adenoma_adv (make sure)
# DROP "other"
# redefine outcome so it's 0/1 only (collapse crc + adenoma_adv)
# need to remove energytot outliers for main analysis, but keep for postharm summary
# Complete Case only - helper function takes care of that
 
metadat <- epidat %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0), 
				 sex = ifelse(sex == "Female", 0, 1),
				 age_ref = as.numeric(age_ref),
				 energytot = as.numeric(energytot), 
				 folate_tot = as.numeric(folate_tot),
				 folate_diet = as.numeric(folate_diet),
				 folate_sup = as.numeric(folate_sup), 
				 #log_folate_tot = log(folate_tot),
				 folate_diet_500 = folate_diet/500,
				 folate_sup_500 = folate_sup/500,
				 folate_tot_500 = folate_tot/500) %>% 
	dplyr::select(outcome, sex, age_ref, energytot, starts_with('folate'), studyname)

x <- metadat %>% 
	dplyr::select(outc, studyname, folate_tot) %>% 
	filter(complete.cases(.))
table(x$studyname)
# because number of NA's varies by specific folate variable, will need to tailor each analyses accordingly!!
summary(metadat$energytot) 
summary(metadat$folate_tot)
summary(metadat$folate_diet)
summary(metadat$folate_sup)

#------ histograms for continuous variables? ------
# png("~/git/FIGI_PostHarmonization/working/figi_postharm_folate_tot_histogram_caco_cohort.png", width = 600, height = 400)
# hist(epidat$folate_tot)
# dev.off()
# 
# png("~/git/FIGI_PostHarmonization/working/figi_postharm_folate_tot_500_histogram_caco_cohort.png", width = 600, height = 400)
# hist(epidat$folate_tot_500)
# dev.off()
# 
# png("~/git/FIGI_PostHarmonization/working/figi_postharm_logfolate_tot_histogram_caco_cohort.png", width = 600, height = 400)
# hist(epidat$log_folate_tot)
# dev.off()
#------ Helper Function ------
# helper function to make code less long
# study_attrib for now is design (caco, cohort, caco_cohort)
# png_name is name of output forestplot
meta_fp_helper <- function(data, exposure, covars, formula_txt, study_attrib, png_name, width = 700, height = 800, title) {
	dat <- metadat %>% 
		filter(studyname %in% study_attrib) %>% 
		dplyr::select(exposure, covars) %>% 
		filter(complete.cases(.))
	r_counts <- counts(data = dat, group = 'studyname')
	r_glm <- glm_by_group_df(data = dat, group = 'studyname', exposure = exposure, formula_txt = formula_txt)
	png(png_name, width = width, height = height)
	meta_fp(r_counts, r_glm, title = title)
	dev.off()
}

# total folate no rescale
meta_fp_helper(data = metadat,
							 exposure = 'folate_tot',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco_cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_fp_caco_cohort.png",
							 width = 700, 
							 height = 800, 
							 title = 'Total Folate Intake\nModel: outc~age_ref+sex+energytot+folate_tot\nCaseControl + Cohort studies')


# total folate
meta_fp_helper(data = metadat,
							 exposure = 'folate_tot_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco_cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_500_fp_caco_cohort.png",
							 width = 700, 
							 height = 800, 
							 title = 'Total Folate Intake\nModel: outc~age_ref+sex+energytot+folate_tot_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_tot_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_500_fp_caco.png",
							 width = 700, 
							 height = 400, 
							 title = 'Total Folate Intake\nModel: outc~age_ref+sex+energytot+folate_tot_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_tot_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_tot_500_fp_cohort.png",
							 width = 700, 
							 height = 600, 
							 title = 'Total Folate Intake\nModel: outc~age_ref+sex+energytot+folate_tot_500\nCaseControl + Cohort studies')



# Dietary folate
meta_fp_helper(data = metadat,
							 exposure = 'folate_diet_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco_cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_500_fp_caco_cohort.png",
							 width = 700, 
							 height = 800, 
							 title = 'Dietary Folate Intake\nModel: outc~age_ref+sex+energytot+folate_diet_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_diet_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_500_fp_caco.png",
							 width = 700, 
							 height = 400, 
							 title = 'Dietary Folate Intake\nModel: outc~age_ref+sex+energytot+folate_diet_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_diet_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_diet_500_fp_cohort.png",
							 width = 700, 
							 height = 600, 
							 title = 'Dietary Folate Intake\nModel: outc~age_ref+sex+energytot+folate_diet_500\nCaseControl + Cohort studies')



# Supplemental folate
meta_fp_helper(data = metadat,
							 exposure = 'folate_sup_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco_cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_supp_500_fp_caco_cohort.png",
							 width = 700, 
							 height = 800, 
							 title = 'Supplemental Folate Intake\nModel: outc~age_ref+sex+energytot+folate_sup_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_sup_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = caco, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_supp_500_fp_caco.png",
							 width = 700, 
							 height = 400, 
							 title = 'Supplemental Folate Intake\nModel: outc~age_ref+sex+energytot+folate_sup_500\nCaseControl + Cohort studies')

meta_fp_helper(data = metadat,
							 exposure = 'folate_sup_500',
							 covars = c('outcome', 'age_ref', 'sex', 'energytot', 'studyname'), 
							 formula_txt = 'outcome ~ age_ref + sex + energytot + ',
							 study_attrib = cohort, 
							 png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_supp_500_fp_cohort.png",
							 width = 700, 
							 height = 600, 
							 title = 'Supplemental Folate Intake\nModel: outc~age_ref+sex+energytot+folate_sup_500\nCaseControl + Cohort studies')



# Heterogeneity ----
# inverse variance weighted p value.. 
exposure = 'folate_tot_500'
hetero <- metadat %>% 
	filter(studyname %in% caco_cohort) %>% 
	dplyr::select(folate_tot_500, c('outcome', 'age_ref', 'sex', 'energytot', 'studyname')) %>% 
	filter(complete.cases(.))

r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = "outcome ~ age_ref + sex + energytot + ") %>% 
	mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/git/FIGI_PostHarmonization/working/heterogeneity_studydesign.rds")

w <- (r_glm$std.error)^2
x <- lm(estimate ~ study_design, data = r_glm, weights = w)
summary(x)



# GGPLOT ----
ggdat <- Epi %>% 
	filter(outc != "Other",
				 gxe == 1,
				 studyname %in% caco_cohort) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0),
				 sex = ifelse(sex == "Female", 0, 1),
				 folate_tot = as.numeric(folate_tot),
				 folate_tot = as.numeric(folate_tot),
				 log_folate_tot = log(folate_tot),
				 folate_tot_500 = folate_tot/500,
				 folate_tot_400 = folate_tot/400,
				 outcomep = ifelse(outcome == 0, 'Control', 'Case')) %>% 
	filter(!is.na(folate_tot)) 


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

p <- ggplot(aes(y = as.numeric(folate_tot_500), x = outcomep, fill = outcomep), data = test) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
p
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot.png", plot = p, device = "png", width = 6, height = 12, scale = 1.5)


# case-control only
test2 <- filter(test, studyname %in% caco)
p <- ggplot(aes(y = as.numeric(folate_tot_500), x = outcomep, fill = outcomep), data = test2) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
p
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot_caco.png", plot = p, device = "png", width = 6, height = 8, scale = 1.5)


# Cohort
test2 <- filter(test, studyname %in% cohort)
p <- ggplot(aes(y = as.numeric(folate_tot_500), x = outcomep, fill = outcomep), data = test2) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom')
p
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot_cohort.png", plot = p, device = "png", width = 6, height = 8, scale = 1.5)





# GGPLOT Folate - controls only, visualize study_design ----

# sort table by median folate_tot levels (grouped by study)
ggdat2 <- ggdat %>% 
	filter(outcome == 0) %>% 
	group_by(studyname) %>% 
	summarise(folate_med = median(folate_tot)) %>% 
	ungroup() %>% arrange(folate_med)

sort_h <- c(rev(unique(ggdat2$studyname)))

test <- ggdat %>% 
	mutate(studyname_f = factor(studyname, levels = sort_h), 
				 study_design = factor(ifelse(studyname %in% caco, 'case-control', 'cohort')))

p <- ggplot(aes(y = as.numeric(folate_tot_500), x = study_design, fill = study_design), data = test) + 
	geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
	coord_flip() +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		legend.position = 'bottom') + 
	scale_fill_manual(values = c("green", "orange"))
p
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot_control_studydesign.png", plot = p, device = "png", width = 6, height = 7, scale = 1.2)




#------ Violin ------







#------ MEGA analysis... ------#
dat <- filter(epidat, studyname %in% caco_cohort) %>% 
	mutate(studyname = as.factor(studyname))

x <- glm(outcome ~ sex + age_ref + energytot + folate_tot_500 + studyname, data = dat, family = 'binomial')

exp(coef(summary(x)))


