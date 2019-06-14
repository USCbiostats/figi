#-----------------------------------------------------------
# FIGI POST HARMONIZATION Analysis 10/04/2018
#
# Summary of E Variable - ALCOHOL
# Meta Analyses
# Study Heterogeneity
#
# use outputs from here to generate rmarkdown reports
#-----------------------------------------------------------
library(tidyverse)
library(broom)
library(data.table)
library(rmeta)

rm(list = ls())
setwd("~/git/FIGI_PostHarmonization/")
load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_10182018.RData")
source("~/git/Functions_AK/FIGI_PostHarmonization_Helper_continuousE.R")


#------------------------------------------------------------
# In my set, I'm missing the following studies: BWHS, DACHS_4
#
#
# Exclude the following: ASTERISK, CCFR_2, CGN, ColoCare_1, Colo2&3, CORSA_1, CORSA_2, Czech_CCS, EPIC, EPICOLON, ESTHER_VERDI, FIRE3, GALEON, HawaiiCCS_AD,  HCES-CRC, HispanicCCS, HPFS_1, HPFS_2, LCCS, MAVERICC, MOFFITT, MSKCC, NFCCR_1, NGCCS, NHSII, PMH-CCFR, PURIFICAR, SLRCCS, SMC_COSM, SMHS, SWHS, Taiwan, TRIBE
#
# also exclude NFCCR_2
#------------------------------------------------------------

exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPIC", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMC_COSM", "SMHS", "SWHS", "Taiwan", "TRIBE")
#cat(unique(sort(Epi$studyname)), sep = '", "')
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "Kentucky", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "REACH_AD", "SMS_AD", "USC_HRT_CRC")
cohort <- c("CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
clintrials <- c("ATBC", "PPS3", "PPS4", "SELECT")

# See Riki email to Yi -- some clinical trials should be treated as cohorts, unless it's the exposure they're trialing
# PPS 3/4 a bit confusing because one of them used aspirin/folate as the trial intervention
# For now, can include as cohort why not 
cohort <- c("ATBC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)


#------ Counts data.frame ------#
nrow(Epi[which(Epi$gxe == 1),]) # number of samples in gxe analysis
# N = 103157

epidat <- samplefile %>% 
	filter(drop == 0,
				 gxe == 1) %>% 
	mutate(outcome = ifelse(outc == "Control", 0, 
				 					 ifelse(adenoma_adv == "Case", 2, 1)))

dat <- filter(epidat, studyname %in% caco_cohort)
counts_df <- counts_control_case_adv(data = dat, group = 'studyname')
saveRDS(counts_df, file = "~/git/FIGI_PostHarmonization/working/counts_etoh.rds")



#------ meta analysis and FP ------#
# complete case analyses
epidat <- Epi %>% 
	filter(outc != "Other",
				 alcohol_ref != "",
				 energytot != "",
				 gxe == 1) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0),
				 sex = ifelse(sex == "Female", 0, 1),
				 age_ref = as.numeric(age_ref),
				 energytot = as.numeric(energytot),
				 alc_ref = ifelse(alcohol_ref == "No", 0, 1)) %>% 
	dplyr::select(outcome, sex, age_ref, energytot, alc_ref, studyname) %>% 
	filter(complete.cases(.))

exposure <- 'alc_ref'

# CaseControl + Cohort
dat <- epidat %>% 
	filter(studyname %in% caco_cohort)
r_counts <- counts(data = dat, group = 'studyname')
r_glm <- glm_by_group_df(data = dat, group = 'studyname',  exposure = exposure)

png("~/git/FIGI_PostHarmonization/working/figi_postharm_etoh_fp_caco_cohort.png", width = 700, height = 800)
meta_fp(r_counts, r_glm, title = 'Alcohol Intake at Ref Time (Yes/No)\nModel: outc~age_ref+sex+energytot+alc_ref\nCaseControl + Cohort studies')
dev.off()



# CaseControl
# CCFR missing alcohol? 
dat <- epidat %>% 
	filter(studyname %in% caco)
r_counts <- counts(data = dat, group = 'studyname')
r_glm <- glm_by_group_df(data = dat, group = 'studyname',  exposure = exposure)

png("~/git/FIGI_PostHarmonization/working/figi_postharm_etoh_fp_caco.png", width = 700, height = 400)
meta_fp(r_counts, r_glm, title = 'Alcohol Intake at Ref Time (Yes/No)\nModel: outc~age_ref+sex+energytot+alc_ref\nCaseControl studies')
dev.off()



# Cohort
dat <- epidat %>% 
	filter(studyname %in% cohort)
r_counts <- counts(data = dat, group = 'studyname')
r_glm <- glm_by_group_df(data = dat, group = 'studyname',  exposure = exposure)

png("~/git/FIGI_PostHarmonization/working/figi_postharm_etoh_fp_cohort.png", width = 700, height = 700)
meta_fp(r_counts, r_glm, title = 'Alcohol Intake at Ref Time (Yes/No)\nModel: outc~age_ref+sex+energytot+alc_ref\nCohort studies')
dev.off()



#------ GGPLOT ------#
ggdat <- Epi %>% 
	filter(outc != "Other",
				 gxe == 1,
				 studyname %in% caco_cohort) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0),
				 sex = ifelse(sex == "Female", 0, 1),
				 age_ref = as.numeric(age_ref),
				 alc_ref = factor(alcohol_ref, labels = c('NA', 'No', 'Yes'))) %>% 
	group_by(studyname, outcome, alc_ref) %>% 
	summarise(count = n()) %>% 
	mutate(perc=count/sum(count),
				 outcomep = ifelse(outcome == 0, "Controls", "Cases"))

p <- ggplot() + 
	geom_bar(aes(y = perc, x = outcomep, fill = alc_ref), data = ggdat, stat='identity') +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank()) +
	scale_fill_manual(values = c("grey", "red", "blue")) +
	facet_wrap( ~ studyname, ncol = 5)

ggsave("~/git/FIGI_PostHarmonization/working/figi_postharm_etoh_barplot.png", plot = p, device = "png", width = 6, height = 8, scale = 1.5)



#------ MEGA analysis... ------#
dat <- filter(epidat, studyname %in% caco_cohort) %>% 
	mutate(studyname = as.factor(studyname))

x <- glm(outcome ~ sex + age_ref + energytot + alc_ref + studyname, data = dat, family = 'binomial')

exp(coef(summary(x)))
