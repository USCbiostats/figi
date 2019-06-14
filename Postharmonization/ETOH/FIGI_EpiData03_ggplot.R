library(tidyverse)
library(data.table)
library(ggplot2)
rm(list = ls())
load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_09242018.RData")
#Epi <- readRDS("~/git/FIGI_EpiData/Epi.rds")

# Exclude from post-harmonization analysis (Ask Yi about this)
exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPIC", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMC_COSM", "SMHS", "SWHS", "Taiwan", "TRIBE")
cat(unique(sort(Epi$studyname)), sep = '", "')
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "Kentucky", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "REACH_AD", "SMS_AD", "USC_HRT_CRC")
cohort <- c("CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
clintrials <- c("ATBC", "PPS3", "PPS4", "SELECT")
# See Riki email to Yi -- some clinical trials should be treated as cohorts, unless it's the exposure they're trialing
# PPS 3/4 a bit confusing because one of them used aspirin/folate as the trial intervention
# For now, can include as cohort why not 
cohort <- c("ATBC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)

tmp <- Epi %>% 
	filter(outc != "Other",
				 #asp_ref != "",
				 gxe == 1,
				 studyname %in% caco_cohort) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0),
				 age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1),
				 asp_ref = ifelse(asp_ref == "No", 0, 
				 					 ifelse(asp_ref == "", 9, 1)),
				 asp_ref_p = factor(asp_ref, levels = c(9, 0, 1), labels = c("NA", "No", "Yes")),
				 studynamef = factor(studyname))

ffs <- tmp %>% 
	group_by(studynamef, outcome, asp_ref_p) %>% 
	summarise(count = n()) %>% 
	mutate(perc=count/sum(count),
				 outcomep = ifelse(outcome == 0, "Co", "Ca"))

p <- ggplot() + 
	geom_bar(aes(y = perc, x = outcomep, fill = asp_ref_p), data = ffs, stat='identity') +
	theme_bw() + 
	theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank()) +
	scale_fill_manual(values = c("grey", "red", "blue")) +
	facet_wrap( ~ studynamef, ncol = 5)

ggsave("~/git/FIGI_EpiData/working/figi_postharm_etoh_barplot.png", plot = p, device = "png", width = 6, height = 8, scale = 1.5)



#======================================
# let's assemble a table for the post-harmonization analyses
# (This works, copy and pasted to the rmarkdown script)


# fuu <- samplefile %>% 
# 	filter(drop == 0,
# 				 gxe == 1) %>% 
# 	mutate(studynamef = factor(studyname),
# 				 outcome = ifelse(outc == "Control", 0, 
# 				 					 ifelse(adenoma_adv == "Case", 2, 1)))
# 
# # N, case/control with percentages
# a <- fuu %>% 
# 	group_by(studynamef, outc) %>% 
# 	summarise(N = n()) %>% 
# 	mutate(N_freq = paste0(N, " (", 100*round(N / sum(N), 3), "%)")) %>% 
# 	ungroup() %>% 
# 	dplyr::select(-N) %>% 
# 	spread(outc, N_freq) %>% 
# 	dplyr::select(studynamef, Control, Case)
# 
# # crc
# c <- fuu %>% 
# 	filter(outcome == 1) %>% 
# 	group_by(studynamef) %>% 
# 	summarise(crc = n())
# 
# # adv adenoma
# d <- fuu %>% 
# 	filter(outcome == 2) %>% 
# 	group_by(studynamef) %>% 
# 	summarise(adenoma_adv = n())
# 
# 
# finalt <- full_join(a, c, by = "studynamef")
# finalt <- full_join(finalt, d, by = "studynamef")
# finalt[is.na(finalt)] <- ""
