#-----------------------------------------------------------
# FIGI POST HARMONIZATION Analysis 10/04/2018
#
# 1) Counts for cases/controls etc
# 2) Summary of E Variable - plots
# 3) Meta Analyses and study heterogeneity
# 4) Forestplots
#-----------------------------------------------------------
library(tidyverse)
library(broom)
library(data.table)
library(rmeta)

rm(list = ls())
load("~/git/FIGI_EpiData/FIGI_Genotype_Epi_09242018.RData")

cat(unique(sort(Epi$studyname)), sep = '", "')
#exclude <- c("ASTERISK", "CCFR_2", "CGN", "ColoCare_1", "CORSA_1", "CORSA_2", "Czech_CCS", "EPIC", "EPICOLON", "ESTHER_VERDI", "FIRE3", "GALEON", "HawaiiCCS_AD", "HCES-CRC", "HispanicCCS", "HPFS_1", "HPFS_2", "LCCS", "MAVERICC", "MOFFITT", "MSKCC", "NFCCR_1", "NGCCS", "NHSII", "PMH-CCFR", "PURIFICAR", "SLRCCS", "SMC_COSM", "SMHS", "SWHS", "Taiwan", "TRIBE")
caco <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "Kentucky", "MECC_1", "MECC_2", "MECC_3", "NCCCSI", "NCCCSII", "REACH_AD", "SMS_AD", "USC_HRT_CRC")
cohort <- c("ATBC", "SELECT", "PPS3", "PPS4", "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1", "HPFS_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MCCS_1", "MCCS_2", "MEC_1", "MEC_2", "NHS_1", "NHS_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD", "UKB_1", "VITAL", "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)

#---------------------------------------------------------
# META-ANALYSIS ALCOHOL_REF
#---------------------------------------------------------
Epiw <- Epi %>% 
	filter(gxe == 1) %>% 
	mutate(age_ref = as.numeric(age_ref),
				 sex = ifelse(sex == "Female", 0, 1),
				 studynamef = factor(studyname))

all_studies_crc <- Epiw %>% 
	filter(outc != "Other",
				 alcohol_ref != "",
				 studyname %in% caco_cohort) %>% 
	mutate(outcome = ifelse(outc == "Case", 1, 0), 
				 alcohol_ref = ifelse(alcohol_ref == "No", 0, 1))
	
caco_crc <- all_studies_crc %>% 
	filter(studyname %in% caco)

cohort_crc <- all_studies_crc %>% 
	filter(studyname %in% cohort)

#------ Run GLM by studyname ------#
by_studyname <- all_studies_crc %>% 
	group_by(studyname)
results_main <- do(by_studyname, tidy(glm(outcome ~ age_ref + sex + alcohol_ref, data = ., family = 'binomial')))
results_ci <- do(by_studyname, confint_tidy(glm(outcome ~ age_ref + sex + alcohol_ref, data = ., family = 'binomial'))) %>% 	filter(!is.na(conf.low))

results <- bind_cols(results_main, results_ci) %>% 
	ungroup %>% 
	dplyr::filter(term == "alcohol_ref") %>% 
	mutate(design = factor(ifelse(studyname %in% caco, 0, 1)))

# Heterogeneity and meta analyses
w <- (results$std.error)^2
x <- lm(estimate ~ design, data = results, weights = w)
coef(summary(x))

mm <- meta.summaries(results$estimate, results$std.error, method = 'random')
sum_mm <- c(mm[[3]], mm[[3]]-(mm[[4]]*1.96), mm[[3]]+(mm[[4]]*1.96))

# forest plot df
fp <- dplyr::select(results, estimate, conf.low, conf.high)
fp <- rbind(rep(NA, 3), fp, rep(NA, 3), sum_mm)

wtf <- all_studies_crc %>% 
	group_by(outcome, studyname) %>% 
	summarize(count = n()) %>% 
	spread(key = outcome, value = count) %>% 
	rename(Control = `0`, Case = `1`)

lab_summary <- c("Summary", nrow(all_studies_crc), sum(wtf$Case), round(mm[[3]], 2), round(mm[[4]], 3), round(exp(mm[[3]]), 2), round(mm[[5]][2], 4))
lab <- inner_join(wtf, results, by = 'studyname') %>% 
	mutate(OR = round(exp(estimate),2),
				 N = Case+Control,
				 Pvalue = round(p.value, 4),
				 SE = round(std.error, 3),
				 beta = round(estimate, 2)) %>% 
	dplyr::select(studyname, N, Case, beta, SE, OR, Pvalue) 

lab <- rbind(colnames(lab), lab, rep(NA, 7), lab_summary)


png("~/git/FIGI_EpiData/working/fp1_alcohol_ref.png", width = 800, height = 1000)
forestplot(label=lab, as.numeric(fp$estimate), as.numeric(fp$conf.low), as.numeric(fp$conf.high), is.summary = c(T, rep(F, nrow(fp)-2), T), col=meta.colors(box="royalblue", line="darkblue", zero='gray0', summary="royalblue"), clip=c(-1.5,1.5), xlog=T, xlab='logOR', zero=0,  title(main = "Alcohol (yes/no) at ref time\nModel: outc~age_ref+sex+alcohol_ref"))
dev.off()


png("~/git/FIGI_EpiData/working/fp1_alcohol_ref.png", width = 800, height = 1000)
forestplot(label=lab[,c(1,2,3,6,7)], as.numeric(fp$estimate), as.numeric(fp$conf.low), as.numeric(fp$conf.high), clip=c(-1.5,1.5), xlog = T)
dev.off()

