#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FIGI POST HARMONIZATION 10/04/2018
#
# Summary of E Variable - HEIGHT
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
cov <- readRDS("~/git/DATA/FIGI_EpiData_rdata/working/FIGI_GxESet_GxEScanR_height_N79453_20181106.rds") %>% 
  dplyr::select(vcfid, starts_with("PC"))


# let's narrow down samples to the ones included in GxEScanR run
dat <- inner_join(Epi, cov, by = 'vcfid')

table(dat$studyname)
table(dat$studyname, dat$outc)

# study groups
ca <- c("ColoCare_1", "EPICOLON", "NFCCR_1", "NFCCR_2")

caco <- c("ASTERISK","CCFR_1", "CCFR_3", "CCFR_4","Colo2&3","CORSA_1", "CORSA_2", "CRCGEN", "CzechCCS", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "EDRN", "ESTHER_VERDI","HawaiiCCS_AD", "Kentucky", "LCCS", "MECC_1", "MECC_2", "MECC_3","NCCCSI", "NCCCSII","REACH_AD", "SEARCH", "SELECT", "SMC_COSM", "SMS_AD", "USC_HRT_CRC")

cohort <- c("ATBC", "CLUEII", "CPSII_1", "CPSII_2", "EPIC", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD","MCCS_1", "MCCS_2","MEC_1", "MEC_2", "NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD", "PHS","PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD","PPS3", "PPS4",  "VITAL",  "WHI_1", "WHI_2", "WHI_3")
caco_cohort <- c(caco, cohort)
ca_caco_cohort <- c(ca, caco, cohort)


# USA vs Others
# American vs others
usa <- c("CCFR_1", "CCFR_3", "CCFR_4", "Colo2&3", "DALS_1", "DALS_2","EDRN","HawaiiCCS_AD", "Kentucky","NCCCSI", "NCCCSII","REACH_AD", "SELECT", "SMC_COSM", "SMS_AD", "USC_HRT_CRC",  "CLUEII", "CPSII_1", "CPSII_2", "HPFS_1_2", "HPFS_3_AD", "HPFS_4", "HPFS_5_AD", "MEC_1", "MEC_2","NHS_1_2", "NHS_3_AD", "NHS_4", "NHS_5_AD",  "PHS", "PLCO_1", "PLCO_2", "PLCO_3", "PLCO_4_AD","PPS3", "PPS4",   "VITAL",  "WHI_1", "WHI_2", "WHI_3")

other <- c("ASTERISK","CORSA_1", "CORSA_2", "CRCGEN", "CzechCCS", "DACHS_1", "DACHS_2", "DACHS_3",  "ESTHER_VERDI", "LCCS", "MECC_1", "MECC_2", "MECC_3","SEARCH","ATBC", "EPIC", "MCCS_1", "MCCS_2")


epidat <- dat %>% 
  mutate(studyname = ifelse(studyname %in% c("HPFS_1", "HPFS_2"), "HPFS_1_2", studyname),
         studyname = ifelse(studyname %in% c("NHS_1", "NHS_2"), "NHS_1_2", studyname)) %>% 
  mutate(outcome = ifelse(outc == "Control", 0, 
                          ifelse(adenoma_adv == "Case", 2, 1))) # creating one variable for crc + adv adenoma outcomes




#------ OUTCOME TABLE ------
dat <- filter(epidat, studyname %in% caco_cohort)
counts_df <- counts_control_case_adv(data = dat, group = 'studyname')
saveRDS(counts_df, file = "~/Dropbox/code/FIGI_PostHarmonization/working/height_counts.rds")


#------ OUTCOME PLOT ------
# plot percentage case/control/adv
# factor level order matters - sorting by # of controls
# tmp <- sort(unique(epidat$studyname), decreasing = T)

ggdat <- epidat %>% 
  filter(studyname %in% ca_caco_cohort) %>% 
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
  geom_bar(aes(fill = outcome_f), stat='identity', alpha = 0.6, width = 0.9) +
  coord_flip() + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c("grey", "red", "blue"))

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_outcome_barplot.png", plot = p, device = "png", width = 6, height = 6, scale = 1.2)


#------ META ANALYSIS ------
# model outc to match Yi shinyapp

metadat <- epidat %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         sex = ifelse(sex == "Female", 0, 1),
         age_ref = as.numeric(age_ref),
         heightcm = as.numeric(heightcm)) %>% 
  dplyr::select(outcome, heightcm, sex, age_ref, starts_with('PC'), studyname)



#------ Helper Function ------
# data <- metadat
# exposure = 'heightcm'
# covars = c('age_ref', 'sex', paste0(rep("PC", 10), seq(1,10)))
# covars = c('age_ref', 'sex')
# study_subset = caco_cohort

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

# Height
meta_fp_helper(data = metadat,
               exposure = 'heightcm',
               covars = c('age_ref', 'sex', paste0(rep("PC", 10), seq(1,10))), 
               study_subset = caco_cohort, 
               png_name = "~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_fp_caco_cohort.png",
               width = 700, 
               height = 1000, 
               title = 'Height (cm)\nModel: outc~age_ref+sex+PC1-10+height\nCaseControl + Cohort studies')



# GGPLOT ----
ggdat <- epidat %>% 
  filter(gxe == 1, # gxe == 1 takes care of outc == "Other"
         studyname %in% caco_cohort) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         sex = factor(sex, labels = c("Female", "Male")),
         heightcm = as.numeric(heightcm),
         outcomep = ifelse(outcome == 0, 'Control', 'Case')) %>% 
  filter(!is.na(heightcm)) 

ggdat2 <- ggdat %>% 
  group_by(studyname) %>% 
  summarise(heightcm_median = median(heightcm),
            heightcm_mean = mean(heightcm)) %>% 
  ungroup() %>% arrange(heightcm_mean)

sort_h <- c(rev(unique(ggdat2$studyname)))

test <- ggdat %>% 
  mutate(studyname_f = factor(studyname, levels = sort_h))

p <- ggplot(aes(y = heightcm, x = outcomep, fill = outcomep), data = test) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_boxplot.png", plot = p, device = "png", width = 5, height = 8, scale = 1.2)


# case-control only
test2 <- filter(test, studyname %in% caco)
p <- ggplot(aes(y = heightcm, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("Case Control") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_boxplot_caco.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)


# Cohort
test2 <- filter(test, studyname %in% cohort)
p <- ggplot(aes(y = heightcm, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_boxplot_cohort.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)



# USA
test2 <- filter(test, studyname %in% usa)
p <- ggplot(aes(y = heightcm, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("USA") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_boxplot_usa.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)


# Other
test2 <- filter(test, studyname %in% other)
p <- ggplot(aes(y = heightcm, x = outcomep, fill = outcomep), data = test2) + 
  geom_boxplot(aes(x = studyname_f), outlier.shape = 20, alpha = 0.4) + 
  coord_flip() +
  ggtitle("Other") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'bottom')
ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_boxplot_other.png", plot = p, device = "png", width = 4, height = 5, scale = 1.5)


# GGPLOT Folate - controls only, visualize study_design ----
# sort table by median folate_tot levels (grouped by study)
ggdat2 <- ggdat %>% 
  filter(outcome == 0) %>% 
  group_by(studyname) %>% 
  summarise(folate_med = median(folate_tot)) %>% 
  arrange(folate_med) %>% 
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

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_folate_boxplot_control_studydesign.png", plot = p, device = "png", width = 5, height = 6, scale = 1.2)





# Heterogeneity ----
# inverse variance weighted p value.. 
exposure = 'heightcm'
hetero <- metadat %>% 
  filter(studyname %in% caco_cohort)

covars = c('age_ref', 'sex', paste0(rep("PC", 10), seq(1,10)))
formula_txt <- paste0("outcome ~ ", paste(covars, collapse = "+"))


# by study design (no heterogeneity)
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(study_design = ifelse(studyname %in% caco, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/height_heterogeneity_studydesign.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ study_design, data = r_glm, weights = w)
summary(x)


# by study country (no heterogeneity)
r_glm <- glm_by_group_df(data = hetero, group = 'studyname', exposure = exposure, formula_txt = formula_txt) %>% 
  mutate(country = ifelse(studyname %in% other, 0, 1))
saveRDS(r_glm, file = "~/Dropbox/code/FIGI_PostHarmonization/working/height_heterogeneity_usa.rds")
w <- (r_glm$std.error)^2
x <- lm(estimate ~ country, data = r_glm, weights = w)
summary(x)




#------ Violin ------
# remember test is defined above!
test2 <- filter(test, studyname %in% caco_cohort)
p1 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("CaseControl + Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% usa)
p2 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("USA") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% other)
p3 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("Other") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_violin_country.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)


#--------------
test2 <- filter(test, studyname %in% caco_cohort)
p1 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("CaseControl + Cohort") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% caco)
p2 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("Case Control") +
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

test2 <- filter(test, studyname %in% cohort)
p3 <- ggplot(aes(y = heightcm, x = sex, fill = sex), data = test2) + 
  geom_violin() +
  geom_boxplot(width = 0.05) + 
  ylim(130, 215) +
  ggtitle("Cohort") + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = 'none')

library(gridExtra)
pp <- grid.arrange(p1, p2, p3, nrow = 1)

ggsave("~/Dropbox/code/FIGI_PostHarmonization/working/figi_postharm_height_violin_studydesign.png", plot = pp, device = "png", width = 6, height = 5, scale = 1.5)





#------ MEGA analysis... ------#
dat <- filter(epidat, studyname %in% caco_cohort) %>% 
  mutate(studyname = as.factor(studyname))

x <- glm(outcome ~ sex + age_ref + heightcm + studyname, data = dat, family = 'binomial')

exp(coef(summary(x)))


