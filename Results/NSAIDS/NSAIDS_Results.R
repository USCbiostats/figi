#================================================================================#
# FIGI NSAID GxE
# (try to make it generalizable for other Es)
#
# Create R objects for:
# 1) Post-Harmonization 
# 2) GWIS Results
# 3) Presentations
#================================================================================#
library(tidyverse)
library(data.table)
library(broom)
library(metafor)
library(rmeta)
library(reshape)
library(rlang)
library(tableone)
library(sjPlot)
rm(list = ls())
source("~/Dropbox/FIGI/code/Functions/postharmonization_GWISresults_useful_functions.R")

# ------ Read File (minimal processing) ------ 
# use study_gxe as grouping variable as shortcut
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190313.RData")
gxe_set <- figi %>% 
	filter(drop == 0 & gxe == 1) %>%  # GxE set - N = 102792 as of 12/13/2018
	mutate(aspirinyrs_b = ifelse(aspirinyrs > 0, 1, 0),
				 nsaidsyrs_b = ifelse(nsaidsyrs > 0, 1, 0),
				 aspyrs_ref_b = ifelse(aspyrs_ref > 0, 1, 0)) # checking availability of different nsaid vars

# ------ Recreate asp_ref variable to ensure you understand coding ------ #
z <- dplyr::select(gxe_set, aspirin, nsaids, asp_ref, asp_ref2) %>% 
	mutate(asp_ref_test = ifelse(aspirin == "Yes" | nsaids == "Yes", "Yes", 
												ifelse(aspirin == "No" & (nsaids == "" | nsaids == "No"), "No", 
												ifelse(nsaids == "No" & (aspirin == "" | aspirin == "No"), "No", ""))))
table(z$asp_ref_test, z$asp_ref) # good


# ------ Counts Table ------
# Make one for each aspirin/nsaid variable for the Rmd doc
# by studyname, case/control, total N
# these will need some tidying before plotting on rmarkdown doc

counts_asp_ref <- getCounts(gxe_set, study_gxe, asp_ref, outc)
counts_asp_ref2 <- getCounts(gxe_set, study_gxe, asp_ref2, outc)
counts_aspirin <- getCounts(gxe_set, study_gxe, aspirin, outc)
counts_nsaids <- getCounts(gxe_set, study_gxe, nsaids, outc)

saveRDS(counts_asp_ref, file = 'rds/counts_asp_ref.rds')
saveRDS(counts_asp_ref2, file = 'rds/counts_asp_ref2.rds')
saveRDS(counts_aspirin, file = 'rds/counts_aspirin.rds')
saveRDS(counts_nsaids, file = 'rds/counts_nsaids.rds')


# ------ Bar Plots ------
# format data with proper factors/labels
# make sure the nsaid variables are factors with order "NA", "No", "Yes"

df <- gxe_set %>% 
	mutate(outc = factor(outc, labels = c("Case", "Control")),
				 asp_ref = factor(asp_ref, exclude = NULL, labels = c("NA", "No", "Yes" )),
				 asp_ref2 =  factor(asp_ref2, exclude = NULL, labels = c("NA", "No", "Yes" )),
				 aspirin = factor(aspirin, exclude = NULL, labels = c("NA", "No", "Yes" )),
				 nsaids = factor(nsaids, exclude = NULL, labels = c("NA", "No", "Yes" )))

# saving ggplot object as rds
saveRDS(getPlots(df, study_gxe, asp_ref, outc), file = "figures/barplot_asp_ref.rds") 
saveRDS(getPlots(df, study_gxe, asp_ref2, outc), file = "figures/barplot_asp_ref2.rds") 
saveRDS(getPlots(df, study_gxe, aspirin, outc), file = "figures/barplot_aspirin.rds") 
saveRDS(getPlots(df, study_gxe, nsaids, outc), file = "figures/barplot_nsaids.rds") 


# ------ Bar Plots - 4 pages ------
# Split Bar Plots into two pages, to fit in presentation
# asp_ref only for now
# includes all studies to highlight case-only/missingness

df <- gxe_set %>% 
  mutate(outc = factor(outc, labels = c("Case", "Control")),
         asp_ref = factor(asp_ref, exclude = NULL, labels = c("NA", "No", "Yes" )),
         asp_ref2 =  factor(asp_ref2, exclude = NULL, labels = c("NA", "No", "Yes" )),
         aspirin = factor(aspirin, exclude = NULL, labels = c("NA", "No", "Yes" )),
         nsaids = factor(nsaids, exclude = NULL, labels = c("NA", "No", "Yes" )))
table(df$study_gxe)

study_list <- sort(unique(df$study_gxe))
length(study_list)

# df_a <- filter(df, study_gxe %in% study_list[1:20])
# df_b <- filter(df, study_gxe %in% study_list[21:40])
# df_c <- filter(df, study_gxe %in% study_list[41:60])
# df_d <- filter(df, study_gxe %in% study_list[61:66])
# 
# saveRDS(getPlots(df_a, study_gxe, asp_ref, outc), file = "./figures/barplot_asp_ref_part1of4.rds")
# saveRDS(getPlots(df_b, study_gxe, asp_ref, outc), file = "./figures/barplot_asp_ref_part2of4.rds")
# saveRDS(getPlots(df_c, study_gxe, asp_ref, outc), file = "./figures/barplot_asp_ref_part3of4.rds")
# saveRDS(getPlots(df_d, study_gxe, asp_ref, outc), file = "./figures/barplot_asp_ref_part4of4.rds")

# also save as png... for powerpoint
# (for this one, modified getPlots to use 9 columns rather than 5
df_a <- filter(df, study_gxe %in% study_list[1:35])
df_b <- filter(df, study_gxe %in% study_list[36:66])

ggsave("./figures/barplot_asp_ref_powerpoint_part1of2.png", getPlots(df_a, study_gxe, asp_ref, outc), device = "png", width = 14, units = 'in')
ggsave("./figures/barplot_asp_ref_powerpoint_part2of2.png", getPlots(df_b, study_gxe, asp_ref, outc), device = "png", width = 14, units = 'in')


# ------ Meta Analysis + Forest Plots ------
# again, format data appropriately.. keep outcome coded as "Case" / "Control"
# remove studies that are case only, control only (if any)
# note the 0/1 coding for outcome and other variables (outc_01), rather than text
# save the getCounts + glm_by_group_df objects, process in rmarkdown

exclude_study <- c("GALEON", "MOFFITT", "NFCCR_1", "NGCCS", "ColoCare_1", "ESTHER_VERDI")

df <- gxe_set %>% 
	filter(!study_gxe %in% exclude_study) %>% 
	mutate(outc_01 = ifelse(outc == "Case", 1, 0),
				 age_ref_imp = as.numeric(age_ref_imp),
				 sex = ifelse(sex == "Female", 1, ifelse(sex == "Male", 0, NA)),
				 asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA)),
				 asp_ref2 = ifelse(asp_ref2 == "Yes", 1, ifelse(asp_ref2 == "No", 0, NA)),
				 aspirin = ifelse(aspirin == "Yes", 1, ifelse(aspirin == "No", 0, NA)),
				 nsaids = ifelse(nsaids == "Yes", 1, ifelse(nsaids == "No", 0, NA)))
	# dplyr::select(outc, outc_01, asp_ref, age_ref_imp, sex, study_gxe) %>% 
	# filter(complete.cases(.)) # doing this in later step

# convenience wrapper function
# PAY ATTENTION TO GROUP, OUTCOME, AND COVARIATES WHEN USING THIS WRAPPER
# also output location of figure
quick_wrap <- function(exposure, title) {
	lol <- enexpr(exposure)
	x <- getCounts_outc(data = df, group = study_gxe, exposure = !! lol, outcome = outc, age_ref_imp, sex)
	y <- glm_by_group_df(data = df, group = study_gxe, exposure = !! lol, outcome = outc_01, age_ref_imp, sex)
	png(file = paste0("./figures/forestplot_", quo_name(lol), ".png"), width = 8, height = 13.5, units = 'in', res = 200)
	z <- meta_fp(x, y, group = study_gxe, title = title)
	dev.off()
}

title <- "Aspirin/NSAID use @ Ref Time \nModel: outc ~ age_ref_imp + sex + asp_ref"
quick_wrap(asp_ref, title)

title <- "Aspirin use @ Ref Time \nModel: outc~age_ref_imp+sex+aspirin"
quick_wrap(aspirin, title)

# title <- "NSAID use @ Ref Time \nModel: outc~age_ref_imp+sex+nsaids"
# quick_wrap(nsaids, title)



# ------ Meta Analysis Heterogeneity ------
source("~/Dropbox/FIGI/code/Functions/postharmonization_GWISresults_useful_functions.R")
cov <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_GxESet_asp_ref_sex_age_pc10_studygxe_72145_GLM.rds")
names(cov); table(cov$study_gxe)

# by study design
studylist <- sort(unique(cov$study_gxe))
casecontrol <- studylist[c(2:4,6,9:15,20,25:30,42,44,46)]
cohort <- studylist[!studylist %in% casecontrol]

model <- glm_by_group_df(data = cov, group = study_gxe, exposure = asp_ref, outcome = outcome, age_ref_imp, sex) %>%
  mutate(study_design = ifelse(study_gxe %in% cohort, 0, 1))
inv_var_weight <- model$std.error^2
model_meta <- lm(estimate ~ study_design, data = model, weights = inv_var_weight)
summary(model_meta)
stargazer(model_meta, type = 'html', title = "Study Design", single.row = T)


# by study country (america vs non_american)
studylist <- sort(unique(cov$study_gxe))
american <- studylist[c(2:8,13:20,23:24,28:44,46:50)]
nonamerican <- studylist[!studylist %in% american]

model <- glm_by_group_df(data = cov, group = study_gxe, exposure = asp_ref, outcome = outcome, age_ref_imp, sex) %>%
  mutate(american_nonAmerican = ifelse(study_gxe %in% american, 0, 1))
inv_var_weight <- model$std.error^2
model_meta <- lm(estimate ~ american_nonAmerican, data = model, weights = inv_var_weight)
summary(model_meta)
stargazer(model_meta, type = 'html', title = "Study Country (American vs Non-American)", single.row = T)



# ------ table one ------
# or... table1. not sure which one i should use.. 
# use df created in the meta analysis step 
# for this package, specify categorical variables
# (dput puts variables in vectors)
df <- gxe_set %>% 
  mutate(outc_01 = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 1, ifelse(sex == "Male", 0, NA)),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA)),
         asp_ref2 = ifelse(asp_ref2 == "Yes", 1, ifelse(asp_ref2 == "No", 0, NA)),
         aspirin = ifelse(aspirin == "Yes", 1, ifelse(aspirin == "No", 0, NA)),
         nsaids = ifelse(nsaids == "Yes", 1, ifelse(nsaids == "No", 0, NA)))

# dput(names(df))
# myVars <- c("outc", "age_ref_imp", "sex", "asp_ref")
# catVars <- c("sex", "asp_ref")
# tab2 <- CreateTableOne(vars = myVars, data = df, factorVars = catVars)
# tab2


# Table1
# (must format table more explicitly)
# (output a smaller table, process rds in rmarkdown)
all_vars <- c("outc", "age_ref_imp", "sex", "asp_ref", "study_gxe")
numeric_vars <- c("age_ref_imp")

df <- gxe_set %>% 
  dplyr::select(all_vars) %>% 
  mutate_at(numeric_vars, as.numeric) %>%
  mutate(outcome = factor(outc, labels = c("Control", "Case")),
         asp_ref = factor(asp_ref, labels = c("No", "Yes")),
         sex = factor(sex, labels = c('Female', 'Male'))) 

label(df$age_ref_imp) <- "Age"
label(df$sex) <- "Sex"
label(df$asp_ref) <- "Regular Aspirin/NSAID use"
label(df$outcome) <- "Outcome (CRC)"

saveRDS(df, file = "rds/df_table1.rds")



# ------ MEGA analysis ------

df <- gxe_set %>% 
  filter(asp_ref != "",
         !study_gxe %in% exclude_study) %>% 
  dplyr::select(all_vars) %>%
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         asp_ref = factor(asp_ref, labels = c("asp_ref_No", "asp_ref_Yes")),
         sex = factor(sex, labels = c('Female', 'Male'))) 

glm1 <- glm(outcome ~ age_ref_imp + sex + study_gxe + asp_ref, data = df, family = 'binomial')
saveRDS(glm1, file = "rds/mega_analysis_glm1.rds")
# tab_model(glm1)
















