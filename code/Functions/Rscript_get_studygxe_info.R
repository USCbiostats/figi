#=============================================================================#
# Create data.frame that contains study_gxe, study design, continent, country (?)
#=============================================================================#

# input dataset
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")


gxe <- figi %>%
  filter(drop == 0 & gxe == 1) %>%
  mutate(study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe))

study_gxe <- sort(unique(gxe$study_gxe))
study_gxe

# [1] "ASTERISK"       "ATBC"           "CCFR_1"         "CCFR_3"         "CCFR_4"         "CLUEII"
# [7] "Colo23"         "ColoCare_1"     "CORSA_1"        "CORSA_2"        "CPSII_1"        "CPSII_2"
# [13] "CRCGEN"         "CzechCCS"       "DACHS_1"        "DACHS_2"        "DACHS_3"        "DALS_1"
# [19] "DALS_2"         "EDRN"           "EPIC"           "EPICOLON"       "ESTHER_VERDI"   "GALEON"
# [25] "HawaiiCCS_AD"   "HPFS_1_2"       "HPFS_3_AD"      "HPFS_4"         "HPFS_5_AD"      "Kentucky"
# [31] "LCCS"           "MCCS_1"         "MCCS_2"         "MEC_1"          "MEC_2"          "MECC_1"
# [37] "MECC_2"         "MECC_3"         "MOFFITT"        "NCCCSI"         "NCCCSII"        "NFCCR_1"
# [43] "NFCCR_2"        "NGCCS"          "NHS_1_2"        "NHS_3_AD"       "NHS_4"          "NHS_5_AD"
# [49] "PHS"            "PLCO_1_Rematch" "PLCO_2"         "PLCO_3"         "PLCO_4_AD"      "PPS3"
# [55] "PPS4"           "REACH_AD"       "SEARCH"         "SELECT"         "SMC_COSM"       "SMS_AD"
# [61] "UKB_1"          "USC_HRT_CRC"    "VITAL"          "WHI_1"          "WHI_2"          "WHI_3"


study_design <- c("Case-Control", "Cohort", "Case-Control", "Case-Control", "Case-Control", "Cohort",
                  "Case-Control", "Case-Series", "Case-Control", "Case-Control", "Cohort", "Cohort",
                  "Case-Control", "Case-Control", "Case-Control", "Case-Control", "Case-Control", "Case-Control",
                  "Case-Control", "Case-Series", "Cohort", "Case-Control", "Case-Control", "Case-Control",
                  "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort", "Case-Control",
                  "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort", "Case-Control",
                  "Case-Control", "Case-Control", "Case-Series", "Case-Control", "Case-Control", "Case-Control",
                  "Case-Control", "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort",
                  "Cohort", "Cohort", "Cohort", "Cohort", "Cohort", "Clinical Trial",
                  "Clinical Trial", "Case-Control", "Case-Control", "Clinical Trial", "Cohort", "Cohort",
                  "Cohort", "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort")

study_design <- c("Case-Control", "Cohort", "Case-Control", "Case-Control", "Case-Control", "Cohort",
									"Case-Control", "Case-Series", "Case-Control", "Case-Control", "Cohort", "Cohort",
									"Case-Control", "Case-Control", "Case-Control", "Case-Control", "Case-Control", "Case-Control",
									"Case-Control", "Case-Control", "Cohort", "Case-Control", "Case-Control", "Case-Control",
									"Case-Control", "Cohort", "Cohort", "Cohort", "Cohort", "Case-Control",
									"Case-Control", "Cohort", "Cohort", "Cohort", "Cohort", "Case-Control",
									"Case-Control", "Case-Control", "Case-Series", "Case-Control", "Case-Control", "Case-Control",
									"Case-Control", "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort",
									"Cohort", "Cohort", "Cohort", "Cohort", "Cohort", "Cohort",
									"Cohort", "Case-Control", "Case-Control", "Cohort", "Cohort", "Cohort",
									"Cohort", "Case-Control", "Cohort", "Cohort", "Cohort", "Cohort")


study_continent <- c("Europe", "Europe", "Australia/America", "America", "America", "America",
                     "America", "Europe/America", "Europe", "Europe", "America", "America",
                     "Europe", "Europe", "Europe", "Europe", "Europe", "America",
                     "America", "America", "Europe", "Europe", "Europe", "Europe",
                     "America", "America", "America", "America", "America", "America",
                     "Europe", "Australia", "Australia", "America", "America", "Israel",
                     "Israel", "Israel", "America", "America", "America", "America",
                     "America", "Europe", "America", "America", "America", "America",
                     "America", "America", "America", "America", "America", "America",
                     "America", "America", "Europe", "America", "Europe", "America",
                     "Europe", "America", "America", "America", "America", "America")

study_info <- data.frame(study_gxe, study_design, study_continent)
saveRDS(study_info, "~/data/Annotations/FIGI_studygxe_info.rds", version = 2)


# figi_gxe <- figi %>%
#   filter(drop == 0 & gxe == 1) %>%
#   inner_join(pc30k, by = c('vcfid' = 'IID')) %>%
#   mutate(outcome = ifelse(outc == "Case", 1, 0),
#          age_ref_imp = as.numeric(age_ref_imp),
#          sex = ifelse(sex == "Female", 0, 1),
#          study_gxe = ifelse(study_gxe == "Colo2&3", "Colo23", study_gxe),
#          aspirin = ifelse(aspirin == "Yes", 1, ifelse(asp_ref == "No", 0, NA)),
#          country = factor(ifelse(study_site %in% c('Alabama', 'Kentucky', 'Los Angeles', 'Seattle', 'Hawaii', 'Mayo Foundation', 'Midwest', 'South', 'West', 'Northeast'), 'USA',
#                           ifelse(study_site %in% c('Barcelona', 'Leon', 'Santiago', 'Asturias'), 'Spain',
#                           ifelse(study_site %in% c('MECC'), 'Israel',
#                           ifelse(study_site %in% c('Melbourne', 'Australia'), 'Australia',
#                           ifelse(study_site %in% c('Heidelberg'), 'Germany',
#                           ifelse(study_site %in% c('888'), 'Italy',
#                           ifelse(study_site %in% c('Central Sweden'), 'Sweden',
#                           ifelse(study_site %in% c('Ontario', 'Newfoundland'), 'Canada',
#                           ifelse(study == 'ATBC', 'Finland',
#                           ifelse(study == 'CzechCCS', 'Czech Republic',
#                           ifelse(study == 'EPICOLON', 'Spain',
#                           ifelse(study %in% c('ASTERISK'), 'France',
#                           ifelse(study %in% c('FIRE3', 'DACHS', 'ESTHER_VERDI', 'NGCCS'), 'Germany',
#                           ifelse(study %in% c('SEARCH', 'UKB', 'LCCS'), 'UK',
#                           ifelse(study %in% c('DALS', 'EDRN', 'MOFFITT', 'VITAL', 'CLUEII', 'Colo23', 'NCCCSI', 'NCCCSII', 'PLCO', 'PPS3', 'PPS4', 'SELECT', 'NHS', 'PHS', 'REACH', 'USC_HRT_CRC', 'HawaiiCCS', 'HPFS', 'SMS'), 'USA',
#                           ifelse(study == 'CORSA', 'Austria', NA))))))))))))))))))


# wtf <- data.frame(table(figi_gxe$study_gxe, figi_gxe$country)) %>%
#   filter(Freq != 0) %>%
#   mutate(continent = recode(Var2,
#                             "USA" = "America",
#                             "Canada" = "America",
#                             "Australia" = "Australia",
#                             "Austria" = "Europe",
#                             "Czech Republic" = "Europe",
#                             "Finland" = "Europe",
#                             "France" = "Europe",
#                             "Germany" = "Europe",
#                             "Spain" = "Europe",
#                             "Sweden" = "Europe",
#                             "UK" = "Europe",
#                             "Israel" = "Israel"))
#
# wtf2 <- wtf %>% dplyr::select(Var1, continent) %>% arrange(Var1) %>%
#   filter(!duplicated(.))
# as.character(wtf2$continent)
#
# cat(paste(shQuote(wtf2$continent, type="cmd"), collapse=", ")) # nice


