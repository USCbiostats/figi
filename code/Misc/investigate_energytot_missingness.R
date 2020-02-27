#=============================================================================#
# Explore the energytot.missing variable 
#
# adjustment for energytot.missing is causing problems with gxescan
# because of collinearity with study_gxe
# 
# let's get a sense of the numbers of studies that would be dropped in each E analysis
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
rm(list = ls())
figi_gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)


# redmeatqc2 table with and without dropped those studies with missing values (with total column)
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  mutate(redmeatqc2 = fct_explicit_na(redmeatqc2, na_level = "NA"),
         redmeatqc2 = fct_relevel(redmeatqc2, "NA", "1", "2", "3", "4"))

tmp <- figi_gxe %>% 
  group_by(study_gxe, outc) %>% 
  count(redmeatqc2)

totals <- ungroup(test) %>%
  group_by(outc, redmeatqc2) %>% 
  summarise(n = sum(n))
totals$study_gxe = "Total"

out <- bind_rows(tmp, totals) %>% 
  mutate(pct = round( (100*n) / sum(n), 1),
         out = paste0(n, " (", pct, ")"),
         redmeatqc2 = fct_relevel(redmeatqc2, "1", "2", "3", "4", "NA")) %>%
  dplyr::select(-n, -pct) %>% 
  arrange(outc, redmeatqc2) %>% 
  pivot_wider(names_from = c(outc, redmeatqc2), values_from = out, ) %>% 
  arrange(study_gxe)

out1 <- rbind(out[out$study_gxe != "Total", ], out[out$study_gxe == "Total", ]) %>% 
  filter(!is.na(Case_2))



# same as above but dropping studies with missing energy information altogether
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  mutate(redmeatqc2 = fct_explicit_na(redmeatqc2, na_level = "NA"),
         redmeatqc2 = fct_relevel(redmeatqc2, "NA", "1", "2", "3", "4")) %>% 
  filter(energytot.missing == 0)

tmp <- figi_gxe %>% 
  group_by(study_gxe, outc) %>% 
  count(redmeatqc2)

totals <- ungroup(test) %>%
  group_by(outc, redmeatqc2) %>% 
  summarise(n = sum(n))
totals$study_gxe = "Total"

out <- bind_rows(tmp, totals) %>% 
  mutate(pct = round( (100*n) / sum(n), 1),
         out = paste0(n, " (", pct, ")"),
         redmeatqc2 = fct_relevel(redmeatqc2, "1", "2", "3", "4", "NA")) %>%
  dplyr::select(-n, -pct) %>% 
  arrange(outc, redmeatqc2) %>% 
  pivot_wider(names_from = c(outc, redmeatqc2), values_from = out, ) %>% 
  arrange(study_gxe)

out2 <- rbind(out[out$study_gxe != "Total", ], out[out$study_gxe == "Total", ]) %>% 
  filter(!is.na(Case_2))




# generate a simpler table - case/control counts before and after removing studies with missing energy information
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  mutate(redmeatqc2 = fct_explicit_na(redmeatqc2, na_level = "NA"),
         redmeatqc2 = fct_relevel(redmeatqc2, "NA", "1", "2", "3", "4")) %>% 
  filter(redmeatqc2 != "NA")


tmp <- figi_gxe %>% 
  group_by(study_gxe) %>% 
  count(outc)

totals <- ungroup(tmp) %>%
  group_by(outc) %>% 
  summarise(n = sum(n))
totals$study_gxe = "Total"

out <- bind_rows(tmp, totals) %>% 
  mutate(pct = round( (100*n) / sum(n), 1),
         out = paste0(n, " (", pct, ")")) %>% 
  dplyr::select(-n, -pct) %>% 
  arrange(outc) %>% 
  pivot_wider(names_from = c(outc), values_from = out )



# generate a simpler table - case/control counts before and after removing studies with missing energy information
# but use the subset that you were initially submitting
x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  filter(vcfid %in% x$vcfid)



tmp <- figi_gxe %>% 
  group_by(study_gxe) %>% 
  count(outc)

totals <- ungroup(tmp) %>%
  group_by(outc) %>% 
  summarise(n = sum(n))
totals$study_gxe = "Total"

out <- bind_rows(tmp, totals) %>% 
  mutate(pct = round( (100*n) / sum(n), 1),
         out = paste0(n, " (", pct, ")")) %>% 
  dplyr::select(-n, -pct) %>% 
  arrange(outc) %>% 
  pivot_wider(names_from = c(outc), values_from = out )





x <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_redmeatqc2_basic_covars_glm.rds")
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  filter(vcfid %in% x$vcfid,
         energytot.missing == 0, 
         energytot_imp != 0)

tmp <- figi_gxe %>% 
  group_by(study_gxe) %>% 
  count(outc)

totals <- ungroup(tmp) %>%
  group_by(outc) %>% 
  summarise(n = sum(n))
totals$study_gxe = "Total"

out2 <- bind_rows(tmp, totals) %>% 
  mutate(pct = round( (100*n) / sum(n), 1),
         out = paste0(n, " (", pct, ")")) %>% 
  dplyr::select(-n, -pct) %>% 
  arrange(outc) %>% 
  pivot_wider(names_from = c(outc), values_from = out )
names(out2) <- c("study_gxe", "Case (no energy)", "Control (no energy)")

wtf <- full_join(out, out2, 'study_gxe')




# energytot_imp
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1)


p <- ggplot(figi_gxe, aes(x = outc, y = energytot_imp)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) + 
  facet_wrap(vars(study_gxe), ncol = 5)
p


p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot()
p



# redmeatqc2
figi_gxe <- dplyr::filter(figi_gwas, gxe == 1, EUR_subset == 1) %>% 
  mutate(redmeatqc2 = fct_explicit_na(redmeatqc2, na_level = "NA"),
         redmeatqc2 = fct_relevel(redmeatqc2, "NA", "1", "2", "3", "4"))


p <- ggplot(figi_gxe, aes(outc)) + 
  geom_bar(aes(fill = redmeatqc2), position = 'fill') + 
  theme_bw() + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) + 
  scale_fill_manual(values = c("black", "red", "orange", "lightgreen", "lightblue")) +
  facet_wrap(vars(study_gxe), ncol = 5)
p



geom_bar(aes(fill = asp_ref), position = 'fill') +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c("black", "red", "cyan")) +
  facet_wrap(vars(study_gxe), ncol = 5)
