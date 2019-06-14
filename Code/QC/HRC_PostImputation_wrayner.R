#=============================================================================#
#
# wrayner IC (postimputation QC tool)
# post-imputation QC summary + plots
# 
# create new plots, separate by platform N samples and platform Type
#
#
# useful - gather files of a certain extension in a folder:
# find ./IC -name \*GW.txt -exec cp {} ./IC_working/ \;
#=============================================================================#
library(tidyverse)
library(data.table)
setwd("~/data/IC_working/")


# ------ INFO Score ------
# running on single batch
# x <- fread("INFO-axiom_acs_aus_nf.GW.txt") %>% 
#   filter(V1 != "Total") %>% 
#   mutate(tmp = ifelse(V1 <= 0.3, 0, 
#                       ifelse(V1 > 0.3 & V1 < 0.8, 1, 2)),
#          Info_cat = factor(tmp, labels = c("<=0.3", ">0.3 & <0.8", ">=0.8"))) %>% 
#   group_by(Info_cat) %>% 
#   summarise(SNPs = sum(V2), Pct = sum(V3)) %>% 
#   mutate(batch = "axiom_acs_aus_nf")
# 
# y <- fread("INFO-axiom_mecc_cfr_ky.GW.txt") %>% 
#   filter(V1 != "Total") %>% 
#   mutate(tmp = ifelse(V1 <= 0.3, 0, 
#                       ifelse(V1 > 0.3 & V1 < 0.8, 1, 2)),
#          Info_cat = factor(tmp, labels = c("<=0.3", ">0.3 & <0.8", ">=0.8"))) %>% 
#   group_by(Info_cat) %>% 
#   summarise(SNPs = sum(V2), Pct = sum(V3)) %>% 
#   mutate(batch = "axiom_mecc_cfr_ky")
# 
# df <- bind_rows(x, y)


# generalize to all imputation batches
quickw <- function(x) {
  fread(paste0("INFO-", x, ".GW.txt")) %>% 
    filter(V1 != "Total") %>% 
    mutate(tmp = ifelse(V1 <= 0.3, 0, 
                        ifelse(V1 > 0.3 & V1 < 0.8, 1, 2)),
           Info_cat = factor(tmp, labels = c("<=0.3", ">0.3 & <0.8", ">=0.8"))) %>% 
    group_by(Info_cat) %>% 
    summarise(SNPs = sum(V2), Pct = sum(V3)) %>% 
    mutate(batch = x)
}


batch_list <- c("axiom_acs_aus_nf", 
                "axiom_mecc_cfr_ky",
                "ccfr_1m_1mduo_reimpute",
                "ccfr_omni",
                "corect_oncoarray", 
                "corsa_axiom",
                "cytosnp_comb",
                "initial_comb_datasets",
                "mecc",
                "newfoundland_omniquad",
                "omni_comb",
                "omniexpress_exomechip",
                "oncoarray_to_usc",
                "plco_3",
                "reach", 
                "ukbiobank")

df <- do.call(rbind, lapply(batch_list, function(x) quickw(x)))

sample_sizes <- c('2766', '7501', '2180', '1600', '36621', '2467', '10908', '5924', '983', '637', '5986', '5503', '20912', '4864', '750', '27594')
df$samplesize <- rep(sample_sizes, each = 3)
df$batch_ss <- paste0(df$batch, " (", df$samplesize, ")")

# add platform information
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")
table(figi$filename, figi$platform)

platform <- c("Axiom", "Axiom", "Illumina_1M_1Mduo", "Illumina_Omni", "OncoArray", "Axiom", "CytoSNP", "InitialGWAS", "Omni2.5", "Omniquad", "OmniExpress", "OmniExpress_ExomeChip", "OncoArray_Custom", "OncoArray", "OncoArray_Custom", "Axiom")
df$platform <- rep(platform, each = 3)

ggplot(df) + 
  geom_bar(aes(x = batch_ss, y = Pct, fill = Info_cat), stat='identity') +
  # geom_text(aes(x = batch_ss, y = rev(Pct), label = round(Pct,1))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


ggplot(df, aes(batch_ss, Pct, label = Pct)) + 
  geom_bar(aes(fill = Info_cat), stat='identity', position = 'stack') +
  # geom_text(aes(label = Pct), position = position_stack(vjust = 0.5, reverse = T)) +
  coord_flip() + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12)) + 
  geom_hline(yintercept = 50, linetype = 'dashed') + 
  geom_hline(yintercept = 25, linetype = 'dashed') 


# separate these graphs by sample size and by platform to explorations
# facet_grid

test <- df %>% 
  arrange(platform, batch_ss) %>% 
  filter(platform == "Axiom")

test <- df %>% 
  arrange(platform, batch_ss)


ggplot(test, aes(Pct, batch_ss)) + 
  geom_barh(aes(fill = Info_cat), stat = "identity") + 
  facet_grid(platform  ~ . , scales = "free", space = "free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 0),
        axis.text = element_text(size = 12)) + 
  geom_vline(xintercept = 50, linetype = "dashed") 





# ------ Allele Frequency ------
x <- fread("AF-axiom_acs_aus_nf.GW.txt") %>%
  filter(V1 != "Total") %>% 
  mutate(group = "wtf")
  
df$group <- rep(batch_list, each = 11)





## generalize to all imputation batches
quickw <- function(x) {
  fread(paste0("AF-", x, ".GW.txt")) %>% 
    filter(V1 != "Total") %>% 
    mutate(group = x,
           afbin = as.numeric(V1),
           maf = as.character(ifelse(V1 > 0.5, 1-afbin, afbin))) %>% 
    group_by(group, maf) %>%
    summarise(SNPs = sum(V2), Pct = sum(V3))
}

x<-"axiom_acs_aus_nf"
wtf <- quickw("axiom_acs_aus_nf")

batch_list <- c("axiom_acs_aus_nf", 
                "axiom_mecc_cfr_ky",
                "ccfr_1m_1mduo_reimpute",
                "ccfr_omni",
                "corect_oncoarray", 
                "corsa_axiom",
                "cytosnp_comb",
                "initial_comb_datasets",
                "mecc",
                "newfoundland_omniquad",
                "omni_comb",
                "omniexpress_exomechip",
                "oncoarray_to_usc",
                "plco_3",
                "reach", 
                "ukbiobank")

# density
df <- do.call(rbind, lapply(batch_list, function(x) quickw(x)))
ggplot(df, aes(V2)) + 
  geom_area()


df <- do.call(rbind, lapply(batch_list, function(x) quickw(x)))
ggplot(df, aes(maf, Pct, group = group)) + 
  geom_line(aes(color = group))


df <- do.call(rbind, lapply(batch_list, function(x) quickw(x))) %>% 
  filter(maf != "0")
ggplot(df, aes(maf, Pct, group = group)) + 
  geom_line(aes(color = group))




###############################################################################
###############################################################################
###############################################################################
###############################################################################

# #=============================================================================#
# # 
# # wrayner IC results
# # post-imputation QC summary + plots
# # create own summary here
# #
# #=============================================================================#
# library(tidyverse)
# library(data.table)
# setwd("~/IC")
# 
# # INFO score
# # let's split into three categories, < 0.3, < 0.8, > 0.8
# x <- fread("axiom_acs_aus_nf/axiom_acs_aus_nf.27-4-2019/INFO-axiom_acs_aus_nf.GW.txt") %>% 
#   filter(V1 != "Total") %>% 
#   mutate(tmp = ifelse(V1 <= 0.3, 0, 
#                       ifelse(V1 > 0.3 & V1 < 0.8, 1, 2)),
#          Info_cat = factor(tmp, labels = c("<=0.3", ">0.3 & <0.8", ">=0.8"))) %>% 
#   group_by(Info_cat) %>% 
#   summarise(SNPs = sum(V2), Pct = sum(V3)) %>% 
#   mutate(batch = "axiom_acs_aus_nf")
# 
# y <- fread("axiom_mecc_cfr_ky/axiom_mecc_cfr_ky.28-4-2019/INFO-axiom_mecc_cfr_ky.GW.txt") %>% 
#   filter(V1 != "Total") %>% 
#   mutate(tmp = ifelse(V1 <= 0.3, 0, 
#                       ifelse(V1 > 0.3 & V1 < 0.8, 1, 2)),
#          Info_cat = factor(tmp, labels = c("<=0.3", ">0.3 & <0.8", ">=0.8"))) %>% 
#   group_by(Info_cat) %>% 
#   summarise(SNPs = sum(V2), Pct = sum(V3)) %>% 
#   mutate(batch = "axiom_mecc_cfr_ky")
# 
# df <- bind_rows(x, y)
# 
# ggplot(df) + 
#   geom_bar(aes(x = batch, y = Pct, fill =), stat='identity')
# 
