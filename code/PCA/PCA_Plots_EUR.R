#=============================================
# PCA Pairwise Plots
# 12/19/2019
#
# with and without subsetting by EUR only 
#
#=============================================
library(tidyverse)
library(data.table)
rm(list = ls())

load("~/data/FIGI_EpiData_rdata/FIGI_HRC_v2.3_GWAS.rds")

header = c("FID", "IID", paste0(rep("PC", 20), seq(1,20)))
pcs <- fread("~/data/PCA/190729/FIGI_GwasSet_KGP_190729.eigenvec", skip = 1, col.names = header)
kgp <- fread("integrated_call_samples_v3.20130502.ALL.panel.fix", stringsAsFactors = F) %>% rename(IID = sample)

figi_pca <- figi %>% 
  mutate(IID = vcfid, 
         race_self_new = factor(race_self, labels = c("Unknown", "AI_AN", "Asian", "Black", "NH_PI", "Other", "White")))

table(figi_pca$race_self_new, figi_pca$EUR_subset, useNA = 'ifany')



x <- full_join(pcs, kgp, by = 'IID') %>%
  mutate(Group = factor(replace(super_pop,is.na(super_pop), 'FIGI'),
                        levels=c("FIGI","AFR", "AMR", "EAS", "EUR", "SAS")))

p <- ggplot(data = x %>% arrange(Group), aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.5) +
  labs(x = 'PC1',
       y = 'PC2',
       title = "PC1 vs PC2") +
  scale_colour_manual(values=c("black", "red", "yellow", "purple", "green", "royalblue")) +
  theme_classic() +
  theme(legend.key.size = unit(0.15, 'inches'))
p



# GWAS set plot, with and without EUR samples

keep_ids <- c(pcs$IID)

x <- full_join(pcs, figi_pca, 'IID')

y <- full_join(x, kgp, by = 'IID') %>%
  mutate(Group = factor(replace(super_pop,is.na(super_pop), 'FIGI'),
                        levels=c("FIGI","AFR", "AMR", "EAS", "EUR", "SAS"))) %>% 
  filter(IID %in% keep_ids)


p <- ggplot(data = y %>% arrange(Group), aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.5) +
  labs(x = 'PC1',
       y = 'PC2',
       title = "PC1 vs PC2") +
  scale_colour_manual(values=c("black", "red", "yellow", "purple", "green", "royalblue")) +
  theme_classic() +
  theme(legend.key.size = unit(0.15, 'inches'))
p
ggsave(p, filename = "~/Dropbox/FIGI/Results/PCA/plots/figi_gwas_pc1_pc2.png", width = 9.5, height = 6)




# EUR ONLY (annoying)
kgp_subset <- filter(pcs, IID %in% kgp$IID)

x <- inner_join(pcs, figi_pca, "IID") %>% 
  filter(EUR_subset == 1)

y <- bind_rows(x, kgp_subset) %>% 
  full_join(kgp, 'IID') %>% 
  mutate(Group = factor(replace(super_pop,is.na(super_pop), 'FIGI'),
                        levels=c("FIGI","AFR", "AMR", "EAS", "EUR", "SAS")))


p <- ggplot(data = y %>% arrange(Group), aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.5) +
  labs(x = 'PC1',
       y = 'PC2',
       title = "PC1 vs PC2") +
  scale_colour_manual(values=c("black", "red", "yellow", "purple", "green", "royalblue")) +
  theme_classic() +
  theme(legend.key.size = unit(0.15, 'inches'))
p

ggsave(p, filename = "~/Dropbox/FIGI/Results/PCA/plots/figi_gwas_pc1_pc2_EUR.png", width = 9.5, height = 6)




