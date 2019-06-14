#========================================================================#
# Compare Allele Frequencies between ukbiobank and corect_oncoarray
# 
# Want to answer the question of minimac3 vs minimac4
# do AAFs differ only for poorly imputed markers? or maybe depends more on
# haplotypes? 
#
#========================================================================#
library(tidyverse)
library(data.table)
library(R.utils)


# we can simply use HRC info files
# I previously created merged info files (all imputation batches)
# use that file to explore MAFs between UKB and corect_oncoarray
# batch_list <- c("axiom_acs_aus_nf",
#                 "axiom_mecc_cfr_ky",
#                 "ccfr_1m_1mduo_reimpute",
#                 "ccfr_omni",
#                 "corect_oncoarray",
#                 "corsa_axiom",
#                 "cytosnp_comb",
#                 "initial_comb_datasets",
#                 "mecc",
#                 "newfoundland_omniquad",
#                 "omni_comb",
#                 "omniexpress_exomechip",
#                 "oncoarray_to_usc",
#                 "plco_3",
#                 "reach",
#                 "ukbiobank")
# sample_sizes <- c(2766, 7501, 2180, 1600, 36621, 2467, 10908, 5924, 983, 637, 5986, 5503, 20912, 4864, 750, 27594)



# start with chromosome 22
chr22 <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr22.rds") %>% 
  dplyr::select(ID, contains("corect_oncoarray"), contains("ukbiobank"))

chr22_maf <- chr22 %>% 
  mutate(corect_oncoarray_MAF = 0.5 - abs(corect_oncoarray_ALT_Frq - 0.5),
         ukbiobank_MAF = 0.5 - abs(ukbiobank_ALT_Frq - 0.5),
         avg_rsq = ((36621*as.numeric(corect_oncoarray_Rsq)) + (27594*as.numeric(ukbiobank_Rsq))) / 64215,
         avg_maf = ((36621*as.numeric(corect_oncoarray_MAF)) + (27594*as.numeric(ukbiobank_MAF))) / 64215,
         a1 = -log10(corect_oncoarray_MAF),
         a2 = -log10(ukbiobank_MAF),
         pct_diff = (abs(corect_oncoarray_MAF-ukbiobank_MAF) / ((corect_oncoarray_MAF+ukbiobank_MAF) / 2)),
         pct_diff_log = (abs(a1-a2) / ((a1 + a2) / 2))) %>% 
  filter(avg_maf > 0.01,
         pct_diff_log > 0.25)

length(which(chr22_maf$pct_diff > 0.5))

ggplot(data=chr22_maf, aes(corect_oncoarray_MAF, ukbiobank_MAF)) +
  geom_point() + geom_abline()


# what if you filtered by Rsq > 0.8?
chr22_maf <- chr22 %>% 
  mutate(corect_oncoarray_MAF = 0.5 - abs(corect_oncoarray_ALT_Frq - 0.5),
         ukbiobank_MAF = 0.5 - abs(ukbiobank_ALT_Frq - 0.5),
         avg_rsq = ((36621*as.numeric(corect_oncoarray_Rsq)) + (27594*as.numeric(ukbiobank_Rsq))) / 64215,
         avg_maf = ((36621*as.numeric(corect_oncoarray_MAF)) + (27594*as.numeric(ukbiobank_MAF))) / 64215,
         a1 = -log10(corect_oncoarray_MAF),
         a2 = -log10(ukbiobank_MAF),
         pct_diff = (abs(corect_oncoarray_MAF-ukbiobank_MAF) / ((corect_oncoarray_MAF+ukbiobank_MAF) / 2)),
         pct_diff_log = (abs(a1-a2) / ((a1 + a2) / 2))) %>% 
  filter(avg_maf > 0.01,
         avg_rsq > 0.8,
         pct_diff_log > 0.5)

length(which(chr22_maf$pct_diff > 0.5))

ggplot(data=chr22_maf, aes(corect_oncoarray_MAF, ukbiobank_MAF)) +
  geom_point() + geom_abline()





# let's expand to all chromosomes (takes a while to run)
chroms <- seq(1,22)
# chroms <- 22
# get info from each chromosome info file, create matrix
get_rsq <- function() {
  do.call(rbind, lapply(chroms, function(chr)
    readRDS(paste0("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr", chr, ".rds")) %>%
      dplyr::select(ID, contains("corect_oncoarray"), contains("ukbiobank")) %>% 
      mutate(corect_oncoarray_MAF = 0.5 - abs(corect_oncoarray_ALT_Frq - 0.5),
             ukbiobank_MAF = 0.5 - abs(ukbiobank_ALT_Frq - 0.5),
             avg_rsq = ((36621*as.numeric(corect_oncoarray_Rsq)) + (27594*as.numeric(ukbiobank_Rsq))) / 64215,
             avg_maf = ((36621*as.numeric(corect_oncoarray_MAF)) + (27594*as.numeric(ukbiobank_MAF))) / 64215,
             a1 = -log10(corect_oncoarray_MAF),
             a2 = -log10(ukbiobank_MAF),
             pct_diff = (abs(corect_oncoarray_MAF-ukbiobank_MAF) / ((corect_oncoarray_MAF+ukbiobank_MAF) / 2)),
             pct_diff_log = (abs(a1-a2) / ((a1 + a2) / 2))) %>% 
      filter(avg_maf > 0.01)))
      # filter(avg_maf > 0.01,
      #        pct_diff_log > 0.25)))
}

x <- get_rsq()
# saveRDS(x, file = "files/corect_oncoarray_ukb_MAF_comparison.rds", version = 2)
names(x)

x2 <- filter(x, pct_diff_log > 0.25)

ggplot(data=x2, aes(corect_oncoarray_MAF, ukbiobank_MAF)) +
  geom_point(size = 0.5) + geom_abline() + ggtitle("UKB vs corect, PctDiff logMAF > 0.25", subtitle = "32,308 SNPs") + xlab("corect_oncoarray MAF") + ylab("UKB MAF")
ggsave("figures/ukb_vs_corectoncoarray_pctdiff.png", width=6, height=4)


x2_filter <- x2 %>% 
  filter(avg_rsq > 0.8) %>% 
  separate(ID, into = c('chr', 'pos', 'ref', 'alt'), sep = ":")

ggplot(data=x2_filter, aes(corect_oncoarray_MAF, ukbiobank_MAF)) +
  geom_point(size = 0.5) + geom_abline() + ggtitle("UKB vs corect, PctDiff logMAF > 0.25, AvgRsq > 0.8", subtitle = "22,750 SNPs") + xlab("corect_oncoarray MAF") + ylab("UKB MAF")
ggsave("figures/ukb_vs_corectoncoarray_pctdiff_rsqfilter.png", width=6, height=4)



table(x_filter$chr)



# are your 126 top hits in there..
# (121 out of 126, probably smaller difference but still driving results, maybe doesn't survive weighted Rsq that includes all imputation batches)

z <- inner_join(results_GxE, x, by = 'ID')
















# Test Chr22 (control samples)
x_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_new/ukbiobank_chr22.rds")
x_aaf <- readRDS("~/ukbiobank_check/ukb_new_controls/ukbiobank_aaf_chr22.rds")

x <- cbind(x_info[[11]], x_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         maf = 0.5 - abs(x_aaf - 0.5))


y_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_old/ukbiobank_chr22.rds")
y_aaf <- readRDS("~/ukbiobank_check/ukb_old_controls/ukbiobank_aaf_chr22.rds")

y <- cbind(y_info[[11]], y_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"),
         maf = 0.5 - abs(y_aaf - 0.5))

z <- inner_join(x, y, by = 'id') %>% 
  mutate(a1 = -log10(x_aaf), 
         a2 = -log10(y_aaf),
         pct_diff = (abs(x_aaf-y_aaf) / ((x_aaf+y_aaf) / 2)),
         pct_diff_log = (abs(a1-a2) / ((a1 + a2) / 2)))

wtf <- filter(z, Location.x == "40651008")

hist(z$pct_diff_log)
length(which(z$pct_diff_log > 0.5))

z_filter <- filter(z, pct_diff_log > 0.5)


ggplot(data=z_filter, aes(a1, a2)) +
  geom_point()


ggplot(data=z_filter, aes(x_aaf, y_aaf)) +
  geom_point()








