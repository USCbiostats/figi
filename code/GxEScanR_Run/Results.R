#========================================
# Process GxEScanR Results
#========================================
library(tidyverse)
library(data.table)
library(qqman)

system("scp hpcc:/staging/dvc/andreeki/figi/results/FIGI_ALL_chr22_NOSTUDY.out ~/")

gxe22 <- fread("~/FIGI_ALL_chr22_NOSTUDY.out") %>% 
  mutate(zP = 2*pnorm(-abs(zG)),
         zlogP = -log10(zP),
         zGxEP = 2*pnorm(-abs(zGxE)), 
         zGxElogP = -log10(zGxEP)) %>% 
  filter(!is.infinite(zlogP))

summary(gxe22$zP)
summary(gxe22$zG)

# N = 47581 cases, 55658 controls
# G Main Effects
qq(gxe22$zP)
manhattan(gxe22, chr="CHR", bp="BP", p="zP")

# GxE 
qq(gxe22$zGxEP)
manhattan(gxe22, chr="CHR", bp="BP", p="zGxEP")

# super significant
gxe22_ss <- filter(gxe22, zlogP > 50)

write.table(gxe22_ss[, c('CHR', 'BP')], file = "~/closerlook.txt", quote = F, row.names = F, col.names = F, sep = '\t')
system("scp ~/closerlook.txt hpcc:/staging/dvc/andreeki/posthoc ")

# I ran the 12 super significant with and without study adjustment
system("scp hpcc:/staging/dvc/andreeki/posthoc/FIGI_Covariates_NOSTUDY_09042018.out ~/")
system("scp hpcc:/staging/dvc/andreeki/posthoc/FIGI_Covariates_09042018.out ~/")


gxe22_ss_nostudy <- fread("FIGI_Covariates_NOSTUDY_09042018.out") %>% 
  mutate(zP = 2*pnorm(-abs(zG)),
         zlogP = -log10(zP),
         zGxEP = 2*pnorm(-abs(zGxE)), 
         zGxElogP = -log10(zGxEP))
gxe22_ss_withstudy <- fread("FIGI_Covariates_09042018.out") %>% 
  mutate(zP = 2*pnorm(-abs(zG)),
         zlogP = -log10(zP),
         zGxEP = 2*pnorm(-abs(zGxE)), 
         zGxElogP = -log10(zGxEP))

manhattan(gxe22_ss_nostudy, chr="CHR", bp="BP", p="zP")
manhattan(gxe22_ss_withstudy, chr="CHR", bp="BP", p="zP") # rs4821082

manhattan(gxe22_ss_nostudy, chr="CHR", bp="BP", p="zGxEP")
manhattan(gxe22_ss_withstudy, chr="CHR", bp="BP", p="zGxEP") # rs4821082








#========================================
# Process GxEScanR Results
#
# without using qqman (too slow)
#========================================
library(tidyverse)
library(data.table)
library(qqman)

rm(list = ls())
source("~/git/FIGI_GxE_FullRun/Results_Functions.R")
tmp <- do.call(rbind, lapply(list.files("~/worktmp", full.names = T), fread, stringsAsFactors = F))

tmp2 <- tmp %>% 
  mutate(zP = 2*pnorm(-abs(zG)),
         zlogP = -log10(zP),
         zGxEP = 2*pnorm(-abs(zGxE)), 
         zGxElogP = -log10(zGxEP),
         CHR = as.integer(CHR),
         BP = as.numeric(BP)) %>% 
  arrange(CHR, BP) %>% 
  filter(!is.infinite(zlogP)) # no more? wtf?



# Marginal G
png('FIGI_ALL_chr1_22_G_manhattan_incomplete.png')
manhattan.plot(tmp2$CHR, tmp2$BP, tmp2$zP, sig.level = 1e-8)
dev.off()

png('FIGI_ALL_chr1_22_G_qq_incomplete.png')
ggd.qqplot(tmp2$zP)
dev.off()

tmp_sig <- filter(tmp2, zlogP > 5) # how many SNPs are highly significant, use cutoff of logp > 8
nrow(tmp_sig)

source("~/git/FIGI_GxE_FullRun/Results_Functions.R")
png('FIGI_ALL_chr1_22_G_qq_subset_incomplete.png')
ggd.qqplot.exclude(tmp2$zP)
dev.off()




# GxE
png('FIGI_ALL_chr1_22_GxE_manhattan_incomplete.png')
manhattan.plot(tmp2$CHR, tmp2$BP, tmp2$zGxEP, sig.level = 8)
dev.off()

png('FIGI_ALL_chr1_22_GxE_qq_incomplete.png')
ggd.qqplot(tmp2$zGxEP)
dev.off()

tmp_sig <- filter(tmp2, zGxElogP > 5) # how many SNPs are highly significant, use cutoff of logp > 5 (suggestive)
nrow(tmp_sig)

source("~/git/FIGI_GxE_FullRun/Results_Functions.R")
png('FIGI_ALL_chr1_22_GxE_qq_subset_incomplete.png')
ggd.qqplot.exclude.gxe(tmp2$zGxEP)
dev.off()





