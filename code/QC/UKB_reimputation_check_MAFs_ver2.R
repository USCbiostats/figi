#========================================================================#
# Compare Allele Frequencies between corect_oncoarray and ukbiobank
# To find discrepancies (flipped/mismatched AAF on ukbiobank..)
#
# maybe it's more legible in manhattan plot format.. ?
#========================================================================#
library(tidyverse)
library(data.table)
library(R.utils)


# Test Chr22 (all samples)
x_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_new/ukbiobank_chr22.rds")
x_aaf <- readRDS("~/ukbiobank_check/ukb_new_all/ukb_aaf_all_chr22.rds")

x <- cbind(x_info[[11]], x_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         maf = 0.5 - abs(x_aaf - 0.5))


y_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_old/ukbiobank_chr22.rds")
y_aaf <- readRDS("~/ukbiobank_check/ukb_old_all/ukb_aaf_all_chr22.rds")

y <- cbind(y_info[[11]], y_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"),
         maf = 0.5 - abs(y_aaf - 0.5))

z <- inner_join(x, y, by = 'id') %>% 
  mutate(a1 = -log10(x_aaf), 
         a2 = -log10(y_aaf),
         pct_diff = (abs(x_aaf-y_aaf) / ((x_aaf+y_aaf) / 2)),
         pct_diff_log = (abs(a1-a2) / ((a1 + a2) / 2)))

wtf <- filter(z, Location.x == "40651008")
head(z)

hist(z$pct_diff_log)
length(which(z$pct_diff_log > 0.5))

z_filter <- filter(z, pct_diff_log > 0.5)

ggplot(data=z_filter, aes(a1, a2)) +
  geom_point() + geom_abline()


ggplot(data=z_filter, aes(x_aaf, y_aaf)) +
  geom_point() + geom_abline()









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










# Test
x_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_new/ukbiobank_chr3.rds")
x_aaf <- readRDS("~/ukbiobank_check/aaf_ukb_controls_new/ukbiobank_aaf_chr3.rds")

x <- cbind(x_info[[11]], x_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"), 
         maf = 0.5 - abs(x_aaf - 0.5))

y_info <- readRDS("~/ukbiobank_check/bdoseinfo_ukb_old/ukbiobank_chr3.rds")
y_aaf <- readRDS("~/ukbiobank_check/aaf_ukb_controls_old/ukbiobank_aaf_chr3.rds")

y <- cbind(y_info[[11]], y_aaf) %>% 
  mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"),
         maf = 0.5 - abs(y_aaf - 0.5))

z <- inner_join(x, y, by = 'id') %>% 
  mutate(a1 = -log10(maf.x), 
         a2 = -log10(maf.y),
         pct_diff_log = (abs(a1-a2) / (a1 + a2) / 2))

hist(z$pct_diff_log)
length(which(z$pct_diff_log > 0.01))

z_filter <- filter(z, pct_diff_log > 0.01)


ggplot(data=z_filter, aes(a1, a2)) +
  geom_point()
  
  
z <- inner_join(x, y, by = 'id') %>% 
  mutate(pct_diff = (abs(corect_oncoarray_MAF - ukbiobank_MAF) / ( (corect_oncoarray_MAF + ukbiobank_MAF) / 2)),
         pct_diff2 = pct_diff * abs(corect_oncoarray_MAF - ukbiobank_MAF),
         a1 = -log10(corect_oncoarray_MAF), 
         a2 = -log10(ukbiobank_MAF), 
         a_diff = ( abs(a1 - a2) / ( (a1 + a2) / 2)), 
         ln_diff = log(corect_oncoarray_MAF/ukbiobank_MAF))

zz <- filter(z, a_diff > 0.20)

test <- gather(zz, key = "batch", value = "maf", -id, -pct_diff) %>% 
  separate(id, into = c("chr", "pos", "ref", "alt")) 

test <- zz %>%  
  separate(id, into = c("chr", "pos", "ref", "alt")) 


ggplot(data = test, aes(corect_oncoarray_MAF, ukbiobank_MAF)) +
  geom_point()


# loop through all chromosomes. (takes a while)
chr <- 22
for(chr in 1:22) {
  x <- fread(paste0("~/data/HRC_InfoFile/corect_oncoarray/chr", chr, ".info.gz")) %>%
    mutate(id = paste(SNP, `REF(0)`, `ALT(1)`, sep = ":"),  
           corect_oncoarray_AAF = ALT_Frq) %>% 
    dplyr::select(id, corect_oncoarray_AAF)
  
  y <- fread(paste0("~/data/HRC_InfoFile/ukbiobank/chr", chr, ".info.gz")) %>%
    rename(id = SNP, 
           ukbiobank_AAF = ALT_Frq) %>% 
    dplyr::select(id, ukbiobank_AAF)
  
  z <- inner_join(x, y, by = 'id')
  
  saveRDS(z, file = paste0("~/corect_oncoarray_vs_ukb_chr", chr, ".rds"))
}



# closer look.. 
quick_wrap <- function(chr) {
  z <- readRDS(paste0("~/corect_oncoarray_vs_ukb_chr", chr, ".rds")) %>% 
    mutate(pct_diff = ukbiobank_AAF/corect_oncoarray_AAF) %>% 
    filter(!(corect_oncoarray_AAF < 0.05 & ukbiobank_AAF < 0.05),
           !(corect_oncoarray_AAF > 0.95 & ukbiobank_AAF > 0.95),
           (pct_diff > 2 | pct_diff < 0.5))
  plot(z$corect_oncoarray_AAF, z$ukbiobank_AAF,
       main = paste0("Comparison of Alt Allele Freq - Chr", chr),
       xlab = "corect_oncoarray AAF", ylab = "ukbiobank AAF")
}



# use wrapper to generate plots

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr1-4.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(1)
quick_wrap(2)
quick_wrap(3)
quick_wrap(4)
dev.off()

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr5-8.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(5)
quick_wrap(6)
quick_wrap(7)
quick_wrap(8)
dev.off()

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr9-12.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(9)
quick_wrap(10)
quick_wrap(11)
quick_wrap(12)
dev.off()

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr13-16.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(13)
quick_wrap(14)
quick_wrap(15)
quick_wrap(16)
dev.off()

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr17-20.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(17)
quick_wrap(18)
quick_wrap(19)
quick_wrap(20)
dev.off()

png("~/Dropbox/AAF_Comparison_corect_oncoarray_vs_ukb_Chr20-22.png", height = 800, width = 800)
par(mfrow = c(2,2))
quick_wrap(21)
quick_wrap(22)
dev.off()



# closer look at individual markers
# start with chromosome 16 (NOT a flip, just inconsistent AAFs.. )
setwd("~/yuru")
z <- readRDS(paste0("~/yuru/corect_oncoarray_ukbiobank_AAF_comparison_chr", 16, ".rds")) %>% 
  filter(grepl("79591072", id))

# Get rsnumber for chr16:79591072
system("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr16\" and chromEnd = 79591072' > chr16_79591072.sql")
system('cat chr16_79591072.sql')
system('mysql -u figi -pgecco -A -D hg19 --skip-column-names < chr16_79591072.sql', ignore.stderr = T)

# Get rsnumber for chr18:3780526
system("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr18\" and chromEnd = 3780526' > chr18_3780526.sql")
system('cat chr18_3780526.sql')
system('mysql -u figi -pgecco -A -D hg19 --skip-column-names < chr18_3780526.sql', ignore.stderr = T)

z <- rbind(paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr2\" and chromEnd = 4964187'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr3\" and chromEnd = 98981063'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr6\" and chromEnd = 138017298'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr7\" and chromEnd = 99995536'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr10\" and chromEnd = 3001065'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr16\" and chromEnd = 79591072'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr17\" and chromEnd = 18528708'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr18\" and chromEnd = 3780526'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr21\" and chromEnd = 34169318'\n"),
           paste0("echo 'select chrom,chromEnd,name from snp138 where chrom =  \"chr22\" and chromEnd = 40651008'\n"))

z <- rbind(paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr2\" and chromEnd = 4964187;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr3\" and chromEnd = 98981063;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr6\" and chromEnd = 138017298;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr7\" and chromEnd = 99995536;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr10\" and chromEnd = 3001065;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr16\" and chromEnd = 79591072;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr17\" and chromEnd = 18528708;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr18\" and chromEnd = 3780526;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr21\" and chromEnd = 34169318;"),
           paste0("select chrom,chromEnd,name from snp138 where chrom =  \"chr22\" and chromEnd = 40651008;"))


write.table(z, file = "~/yuru/bad_snps_10.sql", quote = F, row.names = F, col.names = F)
system('mysql -u figi -pgecco -A -D hg19 --skip-column-names < bad_snps_10.sql', ignore.stderr = T) # copy paste...

grep rs13033564\|rs115593859\|rs6911600\|rs201649203\|rs2454818\|rs420757\|rs60175205\|rs2936649\|rs6506132\|rs34957938\|rs28662924\|rs33987705\|rs115266833


# let's look at these markers... 
quick_wrap <- function(chr) {
  z <- readRDS(paste0("~/yuru/corect_oncoarray_ukbiobank_AAF_comparison_chr", chr, ".rds")) %>% 
    mutate(pct_diff = ukbiobank_AAF/corect_oncoarray_AAF) %>% 
    filter(!(corect_oncoarray_AAF < 0.05 & ukbiobank_AAF < 0.05),
           !(corect_oncoarray_AAF > 0.95 & ukbiobank_AAF > 0.95),
           (pct_diff > 2 | pct_diff < 0.5))
  plot(z$corect_oncoarray_AAF, z$ukbiobank_AAF,
       main = paste0("Comparison of Alt Allele Freq - Chr", chr),
       xlab = "corect_oncoarray AAF", ylab = "ukbiobank AAF")
}