#=============================================================================#
# comparing Rsq and MAF between different imputation/phasing software
#
# (if they're too different, that's bad)
#
#=============================================================================#
library(tidyverse)
library(data.table)

# info files
minimac3_shapeit <- fread("/home/rak/data/HRC_InfoFile/ccfr_1m_1mduo_reimpute/chr22.info.gz") %>% 
  mutate(id = paste(SNP, `REF(0)`, `ALT(1)`, sep = ":"))

minimac4 <- fread("~/data/ccfr_1m_reimpute_minimac3_vs_minimac4/ccfr_1m_reimpute_minimac4.info.gz") %>% 
  mutate(id = SNP)

minimac3 <- fread("~/data/ccfr_1m_reimpute_minimac3_vs_minimac4/chr22.info.gz") %>% 
  mutate(id = paste(SNP, `REF(0)`, `ALT(1)`, sep = ":"))



#-----------------------------------------------------------------#
# shapeit_minimac3 vs eagle_minimac4
z <- inner_join(minimac3_shapeit[, c('id', 'Rsq', 'MAF')], minimac4[, c('id', 'Rsq', 'MAF')], by = 'id')

plot(z$Rsq.x, z$Rsq.y, main = "shapeit_minimac3 vs eagle_minimac4 - Rsq (all markers)", xlab = 'shapeit_minimac3 Rsq', ylab = 'eagle_minimac4 Rsq')
abline(0, 1, col = 'red')

plot(z$MAF.x, z$MAF.y, main = "shapeit_minimac3 vs eagle_minimac4 - MAF (all markers)", xlab = 'shapeit_minimac3 MAF', ylab = 'eagle_minimac4 MAF')
abline(0, 1, col = 'red')



# shapeit_minimac3 vs eagle_minimac4 (MAF > 1%)
minimac3_shapeit_filter <- filter(minimac3_shapeit, MAF > 0.01)
minimac4_filter         <- filter(minimac4, MAF > 0.01)

z <- inner_join(minimac3_shapeit_filter[, c('id', 'Rsq', 'MAF')], minimac4_filter[, c('id', 'Rsq', 'MAF')], by = 'id')

plot(z$Rsq.x, z$Rsq.y, main = "shapeit_minimac3 vs eagle_minimac4 - Rsq (MAF > 1%)", xlab = 'shapeit_minimac3 Rsq', ylab = 'eagle_minimac4 Rsq')
abline(0, 1, col = 'red')

plot(z$MAF.x, z$MAF.y, main = "shapeit_minimac3 vs eagle_minimac4 - MAF (MAF > 1%)", xlab = 'shapeit_minimac3 MAF', ylab = 'eagle_minimac4 MAF')
abline(0, 1, col = 'red')




#-----------------------------------------------------------------#
# eagle_minimac3 vs eagle_minimac4
z <- inner_join(minimac3[, c('id', 'Rsq', 'MAF')], minimac4[, c('id', 'Rsq', 'MAF')], by = 'id')

plot(z$Rsq.x, z$Rsq.y, main = "eagle_minimac3 vs eagle_minimac4 - Rsq (all markers)", xlab = 'eagle_minimac3 Rsq', ylab = 'eagle_minimac4 Rsq') # bad.. 
abline(0, 1, col = 'red')

plot(z$MAF.x, z$MAF.y, main = "eagle_minimac3 vs eagle_minimac4 - MAF (all markers)", xlab = 'eagle_minimac3 MAF', ylab = 'eagle_minimac4 MAF')# bad.. 
abline(0, 1, col = 'red')



# eagle_minimac3 vs eagle_minimac4 MAF > 1%
minimac4_filter <- filter(minimac4, MAF > 0.01)
minimac3_filter <- filter(minimac3, MAF > 0.01)

z <- inner_join(minimac3_filter, minimac4_filter, by = 'id')

plot(z$Rsq.x, z$Rsq.y, main = "eagle_minimac3 vs eagle_minimac4 - Rsq (MAF > 1%)", xlab = 'eagle_minimac3 Rsq', ylab = 'eagle_minimac4 Rsq') # bad.. 
abline(0, 1, col = 'red')

plot(z$MAF.x, z$MAF.y, main = "eagle_minimac3 vs eagle_minimac4 - MAF (MAF > 1%)", xlab = 'eagle_minimac3 MAF', ylab = 'eagle_minimac4 MAF') # bad.. 
abline(0, 1, col = 'red')
