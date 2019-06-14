#=============================================
# Create PC plots
#
# notes:
# FIGI_Sample_Rename_VCF_07302018.txt N = 145194
# ALL_Merged_PCA_KGP.eigenvec N = 140802
#
# differences: 
# - 2504 KGP samples
# - 6899 Asian samples
# - 3 CCFR samples that are not in the samplefile
#
# 140802 - 2504 + 6899 - 3 = 145194
#
#=============================================
library(dplyr)
library(data.table)
library(ggplot2)

rm(list = ls())
setwd("~/git/FIGI_PCA/")
load("~/git/FIGI_EpiData/FIGI_SampleFile_Summary_08132018.RData")

# download results
#system("scp hpcc:/staging/dvc/andreeki/PCA/ALL_Merged_PCA_KGP.eigenvec ~/git/FIGI_PCA/")
system("cp ~/git/FIGI_TestRun_85k/integrated_call_samples_v3.20130502.ALL.panel.fix ~/git/FIGI_PCA")


header = c("FID", "IID", paste0(rep("PC", 10), seq(1,10)))
pcs <- fread("ALL_Merged_PCA_KGP.eigenvec", skip = 1, col.names = header) 
kgp_samples <- fread("integrated_call_samples_v3.20130502.ALL.panel.fix", stringsAsFactors = F) %>% rename(IID = sample)

All_Merged <- mutate(All_Merged, IID = vcfid, 
										 race_self_new = factor(race_self, labels = c("Unknown", "AI_AN", "Asian", "Black", "NH_PI", "Other", "White")))
table(All_Merged$race_self_new)

# Quick Run

				 

									 
# plot combined dataset
x <- full_join(pcs, kgp_samples, by = 'IID') %>%
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

# FILTER KGP from your result
x.kgp <- filter(x, !is.na(super_pop))

p <- ggplot(data = x.kgp %>% arrange(Group), aes(PC1, PC2, color = Group)) +
	geom_point(alpha = 0.5) +
	labs(x = 'PC1',
			 y = 'PC2') +
	scale_colour_manual(values=c("red", "yellow", "purple", "green", "royalblue")) +
	theme_classic() +
	theme(legend.key.size = unit(0.15, 'inches'))
p


# Whites Only
w <- filter(All_Merged, race_self == "White") %>%
	dplyr::select(vcfid) %>%
	dplyr::rename(IID = vcfid) %>%
	bind_rows(kgp_samples[, 'IID', drop = F])

x.w <- inner_join(x, w, by = 'IID')

p <- ggplot(data = x.w %>% arrange(Group), aes(PC1, PC2, color = Group)) +
	geom_point(alpha = 0.5) +
	labs(x = 'PC1',
			 y = 'PC2',
			 title = "PC1 vs PC2") +
	scale_colour_manual(values=c("black", "red", "yellow", "purple", "green", "royalblue")) +
	theme_classic() +
	theme(legend.key.size = unit(0.15, 'inches'))
p



#------------------------------------------------
# Graphs by subsets
# To be Continued... 
#------------------------------------------------
kgp <- fread("~/work/FIGI_TestRun_85k/integrated_call_samples_v3.20130502.ALL.panel.fix", stringsAsFactors = F) %>%
	rename(IID = sample)
pca <- fread("~/work/FIGI_TestRun_85k/ALL_Merged_20KSNPS_KGP.eigenvec", stringsAsFactors = F, col.names = c("FID", "IID", paste0(rep("PC", 10), seq(1,10))))

# data.plot <- full_join(pca, kgp, by = "IID") %>%
# 	mutate(Group = factor(replace(super_pop,is.na(super_pop), 'FIGI'),
# 												levels=c("FIGI","AFR", "AMR", "EAS", "EUR", "SAS")))


plotpc.race <- function(pc1, pc2, group) {
	# pc1 = "PC1"
	# pc2 = "PC2"
	# group = "White"
	z <- bind_rows(All_Merged[which(All_Merged$race_self_new == group), 'IID', drop=F], kgp[, 'IID', drop=F])
	tmp.dat <- full_join(pca, kgp, by = "IID") %>%
		mutate(Group = factor(replace(super_pop, is.na(super_pop), group),
														levels=c(group,"AFR", "AMR", "EAS", "EUR", "SAS"))) %>% filter(IID %in% z$IID)
	N <- nrow(tmp.dat)-2504
	filename = paste0("~/Dropbox/", group, "_", pc1, "_", pc2, ".png")
		
	p <- ggplot(data = tmp.dat %>% 
							arrange(Group), 
							aes(eval(parse(text = pc1)), eval(parse(text = pc2)), 
									color = Group)) + 
		geom_point(alpha = 0.5) + 
		labs(x = pc1,
				 y = pc2,
				 title = paste0(pc1, " vs. ", pc2, " (N = ", N , ")")) + 
		scale_colour_manual(values=c("black", "red", "yellow", "purple", "green", "royalblue")) +
		theme_classic() +
		theme(legend.key.size = unit(0.15, 'inches'))
	p
	
	ggsave(p, filename = filename, width = 9.5, height = 6)
}




# self_race - pc1, pc2, pc3
table(All_Merged$race_self_new)		

plotpc.race('PC1', 'PC2', group = "Unknown")
plotpc.race('PC1', 'PC2', group = "AI_AN")
plotpc.race('PC1', 'PC2', group = "Asian")
plotpc.race('PC1', 'PC2', group = "Black")
plotpc.race('PC1', 'PC2', group = "NH_PI")
plotpc.race('PC1', 'PC2', group = "Other")
plotpc.race('PC1', 'PC2', group = "White")

plotpc.race('PC1', 'PC3', group = "Unknown")
plotpc.race('PC1', 'PC3', group = "AI_AN")
plotpc.race('PC1', 'PC3', group = "Asian")
plotpc.race('PC1', 'PC3', group = "Black")
plotpc.race('PC1', 'PC3', group = "NH_PI")
plotpc.race('PC1', 'PC3', group = "Other")
plotpc.race('PC1', 'PC3', group = "White")

plotpc.race('PC2', 'PC3', group = "Unknown")
plotpc.race('PC2', 'PC3', group = "AI_AN")
plotpc.race('PC2', 'PC3', group = "Asian")
plotpc.race('PC2', 'PC3', group = "Black")
plotpc.race('PC2', 'PC3', group = "NH_PI")
plotpc.race('PC2', 'PC3', group = "Other")
plotpc.race('PC2', 'PC3', group = "White")


# # KGP ONLY
# # (looks much better, i don't think you can calculate PCs by chunks?)
# # use plink2 --pca approx
# rm(list = ls())
# 
# pca <- fread("~/work/FIGI_TestRun_85k/del.kgp.eigenvec", stringsAsFactors = F, col.names = c("FID", "IID", paste0(rep("PC", 10), seq(1,10))))
# 
# kgp <- fread("~/work/FIGI_TestRun_85k/integrated_call_samples_v3.20130502.ALL.panel.fix", stringsAsFactors = F) %>%
# 	rename(IID = sample)
# 
# x <- inner_join(pca, kgp, by = 'IID') %>%
# 	mutate(Group = factor(super_pop,
# 												levels=c("AFR", "AMR", "EAS", "EUR", "SAS")))
# 
# p <- ggplot(data = x %>% arrange(Group), aes(PC1, PC2, color = Group)) + 
# 	geom_point(alpha = 0.5) + 
# 	labs(x = 'PC1',
# 			 y = 'PC2') + 
# 	scale_colour_manual(values=c("red", "yellow", "purple", "green", "royalblue")) +
# 	theme_classic() +
# 	theme(legend.key.size = unit(0.15, 'inches'))
# p
