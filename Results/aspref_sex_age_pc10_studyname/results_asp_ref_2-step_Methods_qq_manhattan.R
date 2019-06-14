#========================================================#
# Manhattan and QQ Plots
# 2-step methods
# (qqman)
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
rm(list = ls())

# keep in mind, some SNP duplicates in results, identical values, remove for cleanliness
# reason for duplicates is slight error in reading in chunks of SNPs in GxEScanR

# G, GxE, Joint
results_main <- do.call(rbind, lapply(list.files(path = "./results/", full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438"), fread, stringsAsFactors = F)) %>% 
  filter(!duplicated(SNP))
names(results_main); any(duplicated(results_main$SNP))

# GE, Case Only, Control Only
colNames = c("SNP", "Chromosome", "Location", "Reference", "Alternate", "Subjects", "Cases", "betaGE", "chiSqGE", "betaCase", "ChiSqCase",  "betaControl", "chiSqControl", "drop")
results_ge <- do.call(rbind, lapply(list.files(path = "./results_geonly/", full.names = T, pattern = "results_asp_ref_age_sex_pc10_studyname_GEOnly"), fread, stringsAsFactors = F, header = F, col.names = colNames)) %>% 
  filter(!duplicated(SNP))
names(results_ge); any(duplicated(results_ge$SNP))


# calculate lambdas ------
getlambda <- function(pvals) {
  chisq <- qchisq(1-pvals, 1)
  lambda <- round(median(chisq)/qchisq(0.5,1),4)
  lambda
}
getlambda2df <- function(pvals) {
  chisq <- qchisq(1-pvals, 2)
  lambda <- round(median(chisq)/qchisq(0.5,2),4)
  lambda
}
getlambda1000 <- function(lambda) {
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  lambda1000
}

results_main_pvals <- results_main %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         chiSqG.pval = pchisq(chiSqG, df = 1, lower.tail = F),
         chiSqGxE.pval = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaG, betaGxE, chiSqG.pval, chiSqGxE.pval)


results_ge_pvals <- results_ge %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         chiSqGE.pval = pchisq(chiSqGE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGE, chiSqGE.pval)


#--------------------------------------------------------#
# Implementing 2-step methods - Kooperberg
#--------------------------------------------------------#
kooperberg <- results_main_pvals %>% 
  filter(chiSqG.pval < 5e-8)

kooperberg_results <- results_main_pvals %>% 
  filter(SNP %in% kooperberg$SNP) %>% 
  mutate(kooperberg.bonferroni.pval = p.adjust(chiSqGxE.pval, method = 'bonferroni'))

manhattan(kooperberg_results, p = 'kooperberg.bonferroni.pval')



#--------------------------------------------------------#
# Implementing 2-step methods - Murcray
#--------------------------------------------------------#
murcray <- results_ge_pvals %>% 
  filter(chiSqGE.pval < 0.001)

0.05/7097

murcray_results <- results_main_pvals %>% 
  filter(SNP %in% murcray$SNP) %>% 
  mutate(murcray.bonferroni.pval = chiSqGxE.pval)
  mutate(murcray.bonferroni.pval = p.adjust(chiSqGxE.pval, method = 'bonferroni'))

manhattan(murcray_results, p = 'murcray.bonferroni.pval')




#--------------------------------------------------------#
# Implementing 2-step methods - EDGxE
# need to join results tables for this step (DG + GE)
#--------------------------------------------------------#

#------ subset testing ------#
edgxe <- inner_join(results_main[, c("SNP", "chiSqG")], results_ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(CC_DGGE = chiSqG + chiSqGE,
         CC_DGGE.pval.2df = pchisq(CC_DGGE, df = 2, lower.tail = F)) %>% 
  filter(CC_DGGE.pval.2df < 0.001)
  
edgxe_results <- results_main_pvals %>%
  filter(SNP %in% edgxe$SNP) %>% 
  mutate(edgxe.bonferroni.pval = chiSqGxE.pval)
  mutate(edgxe.bonferroni.pval = p.adjust(chiSqGxE.pval, method = 'bonferroni'))
  
manhattan(edgxe_results, p = 'edgxe.bonferroni.pval')


#--------------------------------------------------------#
# case-only
#--------------------------------------------------------#
kooperberg <- results_main_pvals %>% 
  filter(chiSqG.pval < 5e-8)

ChiSqCase

results_ge_pvals <- results_ge %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         ChiSqCase.pval = pchisq(ChiSqCase, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaCase, ChiSqCase.pval)

qq(results_ge_pvals$ChiSqCase.pval)
manhattan(results_ge_pvals, p = 'ChiSqCase.pval')


#------ weighted hypothesis testing ------#
# N : sample size of the study 
# Nbio: sample size of 'biological' information 
# m number of SNPs 
# ncp_stat: noncentrality parameter for a 2-df step 1 test statistic based on the statistical info 
# ncp_bio: noncentrality parameter for a 1-df step 1 'score' based on the functional biological info 
# ncp_GxE: noncentrality parameter for the 1-df step 1 GxE test statistic 
# sizeBin0 size of top/first bin 
# alpha: overall alpha level 

m = 5475351
sizeBin0 = 5
alpha=0.05
alpha1=0.001
# N=72438; Nbio=100; ncp_GxE=0

nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction 
nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} 
sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes 
endpointsBin = cumsum(sizeBin) # endpoints of the bins 
alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested 

subs = cumsum(sizeBin) <= alpha1*m 
alphaBin_subset = alpha / sizeBin * subs/ sum(subs==TRUE) 
alphaBin_linear = alpha / sizeBin / nbins 


#### suppose you have a data.frame of p values?
# pvals <- results %>%
#   mutate(Sdg_ge = zG^2 + z_GE^2,
#          Sdg_ge.pval = pchisq(Sdg_ge, df = 2, lower.tail = F),
#          zGxE.pval = 2*pnorm(-abs(zGxE))) %>% 
#   arrange(Sdg_ge.pval) %>% 
#   dplyr::select(SNP, CHR, BP, A1, A2, BetaGxE, Sdg_ge.pval, zGxE.pval)

edgxe <- inner_join(results_main[, c("SNP", "Chromosome", "Location", "Reference", "Alternate", "betaGxE", "chiSqG", "chiSqGxE")], results_ge[, c("SNP", "chiSqGE")], by = "SNP") %>% 
  mutate(CC_DGGE = chiSqG + chiSqGE,
         CC_DGGE.pval.2df = pchisq(CC_DGGE, df = 2, lower.tail = F)) %>% 
  arrange(CC_DGGE.pval.2df) %>% 
  dplyr::select(SNP, Chromosome, Location, Reference, Alternate, betaGxE, CC_DGGE.pval.2df, chiSqGxE)

pv <- data.table(edgxe)

## Function 11: Create Bin and significance level for weighted testing approach (needs to be sorted by p-value)
num = 5 # initial bin size
alpha = 0.05

rk.pv<-c(1:nrow(pv))
grp=ceiling(log(rk.pv/num+1,base=2))
pv[,Bin:=grp]
setkey(pv,Bin)
for(i in 1:max(grp))
{
  pv[J(i),Threshold:=alpha*2^(-i)/nrow(pv[J(i)])]
}


# table(pv$grp)
num = 21
scale = 0.7
min.p = 8
PlotAll = 0
cutoff = 1
last.sig = 7.63
last.sig = -log10(9.575697e-14)

# map info might be helper to sort x axis? probably BP? 
setnames(pv, c(8,9,10,3), c('p', 'grp','wt','MapInfo'))
pv[, y:= pv[,p]]

# create list where each component contains a bin? 
glist<-list()
for(i in 1:num){
  t <- pv[J(i)]
  t[, ref:=-1*log10(min(t[,wt]))]
  t[,x:=(t[,MapInfo]-min(t[,MapInfo]))/((max(t[,MapInfo])+0.0001)-min(t[,MapInfo]))*scale+i-1]
  glist[[i]]<-t
  rm(t)
}


## Scale mapinfo for each Bin to range between 0-1 and add a unit increase for successive Bin

color<-rep(c("blue","olivedrab4"),100)

plot(glist[[1]][,x], glist[[1]][,y], col="blue", xlab="Bin # for step1 p-values", ylab="-log10(step2 p-values)", xlim=c(0,num), ylim=c(0,min.p), axes=F,pch=19, cex=0.5)

lines(glist[[1]][,x], glist[[1]][,ref],col="black",lwd=2)

for(i in 2:num-1){
  points(glist[[i]][,x], glist[[i]][,y], col=color[i],pch=19, cex = 0.5)      
  lines(glist[[i]][,x], glist[[i]][,ref],col="black",lwd=2)
}

points(glist[[num]][,x], glist[[num]][,y], col=color[num],pch=19, cex = 0.5)
lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line


decix<-(num/2)-floor(num/2)
if(decix>0){
  axis(1,at=c(-1.5,seq(.5,num-0.5,2)), label=c(0,seq(1,num,2)),cex.axis=0.6)
} else {
  
  axis(1,at=seq(-.5,num-0.5,2), label=seq(0,num,2),cex.axis=0.5)
}

deciy<-(min.p/2)-floor(min.p/2)

axis(2,at=c(0:floor(min.p)),label=c(0:min.p),cex.axis=0.65)

title (main=title,sub="Bin Size=5")


