#========================================================#
# The problem:
#
# some markers have very different AFs 
# that vary by imputation batch / platform
#
#========================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(broom)
library(reshape)
library(R.utils)
library(lmtest)
rm(list = ls())

# GxEScanR results - asp_ref
results <- do.call(rbind, lapply(list.files(path = "./results/", full.names = T, pattern = "results_asp_ref_sex_age_pcs_studyname_N72438_"), fread, stringsAsFactors = F))
names(results)


results_GxE <- results %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = F) %>% 
  mutate(CHR = as.numeric(CHR), 
         BP = as.numeric(BP), 
         P = pchisq(chiSqGxE, df = 1, lower.tail = F)) %>%
  dplyr::select(SNP, CHR, BP, Subjects, Cases, Reference, Alternate, betaGxE, P) %>% 
  filter(P < 5e-8)




# covariate file
load("~/git/DATA/FIGI_EpiData_rdata/FIGI_Genotype_Epi_181120.RData")
pc30k <- fread("~/Dropbox/code/FIGI_Results/PrincipalComponents/FIGI_GxESet.eigenvec", skip = 1, 
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20))))
cov <- Epi %>%
  filter(drop == 0 & gxe == 1) %>% 
  inner_join(pc30k, by = c('vcfid' = 'IID')) %>% 
  mutate(outcome = ifelse(outc == "Case", 1, 0),
         age_ref_imp = as.numeric(age_ref_imp),
         sex = ifelse(sex == "Female", 0, 1),
         studyname = ifelse(studyname %in% c("NFCCR_1", "NFCCR_2"), "NFCCR", studyname),
         studyname = ifelse(studyname == "Colo2&3", "Colo23", studyname),
         asp_ref = ifelse(asp_ref == "Yes", 1, ifelse(asp_ref == "No", 0, NA))) %>% 
  mutate(V2 = ifelse(V2 == "ccfr_1m_1mduo", "ccfr_1m_1mduo_reimpute", V2)) %>% 
  dplyr::select(vcfid, outcome, age_ref_imp, sex, studyname, paste0(rep("PC", 10), seq(1,10)), asp_ref, V2, platform) %>% 
  filter(complete.cases(.))

exclude_studies <- c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT") # case-only studies.. 
cov <- filter(cov, !studyname %in% exclude_studies)



# First order of business, let's take some examples of markers that have very different AFs by study/platform
# I happened upon because of the GxE results + random top hits on Chromosome 2, so let's start with those examples

# start with GxE 'hits'
snps <- dplyr::select(cov, vcfid)
for(chr in c(2,3,6,7,10, 16, 17, 18, 21, 22)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/Dropbox/code/FIGI_Results/aspref_sex_age_pc10_studyname/extract/GxEScanR_GxE_sig_loci_extract_chr", chr, ".rds"))) %>% rownames_to_column(var = 'vcfid')
  snps <- inner_join(snps, x, by = 'vcfid')
}
df <- inner_join(cov, snps, by = "vcfid")

snps <- dplyr::select(cov, vcfid)
for(chr in c(2)) {
  x <- as.data.frame(readRDS(paste0("/home/rak/git/Dosage_Extract/GxEScanR_asp_ref_extract_chr", chr, ".rds"))) %>% rownames_to_column(var = 'vcfid')
  snps <- inner_join(snps, x, by = 'vcfid') %>% 
    dplyr::select(vcfid, "2:96780986", "2:131698423")
}
df <- inner_join(df, snps, by = "vcfid")

# take a look @ what I mean - mean AFs for each of those snps.. 
af_studyname <- df %>% 
  group_by(studyname) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008", "2:96780986", "2:131698423"), funs( sum(.) / n() / 2))
saveRDS(af_studyname, './gxe_problems/af_studyname.rds')


ggplot(data = af_studyname, aes(studyname, `2:4964187`)) +
  geom_point() + 
  theme_bw()

plot(af_studyname[, c("2:4964187"), drop = T], cex = 0.8)

plot(af_studyname$`2:4964187`, cex = 0.8)

plot(af_studyname$`3:98981063`)
plot(af_studyname$`6:138017298`)
plot(af_studyname$`7:99995536`)
plot(af_studyname$`10:3001065`)
plot(af_studyname$`16:79591072`)
plot(af_studyname$`17:18528708`)
plot(af_studyname$`18:3780526`)
plot(af_studyname$`21:34169318`)
plot(af_studyname$`22:40651008`)

plot(af_studyname$`2:96780986`)
plot(af_studyname$`2:131698423`)

# the last two are interesting
af_platform <- df %>% 
  group_by(platform) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008", "2:96780986", "2:131698423"), funs( sum(.) / n() / 2))

plot(af_platform$`2:96780986`)
plot(af_platform$`2:131698423`)

# by V2 (imputation batch)
af_impBatch <- df %>% 
  group_by(V2) %>% 
  summarise_at(c("2:4964187", "3:98981063",  "6:138017298", "7:99995536",  "10:3001065",  "16:79591072", "17:18528708", "18:3780526",  "21:34169318", "22:40651008", "2:96780986", "2:131698423"), funs( sum(.) / n() / 2))
               
plot(af_impBatch$`2:96780986`)
plot(af_impBatch$`2:131698423`)



# the next step is to get Rsq information (if available...)
# currently don't have that information for: dachs3, omniexpress_exomechip, oncoarray_to_usc, ukbiobank. 

batchList <- as.character(expression(axiom_acs_aus_nf  ,      corect_oncoarray    ,              dachs3   ,              omni_comb  ,            reach,
                                     axiom_mecc_cfr_ky   ,    corect_oncoarray_nonEUR_reimpute,  initial_comb_datasets,  omniexpress_exomechip , ukbiobank,
                                     ccfr_1m_1mduo_reimpute,  corsa_axiom      ,                 mecc      ,             oncoarray_to_usc,
                                     ccfr_omni           ,    cytosnp_comb   ,                   newfoundland_omniquad  ,plco_3))

batchList <- as.character(expression(axiom_acs_aus_nf  ,      corect_oncoarray    ,                  omni_comb  ,            reach,
                                     axiom_mecc_cfr_ky   ,    corect_oncoarray_nonEUR_reimpute,  initial_comb_datasets,  
                                     ccfr_1m_1mduo_reimpute,  corsa_axiom      ,                 mecc      ,          
                                     ccfr_omni           ,    cytosnp_comb   ,                   newfoundland_omniquad  ,plco_3))
colNames <- as.character(expression(SNP  ,REF  ,ALT , ALT_Frq, MAF,     AvgCall ,Rsq    , Genotyped,       LooRsq,  EmpR  ,  EmpRsq,  Dose0 ,  Dose1))
batch <- as.list(batchList)
kwik <- function(b, bp, chr) as.data.frame(matrix(system(paste0('zcat /home/rak/git/HRC_InfoFiles/', b, '/chr', chr, '.info.gz | grep ', bp), intern = T)))

kwik_wrapper <- function(bp, chr) {
  do.call(rbind, lapply(batch, kwik, bp = bp, chr = chr)) %>% 
    separate(V1, into = c(colNames), sep = '\t') %>% 
    mutate(batch = rep(batchList, each = 2))
}

kwik_wrapper <- function(bp, chr) {
  do.call(rbind, lapply(batch, kwik, bp = bp, chr = chr)) %>% 
    separate(V1, into = c(colNames), sep = '\t') %>% 
    mutate(batch = batchList)
}



# varies by platform.
chr2_96780986 <- kwik_wrapper(chr = 2, bp = 96780986)

  # hold on - make sure MAFs match? for those batches you have info for... 
  # yeah they're pretty close (which makes sense, but good sanity check, don't trust anything anymore)
z <- inner_join(af_impBatch[ ,c('V2', "2:96780986")], chr2_96780986, by = c("V2" = "batch"))


# varies by platform.. 
chr2_131698423 <- kwik_wrapper(chr = 2, bp = 131698423)

chr2_4964187 <- kwik_wrapper(chr = 2, bp = 4964187)
chr3_98981063 <- kwik_wrapper(chr = 3, bp = 98981063)
chr6_138017298 <- kwik_wrapper(chr = 6, bp = 138017298)
chr7_99995536 <- kwik_wrapper(chr = 7, bp = 99995536)
chr10_3001065 <- kwik_wrapper(chr = 10, bp = 3001065) #* needs cleaning
chr16_79591072 <- kwik_wrapper(chr = 16, bp = 79591072)
chr17_18528708 <- kwik_wrapper(chr = 17, bp = 18528708)
chr18_3780526 <- kwik_wrapper(chr = 18, bp = 3780526) #* needs cleaning
chr21_34169318 <- kwik_wrapper(chr = 21, bp = 34169318)
chr22_40651008 <- kwik_wrapper(chr = 22, bp = 40651008)


# results.. 
# chr2_96780986 - poorly imputed in most platforms, except when directly genotyped. 
# chr2_131698423 - same story as above

# chr2_4964187 - explore further using glm (with and without UKBIOBANK)
# chr3_98981063 - mostly poorly imputed
# chr6_138017298 - Rsq varies widely
# chr7_99995536 - Rsq varies widely, poor
# chr10_3001065 - explore further.. most of the platforms have good Rsq, with some exceptions
# chr16_79591072 - poorly imputed, though some are good.. 
# chr17_18528708 - poor
# chr18_3780526 - mediocre imputation
# chr21_34169318 - mostly shit imputation
# chr22_40651008 - same as above

# so most of these markers, first of all they're all pretty rare. secondly, many have very poor imputation quality. So without UKBIOBANK messing the dosage values up, none of them would have made it.

# chr2_4964187 - this one is a bit interesting
# confirm we get same results (yes)
df_glm <- df %>% 
  mutate(studyname, factor(studyname))

model_base <- glm(outcome ~ `2:4964187` + asp_ref + age_ref_imp + sex + studyname + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7+ PC8 + PC9 + PC10, data = df_glm, family = 'binomial')

model01 <- glm(outcome ~ `2:4964187`*asp_ref + age_ref_imp + sex + studyname + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7+ PC8 + PC9 + PC10, data = df_glm, family = 'binomial')
coef(summary(model01))
tidy(model01)


lrtest(model_base, model01)



# perform GLM without ukbiobank
df_glm <- df %>% 
  mutate(studyname, factor(studyname)) %>% 
  filter(studyname != "UKB_1")

model_base <- glm(outcome ~ `2:4964187` + asp_ref + age_ref_imp + sex + studyname + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7+ PC8 + PC9 + PC10, data = df_glm, family = 'binomial')

model01 <- glm(outcome ~ `2:4964187`*asp_ref + age_ref_imp + sex + studyname + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7+ PC8 + PC9 + PC10, data = df_glm, family = 'binomial')
coef(summary(model01))

lrtest(model_base, model01) # significance is gone. I bet this is the case for every single marker.. 



# run GLM for all 10 markers, compare before and after estimates (with and without UKB...)
df_glm <- df %>% 
  mutate(studyname, factor(studyname))

glm_base <- function(y) glm(outcome ~ y + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + studyname, data = df_glm, family = binomial)
glm_model <- function(y) glm(outcome ~ y * asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + studyname, data = df_glm, family = binomial)
glm_base_results <- map(df_glm[,19:28], function(x) glm_base(x))
glm_model_results <- map(df_glm[,19:28], function(x) glm_model(x))

# run lrtest from lmtest package.
run_lrtest <- function(x,y) lrtest(x, y)
results_gxe <- mapply(run_lrtest, glm_base_results, glm_model_results, SIMPLIFY = F)


#------------------------- no ukb ---------------------_#
df_glm <- df %>% 
  mutate(studyname, factor(studyname)) %>% 
  filter(studyname != "UKB_1")

glm_base_results <- map(df_glm[,19:28], function(x) glm_base(x))
glm_model_results <- map(df_glm[,19:28], function(x) glm_model(x))

# run lrtest from lmtest package.
run_lrtest <- function(x,y) lrtest(x, y)
results_gxe_noUKB <- mapply(run_lrtest, glm_base_results, glm_model_results, SIMPLIFY = F)

results_gxe_noUKB_clean <- do.call(rbind, results_gxe_noUKB) %>% 
  tibble::rownames_to_column('SNP') %>% 
  filter(!is.na(Chisq)) %>% 
  mutate(SNP = gsub('.{2}$', '', SNP))

z <- inner_join(results_GxE, results_gxe_noUKB_clean, by = 'SNP') %>% 
  mutate(bGxE_gxescan = betaGxE, 
         P_gxescan = P, 
         P_glm_noUKB = `Pr(>Chisq)`) %>% 
  dplyr::select(SNP, Subjects, Cases, Reference, Alternate, bGxE_gxescan, P_gxescan, P_glm_noUKB)

saveRDS(z, file = "./gxe_problems/glm_noUKB.rds")









##### Why did UKBIOBANK cause this 

af_studyname_asp_ref <- df %>% 
  group_by(studyname, outcome, asp_ref) %>% 
  summarise_at(c("2:4964187"), function(x) mean(x))
