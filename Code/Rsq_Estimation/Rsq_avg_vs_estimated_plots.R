#=============================================================================#
# I want to compare John's estimated Rsq with Weighted Average
# and look at individual imputation Rsq squares compared to the estimates
#=============================================================================#

# start with chromosome 22.. maybe take a random sample
# and later use representative SNPs (low/high MAF, low/high Rsq, Typed/Imputed)
# maybe incorporate batch size as part of plot (size of point, symbol)

# let's start with a data.frame of weighted average Rsq I create previously 
# this is filtered @ MAF > 0.005
rsq_filter <- readRDS("~/data/HRC_InfoFile_Merged/HRC_WtAvgRsq_HRCRefPanelMAFs.rds")


# let's start with 10 random SNPs on chromosome 22
rsq_filter_chr22 <- filter(rsq_filter, grepl("^22:", ID)) %>% 
  # filter(mafbin == ">5.0%") %>% 
  filter(Rsq_avg > 0.3) %>% 
  slice(sample(nrow(.), 20))


snp_sample <- rsq_filter_chr22$ID
batch_size <- data.frame(batch = c("axiom_acs_aus_nf_Rsq", "axiom_mecc_cfr_ky_Rsq", "ccfr_1m_1mduo_reimpute_Rsq", "ccfr_omni_Rsq", "corect_oncoarray_Rsq", "corsa_axiom_Rsq", "cytosnp_comb_Rsq", "initial_comb_datasets_Rsq", "mecc_Rsq",  "newfoundland_omniquad_Rsq", "omni_comb_Rsq", "omniexpress_exomechip_Rsq", "oncoarray_to_usc_Rsq", "plco_3_Rsq", "reach_Rsq", "ukbiobank_Rsq"), 
                         sample_size = c(2766, 7501, 2178, 1600, 36621, 2467, 10908, 5924, 983, 637, 5986, 5503, 20912, 4864, 750, 27594))



# batch specific Rsq values
batch_rsq <- readRDS("~/data/HRC_InfoFile_Merged/HRC_InfoFile_Chr22.rds") 
snp_sample_index <- which(batch_rsq$ID %in% snp_sample)

# stick to indices
# add weighted average rsq
x <- batch_rsq[snp_sample_index, ] %>% 
  inner_join(rsq_filter_chr22[, c("ID", "Rsq_avg")], by = "ID")

# Rsq estimate from alt allele probability variance
# z <- readRDS("~/data/HRC_InfoFile_Merged/Rsq_Estimate_AltAlleleProbVar/FIGI_RsqEstimate_chr22.rds")[snp_sample_index]
z <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chr22.rds")[snp_sample_index]

x$Rsq_summary <- z


# to use ggplot, better to have Rsq values, then make data long

xlong <- x %>% 
  dplyr::select(ID, contains("Rsq")) %>% 
  gather(key = "batch", value = "Rsq", -ID, -Rsq_summary, -Rsq_avg) %>% 
  mutate(Rsq = as.numeric(Rsq), 
         myclr = ifelse(batch != "Rsq_Estimate", "black", "red")) %>%
  inner_join(batch_size, by = 'batch')


# add sample size for each imputation batch


ggplot() + 
  geom_point(aes(ID, Rsq, col = batch, size = sample_size), xlong) + 
  geom_point(data = x, aes(ID, Rsq_avg), shape = 2) + 
  geom_line(data = x, aes(ID, Rsq_avg, group = 1), color = "red") + 
  geom_point(data = x, aes(ID, Rsq_summary), shape = 3) + 
  geom_line(data = x, aes(ID, Rsq_summary, group = 1), color = "blue") + 
  theme_bw()
  


