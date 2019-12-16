#=============================================================================#
# process IC (postimputation QC) results
#=============================================================================#
library(tidyverse)
library(data.table)
options(scipen=10000)

# ------ INFO score histograms ------
# loop over each batch

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

for(batch in batch_list) {

  x <- fread(paste0("files/INFO-", batch, ".GW.txt")) %>% 
    filter(V1 != "Total")

  ggplot(data = x) + 
    geom_histogram(aes(V1,V3), stat='identity') +
    theme_bw() + 
    labs(title = paste0("Rsq - ", batch), x = "Rsq", y = "Percent") + 
    ylim(0, 60)
  ggsave(filename = paste0("figures/INFO_histogram_", batch, ".png"), width = 4, height = 3)
}
