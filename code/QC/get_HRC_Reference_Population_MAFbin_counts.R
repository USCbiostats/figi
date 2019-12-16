#=============================================================================#
# Summarizing minor allele frequencies for HRC reference panel
#
# Using HRC reference population MAFs, with the assumption that it 
# does a good job approximating the allele frequencies of the 125k from FIGI
#=============================================================================#

# take HRC tab file, separate by chromosome (why not)
# #!/bin/bash
# for chr in {1..21}
# do
# awk '{if ($1 == '"${chr}"') {print $1, $2, $4, $5, $8} }' HRC.r1-1.GRCh37.wgs.mac5.sites.tab > ~/HRC.r1-1.GRCh37.chr${chr}.txt
# done

# read in HRC reference panel info
setwd("~/data/HRC")

hrc <- do.call(rbind, lapply(list.files(pattern = "HRC.r1-1", full.names = T), fread, stringsAsFactors=F)) %>% 
	mutate(id = paste0(V1, ":", V2, "_", V3, "_", V4),
				 maf = 0.5 - abs(V5 - 0.5),
				 mafbin = cut(maf, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("<0.1%", "0.1%-1.0%", "1.0%-5.0%", ">5.0%")))
saveRDS(hrc, file = "HRC.r1-1.rds")

table(hrc$mafbin)
