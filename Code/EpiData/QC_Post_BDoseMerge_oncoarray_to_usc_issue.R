#====================================================
# oncoarray_to_usc issue:
# chromosomes don't have 100% matching sample names 
#
# find set of samples that ARE in all chromosomes
# identify samples that are not in all chromosomes
#
#====================================================

#system("scp hpcc:/staging/dvc/andreeki/bdose/qc/*.txt ~/Desktop/work")
rm(list = ls())
library(tidyverse)
library(data.table)

setwd("~/Desktop/work/")

# do simple count of unique IDs across chromosomes, compare to chromosome 1 for now
wtf <- function(x) {
	tmp_chr <- fread(paste0(x, "_chr1_samplelist.txt"), sep = '\t', skip = 1, header = F, col.names = "IID")
 	tmp_all <- do.call(rbind, lapply(list.files( full.names = T, pattern = paste0(x, "_chr")), fread, stringsAsFactors = F, sep = '\t', skip = 1, header = F, col.names = "IID"))
	print(length(tmp_chr$IID))
	print(length(unique(tmp_all$II)))
}

wtf("axiom_acs_aus_nf")
wtf("axiom_mecc_cfr_ky")
wtf("axiom_acs_aus_nf")
wtf("ccfr_1m_1mduo_reimpute")
wtf("ccfr_omni")
wtf("corect_oncoarray")
wtf("corect_oncoarray_nonEUR_reimpute")
wtf("corsa_axiom")
wtf("cytosnp_comb")
wtf("dachs3")
wtf("initial_comb_datasets")
wtf("mecc")
wtf("newfoundland_omniquad")
wtf("omni_comb")
wtf("omniexpress_exomechip")
wtf("oncoarray_to_usc") # ** had idea that oncoarray_to_usc has issues
wtf("plco_3")
wtf("reach")
wtf("ukbiobank")


#----------------------------------------------------
# read in oncoarray, inner-join to find IDs in common
# Useful, but didn't use here. 
# my.files <- list.files(full.names = F, pattern = "oncoarray_to_usc")
# for (i in 1:(length(my.files))) {
# 	cur.file <- fread(file = my.files[i], stringsAsFactors = F, sep = '\t', skip = 1, header = F, col.names = "IID")
# 	my.name <- my.files[i]
# 	assign(paste(my.name), cur.file)
# }

#-----------------------------------------------------
# another solution
oncoarray_to_usc <- lapply(list.files(pattern = "oncoarray_to_usc"), fread, stringsAsFactors = F, sep = '\t', skip = 1, header = F, col.names = "IID" )

common_samplenames <- Reduce(function(x,y) inner_join(x,y,by="IID",all=TRUE), oncoarray_to_usc)

problem_samples <- unique(do.call(rbind, lapply(oncoarray_to_usc, function(x) anti_join(x, common_samplenames, by = "IID"))))

# write out list of common samples (to use in the Covariates script...)
write.table(common_samplenames, file = "~/git/FIGI_EpiData/QC_oncoarray_to_usc_includelist.txt", quote = F, row.names = F)

# write out list of problem samples, ask Keith/Yi/Tabitha
write.table(problem_samples, file = "~/git/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt", quote = F, row.names = F)


# send over to rak for convenience
system("scp ~/git/FIGI_EpiData/QC_oncoarray_to_usc_includelist.txt rak:~/git/FIGI_EpiData")
system("scp ~/git/FIGI_EpiData/QC_oncoarray_to_usc_problems.txt rak:~/git/FIGI_EpiData")
