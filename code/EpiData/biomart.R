library(biomaRt)

rsid=c("rs123","rs150")
library(biomaRt)
ensembl_snp=useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

getBM(attributes=c("refsnp_source",'refsnp_id','chr_name','chrom_start','chrom_end',"consequence_type_tv", "clinical_significance"), filters = 'snp_filter', values = rsid,mart = ensembl_snp)




# # add chr, start and end position of the gene
# ensembl = useMart("ensembl")
# ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
# gene.location = getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'),
# 											filters = c('ensembl_gene_id'),
# 											values = list(gene.ids$ensembl_id),
# 											mart = ensembl)
# saveRDS(gene.location, file = "tmp_ak/gene_location.rds");

mart = useMart('ensembl'), followed by listDatasets(mart)
listDatasets(mart)



grch37 = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")
getBM(attributes=c("refsnp_source",'refsnp_id','chr_name','chrom_start','chrom_end',"consequence_type_tv", "clinical_significance"), filters = 'snp_filter', values = rsid, mart = ensembl_snp)


asp <- readRDS("~/Dropbox/FIGI/Results/asp_ref/files/twostep_wht_chiSqG_asp_ref_age_ref_imp_sex_study_gxe_pc1_pc2_pc3_expectation_based.rds")
getBM(attributes=c("refsnp_source",'refsnp_id','chr_name','chrom_start','chrom_end',"consequence_type_tv", "clinical_significance"), filters = 'chromosomal_region', values = "5:40280202", mart = ensembl_snp)
