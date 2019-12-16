library(data.table)
library(tidyverse)
library(kableExtra)
library(R.utils)

options(scipen=999)

quick_wrap = function(x) {
  z <- fread(paste0("~/data/HRC_InfoFile/ukbiobank/chr", x, ".info.gz"))
  
  y <- z %>% 
    filter(Genotyped != "Typed_Only") %>% 
    separate(SNP, into = c("CHROM", "POS"), remove = F) %>% 
    mutate(ALT_Frq = format(as.numeric(ALT_Frq), nsmall=5), 
           MAF = format(as.numeric(MAF), nsmall=5), 
           Rsq = format(as.numeric(Rsq), nsmall=5), 
           EmpRsq = format(as.numeric(EmpRsq), nsmall=5)) %>% 
    mutate(INFO = case_when(Genotyped == "Imputed" ~ paste0("AF=", ALT_Frq, ";", "MAF=", MAF, ";", "R2=", Rsq, ";", "IMPUTED"),
                            Genotyped == "Genotyped" ~ paste0("AF=", ALT_Frq, ";", "MAF=", MAF, ";", "R2=", Rsq, ";", "ER2=", EmpRsq, ";", "TYPED")),
           QUAL = ".", 
           FILTER = "PASS") %>% 
    dplyr::select(CHROM, POS, SNP, `REF(0)`, `ALT(1)`, QUAL, FILTER, INFO)
  
  writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), collapse = "\t"), con = paste0("~/chr", x, ".vcf.cut"))
  write.table(y, file = paste0("~/chr", x, ".vcf.cut"), quote = F, sep = "\t", row.names = F, col.names = F, append = T)
  
}


for(chr in 1:22) {
  quick_wrap(chr)
}




## let's see if we can run IC on a subset of markers
# (in progress)
z <- fread(paste0("~/data/HRC_InfoFile/axiom_acs_aus_nf/chr", x, ".info.gz"))

y <- z %>% 
  filter(Genotyped != "Typed_Only") %>% 
  separate(SNP, into = c("CHROM", "POS"), remove = F) %>% 
  mutate(ALT_Frq = format(as.numeric(ALT_Frq), nsmall=5), 
         MAF = format(as.numeric(MAF), nsmall=5), 
         Rsq = format(as.numeric(Rsq), nsmall=5), 
         EmpRsq = format(as.numeric(EmpRsq), nsmall=5)) %>% 
  mutate(INFO = case_when(Genotyped == "Imputed" ~ paste0("AF=", ALT_Frq, ";", "MAF=", MAF, ";", "R2=", Rsq, ";", "IMPUTED"),
                          Genotyped == "Genotyped" ~ paste0("AF=", ALT_Frq, ";", "MAF=", MAF, ";", "R2=", Rsq, ";", "ER2=", EmpRsq, ";", "TYPED")),
         QUAL = ".", 
         FILTER = "PASS") %>% 
  dplyr::select(CHROM, POS, SNP, `REF(0)`, `ALT(1)`, QUAL, FILTER, INFO)

writeLines(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), collapse = "\t"), con = paste0("~/chr", x, ".vcf.cut"))
write.table(y, file = paste0("~/chr", x, ".vcf.cut"), quote = F, sep = "\t", row.names = F, col.names = F, append = T)

