


#' create_manhattanplot
#'
#' Function to create manhattanplot. Uses package EasyStrata - need to write 2 files: results flat text (with appropriate columns), and .ecf file which contains EasyStrata parameters.
#'
#' NOTE: using hardcoded paths to store temporary results text files, and to read annotation file. Be careful!
#'
#' @param dat Input data
#' @param exposure Exposure variable (string)
#' @param covars Covariate string vector
#' @param stat GxEScan results chi-square statistic
#' @param df Degrees of freedom
#' @param filename_suffix For convenience when you're creating test plots
#'
#' @return Writes a .png file into location ./figures/
#' @export
#'
#' @examples create_manhattanplot(x, 'aspirin', c('sex', 'age'), 'chiSqGxE', df = 1, "_test")
create_manhattanplot <- function(dat, exposure, covars, stat, df, filename_suffix = "") {
  
  # get case/control/total counts (if file structure changes... bad)
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])
  total <- cases + controls
  
  # annotation (only used to make sure plot retains markers p > 0.05 but in LD with gwas hit)
  # probably bad practice
  # fh_annotations <- fread("~/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv") %>%
  #   mutate(SNP = paste(Chr, Pos, sep = ":"))
  
  # format data for EasyStrata package, output to file
  # dat <- dat %>%
  #   mutate(P = calculate_pval(dat[,stat], df),
  #          annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>%
  #   filter(!(P > 0.05 & annot == 0)) %>%
  #   dplyr::rename(CHR = Chromosome,
  #                 BP = Location) %>%
  #   dplyr::select(ID, CHR, BP, P)
  
  dat <- dat %>%
    mutate(P = calculate_pval(dat[,stat], df)) %>%
    filter(!(P > 0.05)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(ID, CHR, BP, P)
  
  
  write.table(dat, file = paste0("/media/work/tmp/EasyStrata_", exposure, "_", stat, "_", paste0(covars, collapse = "_"), "_N_", total, filename_suffix, ".txt"), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(getwd(), "/figures")
  ecf2 <- "ID;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- paste0("/media/work/tmp/EasyStrata_", exposure, "_", stat, "_", paste0(covars, collapse = "_"), "_N_", total, filename_suffix, ".txt")
  ecf_file_name <- paste0("files/EasyStrata_", exposure, "_", stat, "_", paste0(covars, collapse = "_"), "_N_", total, filename_suffix, ".ecf")
  
  cat(paste0("DEFINE	--pathOut ", ecf1, "
      --acolIn ", ecf2, "
      --acolInClasses ", ecf3, "
      --strMissing NA
      --strSeparator TAB

      EASYIN	--fileIn ", ecf4, "

      START EASYX

      ################
      ## MHplot
      ################

      MHPLOT
      --colMHPlot P
      --colInChr CHR
      --colInPos BP
      --numWidth 1280
      --numHeight 720
      --astrDefaultColourChr gray70;gray80
      --blnYAxisBreak 1
      --numYAxisBreak 22
      # Annotation
      #--fileAnnot /home/rak/data/Annotations/crc_gwas_indep_signals_140_EasyStrata_LDBased.tsv
      #--fileAnnot /home/rak/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv
      --fileAnnot /home/rak/data/Annotations/temp_annotation_ver2.txt
      --numAnnotPosLim 1
      # Horizontal lines
      --anumAddPvalLine 5e-8
      --anumAddPvalLineLty 6
      --astrAddPvalLineCol coral3
      # Other Graphical Params
      --anumParMar 7;4;7;4
      --numDefaultSymbol 20
      --numDefaultCex 0.6
      --numCexAxis 1.2
      --numCexLab 1.2
      --arcdSymbolCrit P<5e-8
      --anumSymbol 17
      --arcdColourCrit P<5e-8
      --astrColour grey31
      --arcdCexCrit P<5e-8
      --anumCex 1.2

      STOP EASYX"), file = ecf_file_name, append = F)
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}



