#=============================================================================#
# Rscript to create EasyStrata .ecf files
#
# make sure you define the variables in the script you're sourcing from
# (ecf1-ecf4)
# (not sure if best practices but move on)
#=============================================================================#

# ecf1 <- "/home/rak/Dropbox/FIGI/Results/ASPIRI/figures"
# ecf2 <- "ID;CHR;BP;A1;A2;betaG;P"
# ecf3 <- "character;numeric;numeric;character;character;numeric;numeric"
# ecf4 <- "/media/work/tmp/EasyStrata_MarginalG_folate_dietqc2_sex_age_pc_energytot_studygxe_N52447.txt"

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
--astrDefaultColourChr gray51;gray66
--blnYAxisBreak 1
--numYAxisBreak 22
# Annotation
--fileAnnot /home/rak/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv
--numAnnotPosLim 1
# Horizontal lines
--anumAddPvalLine 5e-8
--anumAddPvalLineLty 6
--astrAddPvalLineCol coral3
# Other Graphical Params
--numDefaultSymbol 20
--numDefaultCex 0.1
--numCexAxis 1.5
--numCexLab 2
--arcdSymbolCrit P<5e-8
--anumSymbol 17
--arcdColourCrit P<5e-8
--astrColour gray30
--arcdCexCrit P<5e-8
--anumCex 0.6

STOP EASYX"), file = ecf_file_name, append = F)