#=============================================================================#
# Rscript to create EasyStrata .ecf files
#
# make sure you define the variables in the script you're sourcing from
# (ecf1-ecf4)
# (not sure if best practices but move on)
#
# Script specific for miami plot
# for convenience just make sure 
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
## Miami Plot 
################

MIAMIPLOT
--colMIAMIPlotUp P1
--colMIAMIPlotDown P2
--colInChr CHR
--colInPos BP
--numPvalOffset 0.05
--numWidth 1280
--numHeight 720	
--blnYAxisBreak 1
--numYAxisBreak 22
## Annotations
--fileAnnot /home/rak/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv
--numAnnotPosLim 1
## Colors + axis breaks
#--astrDefaultColourChrUp dodgerblue3;dodgerblue4
#--astrDefaultColourChrDown dodgerblue3;dodgerblue4
## horizontal line
--anumAddPvalLine 5e-8
--astrAddPvalLineCol coral3
--anumAddPvalLineLty 6
## others
--numDefaultSymbol 20
--numDefaultCex 0.1
--numCexAxis 1.5
--numCexLab 2
--arcdSymbolCritUp P1 <5e-8
--arcdSymbolCritDown P2 <5e-8
--anumSymbolUp 17
--anumSymbolDown 17
--arcdColourCritUp P1 <5e-8
--arcdColourCritDown P2 <5e-8
--astrColourUp gray30
--astrColourDown gray30
--arcdCexCritUp P1 <5e-8
--arcdCexCritDown P2 <5e-8
--anumCexUp 0.6
--anumCexDown 0.6

STOP EASYX"), file = ecf_file_name, append = F)