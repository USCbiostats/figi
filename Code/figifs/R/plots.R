#-----------------------------------------------------------------------------#
# QQ and Manhattan Plots
# Two-Step Weighted Hypothesis Test Plots
#-----------------------------------------------------------------------------#

#' calculate_pval
#'
#' Takes a vector of ChiSq statistics, output a vector of p values
#'
#' @param chiSq Vector of chiSq statistics - this is most applicable to GxEScanR results (e.g. a single column)
#' @param df ChiSq distribution degrees of freedom
#' @return Returns a vector of chi-square p values
#' @export
#'
#' @examples calculate_pval(3.6, 1)
#' ex <- rchisq(10, df=1)
#' calculate_pval(ex, df=1)
calculate_pval <- function(chiSq, df) {
  pchisq(chiSq, df = df, lower.tail = F)
}


#' write_plot_title
#'
#' Function to create figure titles based on the results statistic plotted
#'
#' @param stat GxEScan chi-square statistic (look up how to set/list parameter choices)
#' @param exposure Exposure variable (string)
#' @param covars Covariate string vector
#' @param N Sample size
#'
#' @return Plot title string
#' @export
#'
#' @examples write_plot_title('chiSqG', 'aspirin', c('sex', 'age'), 100000)
write_plot_title <- function(stat, exposure, covars, N) {
  gxescan_tests <- c(paste0("G Main Effects Results (N = ", N, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                     paste0("GxE Results (N = ", N, ")\noutc ~ G*", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("2DF Results (N = ", N, ")\noutc ~ G+G*", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("G|E Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("Case-Only Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("Control-Only Results (N = ", N, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                     paste0("3DF Results (N = ", N, ")\nchiSqG+chiSqGxE+chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGxE", "chiSq2df", "chiSqGE", "chiSqCase", "chiSqControl", "chiSq3df")
  return(gxescan_tests[stat])
}


#' create_qqplot
#'
#' Function creates a QQ plot and outputs *.png file (default location is a './figures' folder).
#'
#' @param dat Input data
#' @param exposure Exposure variable
#' @param stat GxEScan results chi-square statistic to plot
#' @param df Degrees of Freedom
#' @param filename_suffix For convenience when you're creating test plots
#'
#' @return Writes a .png file into location ./figures/
#' @export
#'
#' @examples create_qqplot(x, 'aspirin', c('sex', 'age'), 'chiSqGxE', df = 1, "_test")
create_qqplot <- function(dat, exposure, covars, stat, df, filename_suffix = "") {

  # calculate p value + lambda
  pvals <- calculate_pval(dat[, stat], df)
  lambda <- round( (median(qchisq(1-pvals, df)) / qchisq(0.5, df)), 4)

  # calculate lambda1000
  cases <- unique(dat[, 'Cases'])
  controls <- unique(dat[, 'Subjects']) - unique(dat[, 'Cases'])
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000)))

  # plotting function
  png(paste0("figures/QQ_Plot_", exposure, "_", stat, "_", paste0(covars, collapse = "_"), "_N_", total, filename_suffix, ".png"), height = 720, width = 1280)
  qqman::qq(pvals,
            xlab = "Expected -log10(p)",
            ylab = "Observed -log10(p)",
            main = write_plot_title(stat, exposure, covars, total),
            cex.main = 1.6,
            cex.axis = 1.3,
            cex.lab = 1.3,
            cex.sub = 1.3,
            cex = 1.4,
            pch = 1,
            col = 'blue4')
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3) # FYI ~~ adds spaces when using signif
  dev.off()
}





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
  fh_annotations <- fread("~/data/Annotations/crc_gwas_125k_indep_signals_95_EasyStrata_LDBased.tsv") %>%
    mutate(SNP = paste(Chr, Pos, sep = ":"))

  # format data for EasyStrata package, output to file
  dat <- dat %>%
    mutate(P = calculate_pval(dat[,stat], df),
           annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>%
    filter(!(P > 0.05 & annot == 0)) %>%
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
      --anumParMar 7;4;7;4
      --numDefaultSymbol 20
      --numDefaultCex 0.6
      --numCexAxis 1.2
      --numCexLab 1.2
      --arcdSymbolCrit P<5e-8
      --anumSymbol 17
      --arcdColourCrit P<5e-8
      --astrColour gray30
      --arcdCexCrit P<5e-8
      --anumCex 1.2

      STOP EASYX"), file = ecf_file_name, append = F)

  # run EasyStrata
  EasyStrata(ecf_file_name)
}




#-----------------------------------------------------------------------------#
# Weighted Hypothesis Testing for 2-step methods ------
# kooperberg, murcray, edge
#-----------------------------------------------------------------------------#

#' format_twostep_data
#'
#' Function to format step 1 (arrange by p value) --- assume data is the gxe object from GxEScanR
#'
#' @param dat Input data
#' @param stats_step1 Step1 chiSq statistic
#' @param sizeBin0 Size of initial bin (this value needs to be power optimized)
#' @param alpha Overall alpha level
#'
#' @return A DATA.TABLE (!) - mainly because the code I copied pasted used a dt. Contains bin number, bin logp threshold, normalized X axis information for plotting.
#' @export
#' @import data.table
#'
#' @examples format_data_twostep_data(gxe, 'chiSqG', 5, 0.05)


# dat <- gxe
# stats_step1 <- 'chiSqG'
# sizeBin0 = 5
# alpha = 0.05


format_twostep_data <- function(dat, stats_step1, sizeBin0, alpha) {

  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(dat, stats_step1, df=1) {
    data.table(dat)[, step1p := pchisq(dat[, eval(quote(stats_step1))], df = df, lower.tail = F)
      ][
      , step2p := pchisq(dat[,'chiSqGxE'],  df = 1, lower.tail = F)
      ][
      , y := -log10(step2p)
      ][
      order(step1p)
      ][
      , MapInfo := Location
      ]
  }

  if(stats_step1 == 'chiSqEDGE') {
    pv <- create_pval_info(dat, stats_step1, df = 2)
  } else {
    pv <- create_pval_info(dat, stats_step1, df = 1)
  }


  # format output for plotting..
  m = nrow(pv)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1}
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha elevel at which a SNP landing in bin 1,2,...nbins is tested

  # add grp, wt (p value threshold for each bin), rename p value to 'y'
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2)) # math shortcut to put bin size to every single marker by c(1:nrow(pv))
  pv[,grp:=grp] # assigning group to the p values..
  setkey(pv,grp)
  for(i in 1:max(grp))
  {
    pv[J(i),wt:=alpha*2^(-i)/nrow(pv[J(i)])] # data.table syntax, create threshold value
  }

  # return the data.table
  return(pv)
}


#' create_twostep_weighted_plot
#'
#' Create weighted hypothesis testing plot. First step is creating a list object based on the data.table object output from "format_data_twostep_data" function. Second step is actual plotting. Remember that it only plots the first 10 bins for legibility.
#'
#' @param dat Input data
#' @param sizeBin0 Initial bin size (seeds size of remaining bins)
#' @param alpha Overall alpha value (default should be 0.05)
#' @param binsToPlot Number of bins to include in plot
#' @param statistic Step 1 filtering chi-square statistic (depends on method desired)
#'
#' @return A weighted hypothesis plot (png file)
#' @export
#'
#' @examples create_twostep_weighted_plot(gxe_twostep, exposure = plot_exposure, covars = plot_covariates, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')


# dat = gxe_twostep
# exposure = plot_exposure
# covars = plot_covariates
# sizeBin0 = 5
# alpha = 0.05
# binsToPlot = 10
# statistic = 'chiSqG'
# filename_suffix = ""


create_twostep_weighted_plot <- function(dat, exposure, covars, sizeBin0, alpha, binsToPlot, statistic, filename_suffix = "") {

  cases <- unique(data.frame(dat[, 'Cases']))
  controls <- unique(data.frame(dat[, 'Subjects'])) - unique(data.frame(dat[, 'Cases']))
  total <- cases + controls

  # some plot title and file names
  write_twostep_weightedHT_plot_title <- function(statistic, exposure, covars, total) {
    gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", total, ")\noutc ~ G+", paste0(covars, collapse = "+"),"+", exposure),
                       paste0("G|E 2-step Procedure Results (N = ", total, ")\nG ~ ", exposure, "+", paste0(covars, collapse = "+")),
                       paste0("EDGE 2-step Procedure Results (N = ", total, ")\nchiSqG + chiSqGE"))
    names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
    return(gxescan_tests[statistic])
  }

  # write_twostep_weightedHT_plot_title("chiSqG", exposure, covars, total)

  # bin information
  m = nrow(dat)
  nbins = floor(log2(m/sizeBin0 + 1)) # number of bins for second-step weighted Bonferroni correction
  nbins = if (m > sizeBin0 * 2^nbins) {nbins = nbins + 1} # add +1 bin if condition met
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) ) # bin sizes
  endpointsBin = cumsum(sizeBin) # endpoints of the bins
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin # alpha level for each bbin 1,2, ... N bin tested

  # create list and vars for plotting
  min.p = 12 # this might be plot upper limit in -log10 scale, not sure why called 'min.p'
  last.sig = alphaBin[binsToPlot]

  # create list where each component contains a bin
  # log transform bin alpha value
  # create 'x', normalized position information for each bin
  glist<-list()
  for(i in 1:binsToPlot){
    t <- dat[J(i)]
    t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
    t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
    glist[[i]]<-t
    rm(t)
  }

  # CREATE PLOT
  head(glist[[1]]) # for reference

  png(paste0("figures/TwoStep_WeightedHypothesis_", deparse(substitute(dat)), "_", exposure, "_", statistic, "_", paste0(covars, collapse = "_"), "_N_", total, filename_suffix, ".png"), height = 720, width = 1280)
  # png(write_weighted_plot_filename(deparse(substitute(dat)), statistic), width = 1280, height = 720)
  color <- rep(c("blue","olivedrab4"),100)
  plot(glist[[1]][,x], glist[[1]][,y],
       col = "blue1",
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0,binsToPlot),
       ylim=c(0,min.p),
       axes=F, pch=19,
       cex.main = 1.6,
       cex.axis = 1.3,
       cex.lab = 1.3,
       cex.sub = 1.3,
       cex = 1)
  lines(glist[[1]][,x], glist[[1]][,ref],
        col = "black",lwd=2)

  # the rest of the points ...  =|
  # (adding to current plot..)
  for(i in 2:binsToPlot){
    points(glist[[i]][,x], glist[[i]][,y],
           col = color[i], pch = 19,
           cex.main = 1.6,
           cex.axis = 1.3,
           cex.lab = 1.3,
           cex.sub = 1.3,
           cex = 1)
    lines(glist[[i]][,x], glist[[i]][,ref],
          col = "black",lwd = 2)
  }

  ## the last bin..
  ## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
  ## let's only plot the first 15 bins for now, so change code a bit above.
  # points(glist[[num]][,x], glist[[num]][,y],
  #        col= color[num], pch = 19, cex = 0.5)
  # lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line

  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=1)
  title(main = write_twostep_weightedHT_plot_title(statistic, exposure, covars, total), sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.6, cex.sub = 1.3)

  dev.off()

}
