#=============================================================================#
# GxEScan post-hoc analysis
# 06/01/2019
# 
# Commonly used functions to generate plots etc
# (I know you tried this before, but this time it might be worth it)
# (make sure john's output absolutely never changes again)
# BAD PRACTICE - RELIES ON GLOBAL ENVIRONMENTS: 
# - E: exposure (character)
# - covs: covariates (vector)
# - N: samplesize
#
# keep in mind that for now, you're only creating a set number of plots:
# - G, GxE, 2DF, 3DF, GE, Case, Control. don't go nuts with overcomplication
#=============================================================================#
library(tidyverse)
library(data.table)
library(EasyStrata)
library(qqman)


#-----------------------------------------------------------------------------#
# calculate lambdas ------
# quick function to calculate lambdas and lambda1000
# make sure you use the appropriate DF for the result
# variables cases/controls/cases1000/controls1000 
#-----------------------------------------------------------------------------#
# getlambda <- function(pvals) {
#   chisq <- qchisq(1-pvals, 1)
#   lambda <- round(median(chisq)/qchisq(0.5,1),4)
# }
# 
# getlambda2df <- function(pvals) {
#   chisq <- qchisq(1-pvals, 2)
#   lambda <- round(median(chisq)/qchisq(0.5,2),4)
# }
# 
# getlambda3df <- function(pvals) {
#   chisq <- qchisq(1-pvals, 3)
#   lambda <- round(median(chisq)/qchisq(0.5,3),4)
# }
# 
# # takes gxescan results data + lambda calculated above
# getlambda1000 <- function(data, lambda) {
#   cases <- unique(data[, 'Cases'])
#   controls <- unique(data[, 'Subjects']) - unique(data[, 'Cases'])
#   total <- cases + controls
#   cases1000 <- (cases/total) * 1000
#   controls1000 <- (controls/total) * 1000
#   lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
# }




#-----------------------------------------------------------------------------#
# Results QQ and Manhattan Plots ------
# uses lambda functions above
# needs script called "~/Dropbox/FIGI/code/Functions/Rscript_create_ecf.R"
#-----------------------------------------------------------------------------#
# just add a variable called 'P'
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}

calculate_pval_miami <- function(data, statistic1, statistic2, df1, df2) {
  data$P1 <- pchisq(data[,statistic1], df = df1, lower.tail = F)
  data$P2 <- pchisq(data[,statistic2], df = df2, lower.tail = F)
  data
}



# dumb little helper functions for plot titles and file names
# NOTE SAVE LOCATION
# all arguments should be quote (Std evaluation...)
write_plot_filename <- function(data, statistic) {
  return(paste0("figures/QQ_Plot_", data, "_", statistic, "_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, ".png"))
}

write_easystrata_filename <- function(data, statistic) {
  return(paste0(paste("/media/work/tmp/EasyStrata", data, statistic, global_E, paste0(global_covs, collapse = "_"), "N", global_N, sep = "_"), ".txt"))
}

write_easystrata_filename_ecf <- function(data, statistic) {
  return(paste0(paste("files/EasyStrata", data, statistic, global_E, paste0(global_covs, collapse = "_"), "N", global_N, sep = "_"), ".ecf"))
}

write_plot_title <- function(statistic) {
  gxescan_tests <- c(paste0("G Main Effects Results (N = ", global_N, ")\noutc ~ G+", paste0(global_covs, collapse = "+"),"+", global_E), 
                     paste0("GxE Results (N = ", global_N, ")\noutc ~ G*", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("2DF Results (N = ", global_N, ")\noutc ~ G+G*", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("G|E Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("Case-Only Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("Control-Only Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("3DF Results (N = ", global_N, ")\nchiSqG+chiSqGxE+chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGxE", "chiSq2df", "chiSqGE", "chiSqCase", "chiSqControl", "chiSq3df")
  return(gxescan_tests[statistic])
}

write_plot_filename("gxe", "chiSqGxE")
write_easystrata_filename("gxe", "chiSqGxE")
write_easystrata_filename_ecf("gxe", "chiSqGxE")
# write_plot_title("chiSqG")
# write_plot_title("chiSqGxE")
# write_plot_title("chiSq2df")
# write_plot_title("chiSqGE")
# write_plot_title("chiSqCase")
# write_plot_title("chiSqControl")


# qq plot function
create_qqplot <- function(data, statistic, df) {

  # testing
  # statistic = 'chiSqGxE'
  # data = gxe
  # df = 1
  
  # calculate p value
  # calculate lambda
  tmpdata <- calculate_pval(data, statistic, df)
  lambda <- round( (median(qchisq(1-tmpdata[,'P'], df)) / qchisq(0.5, df)), 4)
  
  # calculate lambda1000
  cases <- unique(tmpdata[, 'Cases'])
  controls <- unique(tmpdata[, 'Subjects']) - unique(tmpdata[, 'Cases'])
  total <- cases + controls
  cases1000 <- (cases/total) * 1000
  controls1000 <- (controls/total) * 1000
  lambda1000 <- 1 + (lambda - 1) * ( (1/cases + 1/controls) / (1/(2*cases1000) + 1/(2*controls1000))) 
  
  # plotting function
  png(write_plot_filename(deparse(substitute(data)), statistic), height = 720, width = 1280)
  qqman::qq(tmpdata[, 'P'], 
            xlab = "Expected -log10(p)", 
            ylab = "Observed -log10(p)",
            main = write_plot_title(statistic),
            cex.main = 1.6, 
            cex.axis = 1.3, 
            cex.lab = 1.3,
            cex.sub = 1.3,
            col = 'blue4') 
  par(adj = 1)
  title(sub = bquote(lambda ~ '=' ~ .(signif(lambda, 4)) ~~ lambda[1000] ~ '=' ~.(signif(lambda1000, 4))), cex.sub = 1.3) # FYI ~~ adds spaces when using signif
  dev.off()
}


# Manhattan Plot
# work pending - modify annotion to match chr:pos:ref:alt instead of just chr:pos
create_manhattanplot <- function(data, statistic, df) {
  
  # testing
  # statistic = 'chiSqG'
  # data = gxe
  # df = 1
  
  # format data for easystrata (just be consistent)
  # remember 'calculate_pval' creates variable 'P'
  tmpdata <- calculate_pval(data, statistic, df)

  # need to write results to table. using same location for now (temporary directory in /media/work)
  data_easystrata <- tmpdata %>%
    mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>%
    filter(!(P > 0.05 & annot == 0)) %>%
    dplyr::rename(CHR = Chromosome,
                  BP = Location) %>%
    dplyr::select(ID, CHR, BP, P)
  write.table(data_easystrata, file = write_easystrata_filename(deparse(substitute(data)), statistic), quote = F, row.names = F, sep = '\t')

  # create ecf file
  ecf1 <- paste0(getwd(), "/figures")
  ecf2 <- "ID;CHR;BP;P"
  ecf3 <- "character;numeric;numeric;numeric"
  ecf4 <- write_easystrata_filename(deparse(substitute(data)), statistic)
  ecf_file_name <- write_easystrata_filename_ecf(deparse(substitute(data)), statistic)
  source("/home/rak/Dropbox/FIGI/Code/Functions/Rscript_create_ecf.R", local = T) # edit script as necessary

  # run EasyStrata
  EasyStrata(ecf_file_name)
}




# MIAMI Plot
# work pending - modify annotion to match chr:pos:ref:alt instead of just chr:pos
create_miamiplot <- function(data, statistic1, statistic2, df1, df2) {

  # format data for easystrata (just be consistent)
  # remember 'calculate_pval' creates variable 'P'
  tmpdata <- calculate_pval_miami(data, statistic1, statistic2, df1, df2)
  
  # need to write results to table. using same location for now (temporary directory in /media/work)
  data_easystrata <- tmpdata %>% 
    mutate(annot = ifelse(SNP %in% fh_annotations$SNP, 1, 0)) %>% 
    filter(!((P1 > 0.05 | P2 > 0.05) & annot == 0)) %>% 
    dplyr::rename(CHR = Chromosome, 
                  BP = Location) %>% 
    dplyr::select(ID, CHR, BP, P1, P2)
  write.table(data_easystrata, file = write_easystrata_filename(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_")), quote = F, row.names = F, sep = '\t')
  
  # create ecf file
  ecf1 <- paste0(getwd(), "/figures")
  ecf2 <- "ID;CHR;BP;P1;P2"
  ecf3 <- "character;numeric;numeric;numeric;numeric"
  ecf4 <- write_easystrata_filename(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_"))
  ecf_file_name <- write_easystrata_filename_ecf(deparse(substitute(data)), paste("MIAMI", statistic1, statistic2, sep = "_"))
  source("/home/rak/Dropbox/FIGI/Code/Functions/Rscript_create_ecf_miami.R", local = T) # edit script as necessary
  
  # run EasyStrata
  EasyStrata(ecf_file_name)
}




#-----------------------------------------------------------------------------#
# Weighted Hypothesis Testing for 2-step methods ------
# kooperberg, murcray, edge
#-----------------------------------------------------------------------------#

# variables required
# step1 statistic
# size of initial bin 
# number of SNPs
# overall alpha level = 0.05


# function to format step 1 (arrange by p value) --- assume data is the gxe object from GxEScanR
# output: data.table (because original code used data.tables) with bin number, bin logp threshold, normalized X axis info for plotting
format_2step_data <- function(data, step1_statistic, sizeBin0, alpha) {
  
  # quick function to calculate p values from chisq stats depending on method
  create_pval_info <- function(data, statistic, df=1) {
    data.table(data)[
      , step1p := pchisq(data[,statistic], df = df, lower.tail = F)
      ][
      , step2p := pchisq(data[,'chiSqGxE'],  df = 1, lower.tail = F)
      ][
      , y := -log10(step2p)
      ][
        order(step1p)
      ][
      , MapInfo := Location
      ]
  }
  
  if(step1_statistic == 'chiSqEDGE') {
    pv <- create_pval_info(data, step1_statistic, df = 2)
  } else {
    pv <- create_pval_info(data, step1_statistic, df = 1)
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
  



# file name =D
write_weighted_test_plot_title <- function(statistic) {
  gxescan_tests <- c(paste0("D|G 2-step Procedure Results (N = ", global_N, ")\noutc ~ G+", paste0(global_covs, collapse = "+"),"+", global_E), 
                     paste0("G|E 2-step Procedure Results (N = ", global_N, ")\nG ~ ", global_E, "+", paste0(global_covs, collapse = "+")),
                     paste0("EDGE 2-step Procedure Results (N = ", global_N, ")\nchiSqG + chiSqGE"))
  names(gxescan_tests) <- c("chiSqG", "chiSqGE", "chiSqEDGE")
  return(gxescan_tests[statistic])
}

write_weighted_plot_filename <- function(data, statistic) {
  return(paste0("figures/TwoStep_WeightedHypothesis_", data, "_", statistic, "_", global_E, "_", paste0(global_covs, collapse = "_"), "_N_", global_N, ".png"))
}



# create weighted hypothesis plot
# first step is creating a list object based on data.table created by the 'format_2step_data' function
# second step is actually plotting. remember that it only plots first 15 bins for legibility 
create_2step_weighted_plot <- function(data, sizeBin0, alpha, binsToPlot, statistic) {

  # m = nrow(koop_filter)
  # sizeBin0 = 5
  # alpha = 0.05
  # binsToPlot = 15
  # data = test
  
  # bin information
  m = nrow(data)
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
    t <- data[J(i)]
    t[, ref := -1*log10(min(t[,wt]))] # -log10 of bin specific alpha
    t[, x := 0.8*((t[,MapInfo]-min(t[,MapInfo])) / (max(t[,MapInfo])-min(t[,MapInfo]))) + 0.1 + i - 1]
    glist[[i]]<-t
    rm(t)
  }
  
  # trying to understand code above (just arrange BP by bins looks like)
  # (Scale mapinfo for each Bin to range between 0.1-0.9 for neatness, and add a unit increase for successive Bin)
  # x <- pv[1:5, MapInfo]
  # normalized = (x-min(x))/(max(x)-min(x))
  # normalized_scaled = 0.8 * normalized + 0.1
  # x;normalized;normalized_scaled
  
  # CREATE PLOT
  head(glist[[1]]) # for reference
  
  png(write_weighted_plot_filename(deparse(substitute(data)), statistic), width = 1280, height = 720)
  color <- rep(c("blue","olivedrab4"),100)
  plot(glist[[1]][,x], glist[[1]][,y], 
       col = "blue", 
       xlab="Bin # for step1 p-values", 
       ylab="-log10(step2 p-values)", 
       xlim=c(0,binsToPlot), 
       ylim=c(0,min.p), 
       axes=F, pch=19, cex=0.5)
  lines(glist[[1]][,x], glist[[1]][,ref],
        col = "black",lwd=2)
  
  # the rest of the points ...  =|
  # (adding to current plot..)
  for(i in 2:binsToPlot){
    points(glist[[i]][,x], glist[[i]][,y], 
           col = color[i], pch = 19, cex = 0.5)      
    lines(glist[[i]][,x], glist[[i]][,ref],
          col = "black",lwd = 2)
  }
  
  ## the last bin..
  ## it's this way because the last bin has smaller number of samples compared to bins before, thus lower bar
  ## let's only plot the first 15 bins for now, so change code a bit above. 
  # points(glist[[num]][,x], glist[[num]][,y], 
  #        col= color[num], pch = 19, cex = 0.5)
  # lines(glist[[num]][,x], rep(last.sig, nrow(glist[[num]])) ,col="black",lwd=2) # need to fix to create last horizontal line
  
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.5, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 0.8)
  axis(2, at = c(0:floor(min.p)), label = c(0:min.p), cex.axis=0.8)
  title(main = write_weighted_test_plot_title(statistic), sub = "iBin Size = 5, alpha = 0.05", cex.main = 1.5, cex.sub = 1.2)
  
  dev.off()
  
}








