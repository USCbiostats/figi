library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
source("~/Dropbox/FIGI/Code/Functions/GxEScan_PostHoc_Analyses.R")

#-----------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------#

# Manhattan Plot
create_manhattanplot(gxe, 'chiSqG', df = 1)

#-----------------------------------------------------------------------------#
# GxE results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSqGxE', df = 1)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSqGxE', df = 1)

#-----------------------------------------------------------------------------#
# 2DF results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSq2df', df = 2)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSq2df', df = 2)


#-----------------------------------------------------------------------------#
# 3DF results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSq3df', df = 3)
# Manhattan Plot
create_manhattanplot(gxe, 'chiSq3df', df = 3)


#-----------------------------------------------------------------------------#
# GE, Case, Control ----
#-----------------------------------------------------------------------------#

# GE QQ Plot
create_qqplot(gxe, 'chiSqGE', df = 1)

# Control-Only
create_qqplot(gxe, 'chiSqControl', df = 1)

# Case-Only QQ + manhattan
create_qqplot(gxe, 'chiSqCase', df = 1)
create_manhattanplot(gxe, 'chiSqCase', df = 1)