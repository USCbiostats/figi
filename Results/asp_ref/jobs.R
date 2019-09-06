library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(animation)


#-----------------------------------------------------------------------------#
# Marginal G Results ----
#-----------------------------------------------------------------------------#
# QQ Plot
create_qqplot(gxe, 'chiSqG', df = 1)
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

#-----------------------------------------------------------------------------#
# D|G 2-step Kooperberg ----
#-----------------------------------------------------------------------------#
gxe_twostep <- format_2step_data(data = gxe, 'chiSqG', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqG')

#-----------------------------------------------------------------------------#
# G|E 2-step Murcray ----
#-----------------------------------------------------------------------------#
gxe_twostep <- format_2step_data(data = gxe, 'chiSqGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqGE')

#-----------------------------------------------------------------------------#
# EDGE 2-step Gauderman ----
#-----------------------------------------------------------------------------#
gxe_twostep <- format_2step_data(data = gxe, 'chiSqEDGE', 5, 0.05)
create_2step_weighted_plot(gxe_twostep, sizeBin0 = 5, alpha = 0.05, binsToPlot = 10, statistic = 'chiSqEDGE')