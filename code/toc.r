#Libraries
require(tidyverse)
require(foreach)
require(emdbook)
require(cowplot)
require(doParallel)
require(parallel)
require(RColorBrewer)
require(gridExtra)
require(gtable)
require(grid)
require(xtable)
library(stringr)

#Read in shared functions for plotting
source('shared-functions.r')
source('justsims.r')
source('abc-functions-2.r')
source('ggplot_theme.r')
#Read in graphical parameters
source('cols.r')

#Before running anything below, the set up code in the cluster_code directory must be run once
#
#Note, the paper contains millions of trials and very low tolerances for acceptance
# whereas this code is set up with a thousand trials and an extremely high tolerance. Therefore, 
# basically all the example plots are very bad. However,  the values run in the paper are commented in code
# and it is tractable to run on a local machine without parallelization. 
setwd('cluster_code')
source('run_precompute.r')
setwd('..')

#Varying FST trajectories 
#Figure F4
source('plot_fst_condensed.r')

#Figure: Vary n
#Figure F5, S2
source('fig-varyn-grob.r')

#Figure: Vary m and s (grid)
#Figure F4
source('fig-vary-m-and-s.r')

#Figure: Vary timepoints
#Figure S3, S4
source('fig-varytp-grob-2.r')

#Figure: Vary tolerance
#Figure S5
source('fig-varytol-grob.r')

#Read in for the SHIV data
source('make-shiv-data-new-3.r')
source('read-in-big-files-perc.r')
source('shiv-plotter.r')

