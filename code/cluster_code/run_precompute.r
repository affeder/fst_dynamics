#Note: this precomputation was all run on a cluster
#I've created a version that can be run locally, but it is not very practical to run locally.
#That this was written originally for the cluster is why files are often indexed by job ids 
library(foreach)
library(doParallel)
library(emdbook)
library(dplyr)
library(tidyverse)
library(untb)
library(abc)
library(emdbook)


#Step 1: Import necessary functions
source('ABC_proc.r')

#Step 2: Set up directories:
system("mkdir ../../raw")
system("mkdir ../../raw/prior")
system("mkdir ../../raw/stats")
system("mkdir ../../raw/stats/n30")
system("mkdir ../../raw/stats/n100")
system("mkdir ../../raw/stats/n500")
system("mkdir ../../raw/svals")
system("mkdir ../../raw/svals/n30")
system("mkdir ../../raw/svals/n100")
system("mkdir ../../raw/svals/n500")
system("mkdir ../../out")
system("mkdir ../../graphs")
system("mkdir ../../output")

#Step 3: Generate priors
source('gen-abc-prior-trials.r')

#Step 4: Process priors into summary statistics
source('parallel-process-stats.r')
source('computeTP_s_forSherlock.r')

#Step 5: Compress files across different runs
system("cat ../../raw/svals/n30/sout_*-gen2.txt > ../../raw/svals/svals_n30_gen2.txt")
system("cat ../../raw/svals/n100/sout_*-gen2.txt > ../../raw/svals/svals_n100_gen2.txt")
system("cat ../../raw/svals/n500/sout_*-gen2.txt > ../../raw/svals/svals_n500_gen2.txt")

system("cat ../../raw/svals/n30/sout_*-gen5.txt > ../../raw/svals/svals_n30_gen5.txt")
system("cat ../../raw/svals/n100/sout_*-gen5.txt > ../../raw/svals/svals_n100_gen5.txt")
system("cat ../../raw/svals/n500/sout_*-gen5.txt > ../../raw/svals/svals_n500_gen5.txt")


system("cat ../../raw/stats/n30/stats_output_*.txt > ../../raw/stats/n30_prior.txt")
system("cat ../../raw/stats/n100/stats_output_*.txt > ../../raw/stats/n100_prior.txt")
system("cat ../../raw/stats/n500/stats_output_*.txt > ../../raw/stats/n500_prior.txt")

#Step 5: Run ABC process to create various figures
########Figure 4, S3
source('nrun.r')
system("cat ../../out/sumStats-varyn-n*.txt > ../../out/sumStats-varyn.txt")
########Figure 5
source('svm.r')
########Figure S4, S5
source('run_tp.r')
########Figure S6
source('tolrun.r')
########Figure 4, S3
source('nrun-asym.r')
