#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

suppressMessages(library(argparse))
suppressMessages(library(hopach))
suppressMessages(library(reshape2))

##########################################################################################
#
# v0.0.1
# - take rna-seq or na-binding data input from genelists heatmap script (run_genelists.sh) and output hopach clustering results
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]       # ex. input_file <- 'genelist_0-cluster_data.tmp'
acast_value <- args[2]      # ex. acast_value <- 'VST' ('Rpkm' or 'VST' for rna-seq, or 'avg_depth_sum' for na-binding data)

# read in files to data frames (counts need to be in matrix)
df <- read.table(input_file, header=T, sep='\t')

#   acast (to matrix) opposed to dcast (as dataframe), to geneid-by-sample_name
data <- acast(df, geneid ~ sample_name, value.var=acast_value)

# cluster with hopach
hopach_results <- hopach(data)

# save results to file for parsing
makeoutput(data, hopachobj=hopach_results, file="hopach_results.out")