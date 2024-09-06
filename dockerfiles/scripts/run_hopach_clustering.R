#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)
#options(echo=TRUE)

suppressMessages(library(argparse))
suppressMessages(library(hopach))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))

##########################################################################################
#
# v0.0.2
# - added step to remove any rows with 1 or more NA values
#
# v0.0.1
# - take rna-seq or na-binding data input from genelists heatmap script (run_genelists.sh) and output hopach clustering results
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]       # ex. input_file <- 'genelist_0-cluster_data.tmp'
acast_value <- args[2]      # ex. acast_value <- 'VST' ('TotalReads', 'Rpkm', 'VST', or 'zscore' for rna-seq, or 'avg_depth_sum' for na-binding data)

# read in files to data frames (counts need to be in matrix)
df <- read.table(input_file, header=T, sep='\t')


print("run_hopach_clustering.R stdout: Formatting data...")
#   acast (to matrix) opposed to dcast (as dataframe), to geneid-by-sample_name
data <- acast(df, geneid ~ sample_name, value.var=acast_value)
# replace NaN, Inf, -Inf with 0
data[!is.finite(data)] <- 0


# remove any rows with 1 or more NA values
#no_na <- data.frame(data)
#no_na <- drop_na(no_na)

#       removing not optimal, will use cap values in place of NA, Inf, -Inf
# Cap extreme values
cap_value <- 2*max(data)
# recast df to matrix (containing Inf values)
data <- acast(df, geneid ~ sample_name, value.var=acast_value)
# replace with cap value
data_capped <- data
data_capped[is.infinite(data_capped)] <- cap_value
data_capped[is.nan(data_capped)] <- 0  # Replace NaN with 0 (if any)




# cluster with hopach
print("run_hopach_clustering.R stdout: Running clusting with hopach...")
hopach_results <- hopach(data_capped, verbose=TRUE)

# save results to file for parsing
makeoutput(data_capped, hopachobj=hopach_results, file="hopach_results.out")