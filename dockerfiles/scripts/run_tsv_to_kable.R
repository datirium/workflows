#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)
#options(echo=TRUE)

suppressMessages(library(argparse))
suppressMessages(library(kableExtra))

##########################################################################################
#
# v0.0.1
# - will transform any tsv file (should have a header) into a kable markdown table
# - initially created for GSEApy summary table
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]       # ex. input_file <- 'reportsummary.tsv'
output_file <- args[2]      # ex. output_file <- 'reportsummary.md'


df <- read.table(input_file, sep="\t", header=T)

# make kable
summary_output_md <- kableExtra::kable(df, format = "html", escape = FALSE)

# Write the content to the output file
writeLines(summary_output_md, con = output_file)