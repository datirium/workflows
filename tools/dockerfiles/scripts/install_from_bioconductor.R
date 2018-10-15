#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library("BiocInstaller")
biocLite(args[1], ask=FALSE)