#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

install.packages("BiocInstaller", repos=args[1]);
library("BiocInstaller");
biocLite("GenomicRanges")