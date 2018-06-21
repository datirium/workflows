#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

install.packages("BiocInstaller", repos=args[1]);
install.packages("optparse", repo = "https://cloud.r-project.org/")
install.packages("sqldf", repo = "https://cloud.r-project.org/")

library("BiocInstaller");

biocLite("GenomicRanges", ask=FALSE)
biocLite("Rsamtools", ask=FALSE)
biocLite("BiocParallel", ask=FALSE)



