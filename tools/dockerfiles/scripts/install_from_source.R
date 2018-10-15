#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
install.packages(args[1], repos=NULL, type="source")