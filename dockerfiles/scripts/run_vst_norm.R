#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

print("loading libraries")
suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
suppressMessages(library(reshape2))

##########################################################################################
#
# v0.0.1
# - input rna-seq count matrix is VST normalized, output is VST count matrix ''
#
##########################################################################################

# character vector of positional arguments
print("getting args")
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]       # ex. input_file <- 'expression_matrix_totalreads.tsv'

print("reading table")
counts <- read.table(input_file, row.names=1, header=T, sep='\t')
# remove rows with any NA values
counts_sans_na <- na.omit(counts)

print("making colData")
colData <- data.frame(
  row.names = colnames(counts_sans_na),
  condition = factor(c(replicate(ncol(counts_sans_na)-1, "placeholder"),"outlier"))
)

print("making DESeq matrix")
dse <- DESeqDataSetFromMatrix(countData = counts_sans_na, colData = colData, design = ~ condition)

print("running VST")
vst <- varianceStabilizingTransformation(dse)
vst_counts <- assay(vst)

print("exporting normalized count matrix")
write.table(vst_counts, file = "vst_normalized_counts_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
print("VST normalized count matrix exported to 'vst_normalized_counts_matrix.tsv'")

print("exporting normalized count table")
vst_table <- melt(vst_counts, id = 'columns')
colnames(vst_table) <- c("geneid","sample_name","vst_count")
write.table(vst_table, file = "vst_normalized_counts_table.tsv", sep = "\t", quote = FALSE, row.names=FALSE)
print("VST normalized count matrix exported to 'vst_normalized_counts_table.tsv'")