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

print("exporting VST count matrix")
write.table(vst_counts, file = "vst_counts_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
print("VST matrix exported to 'vst_counts_matrix.tsv'")

print("exporting VST count table")
table <- melt(vst_counts, id = 'columns')
colnames(table) <- c("geneid","sample_name","vst_count")
write.table(table, file = "vst_counts_table.tsv", sep = "\t", quote = FALSE, row.names=FALSE)
print("VST table exported to 'vst_counts_table.tsv'")



#   Z-SCORE calculations for each gene across all samples
# Calculate row means and standard deviations
print("running z-score calcs")
gene_means <- rowMeans(vst_counts)
gene_sds <- apply(vst_counts, 1, sd)
# Calculate z-scores
z_scores <- t((t(vst_counts) - gene_means) / gene_sds)

# replace Inf or NaN with 0 (not appropriate to do, so will leave this step commented out)
#z_scores[!is.finite(z_scores)] <- 0
# but we ccould simply remove them, we must note that in the output file for z-score heatmap
# leaving commented out for now, add this quote to workflow output if uncommented "(-Inf and Inf values are removed from the heatmap, these are due to near-zero or huge [relative to the mean] standard deviations, respectively)"
#z_scores <- z_scores[!is.infinite(rowSums(z_scores)),]


print("exporting VST z-score matrix")
write.table(z_scores, file = "vst_z-score_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
print("VST z-score matrix exported to 'vst_zscore_matrix.tsv'")

print("exporting VST z-score table")
table <- melt(z_scores, id = 'columns')
colnames(table) <- c("geneid","sample_name","vst_count")
write.table(table, file = "vst_zscore_table.tsv", sep = "\t", quote = FALSE, row.names=FALSE)
print("VST z-score table exported to 'vst_zscore_table.tsv'")
