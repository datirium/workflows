#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)


suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))


##########################################################################################
# v0.0.3
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified based on file's extension
# *.csv - CSV
# *.tsv - TSV
# Otherwise used CSV by default
#
# The output file's rows order corresponds to the rows order of the first CSV/TSV file in
# the untreated group. Output is always saved in TSV format
#
# Output file includes only intersected rows from all input files. Intersected by
# RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
#
# DESeq/DESeq2 always compares untreated_vs_treated groups
##########################################################################################


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
U_PREFIX <- "u"
T_PREFIX <- "t"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}

load_isoform_set <- function(filenames, prefix, read_colname, rpkm_colname, gname_suffix, intersect_by, collected_isoforms=NULL) {
    for (i in 1:length(filenames)) {
        isoforms <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, stringsAsFactors=FALSE)
        print(paste("Load ", nrow(isoforms), " rows from ", filenames[i], sep=""))
        colnames(isoforms)[colnames(isoforms) == read_colname] <- paste(prefix, i, read_colname, sep="")
        colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(prefix, i, rpkm_colname, sep="")
        if (is.null(collected_isoforms)){
            collected_isoforms <- isoforms
        } else {
            collected_isoforms <- merge(collected_isoforms, isoforms, by=intersect_by, sort = FALSE)
        }
    }
    rpkm_columns = grep(paste("^", prefix, "[0-9]", rpkm_colname, sep=""), colnames(collected_isoforms), value = TRUE, ignore.case = TRUE)
    collected_isoforms[paste(rpkm_colname, gname_suffix, sep="")] = rowSums(collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames)
    collected_isoforms <- collected_isoforms[, !colnames(collected_isoforms) %in% rpkm_columns]
    return (collected_isoforms)
}


# Parser
parser <- ArgumentParser(description='Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups')
parser$add_argument("-u", "--untreated", help='Untreated CSV/TSV isoforms expression files',    type="character", required="True", nargs='+')
parser$add_argument("-t", "--treated",   help='Treated CSV/TSV isoforms expression files',      type="character", required="True", nargs='+')
parser$add_argument("-un","--uname",     help='Suffix for untreated RPKM column name',      type="character", default="_u")
parser$add_argument("-tn","--tname",     help='Suffix for treated RPKM column name',        type="character", default="_t")
parser$add_argument("-o", "--output",    help='Output TSV filename',    type="character", default="./deseq_results.tsv")
parser$add_argument("-p", "--threads",   help='Threads',            type="integer",   default=1)
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


# Set threads
register(MulticoreParam(args$threads))


# Define conditions for DESeq
conditions <- c(rep("untreated", length(args$untreated)), rep("treated", length(args$treated)))
column_data <- data.frame(conditions, row.names=c(paste(U_PREFIX, 1:length(args$untreated), READ_COL, sep=""), paste(T_PREFIX, 1:length(args$treated), READ_COL, sep="")))


# Load isoforms/genes/tss files. Select columns with read data
collected_isoforms <- load_isoform_set(args$treated, T_PREFIX, READ_COL, RPKM_COL, args$tname, INTERSECT_BY, load_isoform_set(args$untreated, U_PREFIX, READ_COL, RPKM_COL, args$uname, INTERSECT_BY))
read_count_cols = grep(paste("^", U_PREFIX, "[0-9]+", READ_COL, "|^", T_PREFIX, "[0-9]+", READ_COL, sep=""), colnames(collected_isoforms), value = TRUE, ignore.case = TRUE)
print(paste("Number of rows common for all input files ", nrow(collected_isoforms), sep=""))


# Run DESeq or DESeq2
if (length(args$treated) > 1 && length(args$untreated) > 1){
    print("Run DESeq2")
    suppressMessages(library(DESeq2))
    dse <- DESeqDataSetFromMatrix(countData=collected_isoforms[read_count_cols], colData=column_data, design=~conditions)
    dsq <- DESeq(dse)
    DESeqRes <- as.data.frame(results(dsq, contrast=c("conditions", "treated", "untreated"))[,c(2,5,6)])
} else {
    print("Run DESeq")
    suppressMessages(library(DESeq))
    cds <- newCountDataSet(collected_isoforms[read_count_cols], conditions)
    cdsF <- estimateSizeFactors(cds)
    cdsD <- estimateDispersions(cdsF, method="blind", sharingMode="fit-only", fitType="local")
    DESeqRes <- nbinomTest(cdsD, "untreated", "treated")
    isinfl <- is.infinite(DESeqRes$log2FoldChange)
    DESeqRes$log2FoldChange[isinfl] <- log2((DESeqRes$baseMeanB[isinfl]+0.1)/(DESeqRes$baseMeanA[isinfl]+0.1))
    DESeqRes <- DESeqRes[,c(6,7,8)]
}


# Filter DESeq/DESeq2 output
DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0;
DESeqRes[is.na(DESeqRes)] <- 1;


# Export results to file
collected_isoforms <- data.frame(cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_count_cols], DESeqRes), check.names=F, check.rows=F)
write.table(collected_isoforms,
            file=args$output,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
print(paste("Export results to ", args$output, sep=""))
