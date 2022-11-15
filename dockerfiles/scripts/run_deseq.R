#!/usr/bin/env Rscript
options(warn=-1)
options("width"=300)


suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))


##########################################################################################
#
# v0.0.21
#
# - Add ggrepel for proper label positioning
#
# v0.0.20
#
# - Add --batchfile parameter to run_deseq.R to compensate batch effect
#
# v0.0.19
#
# - Add --batchfile to compensate batch effect. Works only for DESeq2
#
# v0.0.18
#
# - Add --digits parameter to set a precision in output table
#
# v0.0.17
#
# - Update labels in cls, replaces n/a with na in gct files
#
# v0.0.16
#
# - Add max(rpkm) cutoff filtering
#
# v0.0.15
#
# - fix bug with " and ' in arguments. Replace all with ""
#
# v0.0.14
#
# - add PCA plot
#
# v0.0.13
#
# - Fix bug in phenotype.cls column order
# - Fix bug in logFC sign for DESeq2
#
# v0.0.8
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
# 
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Use -un and -tn to set custom names for treated and untreated conditions
#
# Use -ua and -ta to set aliases for input expression files. Should be unique
# Exports GCT and CLS files to be used by GSEA. GCT files is always with uppercase GeneId
# 
##########################################################################################


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


load_isoform_set <- function(filenames, prefixes, read_colname, rpkm_colname, conditions, intersect_by, digits, batch_metadata, collected_data=NULL) {
    for (i in 1:length(filenames)) {
        isoforms <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, stringsAsFactors=FALSE)
        new_read_colname = paste(prefixes[i], " [", conditions, "]", sep="")
        colnames(isoforms)[colnames(isoforms) == read_colname] <- new_read_colname
        colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(conditions, i, rpkm_colname, sep=" ")
        if (!is.null(batch_metadata)){
            batch <- batch_metadata[prefixes[i], "batch"]
            print(paste("Load ", nrow(isoforms), " rows from '", filenames[i], "' as '", new_read_colname, "', batch '", batch, "'", sep=""))
            column_data_frame <- data.frame(conditions, batch, row.names=c(new_read_colname))
        } else {
            print(paste("Load ", nrow(isoforms), " rows from '", filenames[i], "' as '", new_read_colname, "'", sep=""))
            column_data_frame <- data.frame(conditions, row.names=c(new_read_colname))
        }
        if (is.null(collected_data)){
            collected_data = list(collected_isoforms=isoforms, read_colnames=c(new_read_colname), column_data=column_data_frame)
        } else {
            collected_data$collected_isoforms <- merge(collected_data$collected_isoforms, isoforms, by=intersect_by, sort = FALSE)
            collected_data$read_colnames = c(collected_data$read_colnames, new_read_colname)
            collected_data$column_data <- rbind(collected_data$column_data, column_data_frame)
        }
    }
    rpkm_columns = grep(paste("^", conditions, " [0-9]+ ", rpkm_colname, sep=""), colnames(collected_data$collected_isoforms), value = TRUE, ignore.case = TRUE)
    new_rpkm_colname = paste(rpkm_colname, " [", conditions, "]", sep="")
    collected_data$collected_isoforms[new_rpkm_colname] = format(rowSums(collected_data$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames), digits=digits)
    collected_data$rpkm_colnames = c(collected_data$rpkm_colnames, new_rpkm_colname)
    collected_data$collected_isoforms <- collected_data$collected_isoforms[, !colnames(collected_data$collected_isoforms) %in% rpkm_columns]
    return( collected_data )
}


write.gct <- function(gct, filename) {
	rows <- dim(gct$data)[1]
	columns <- dim(gct$data)[2]
	rowDescriptions <- gct$rowDescriptions
	m <- cbind(row.names(gct$data), rowDescriptions, gct$data)
	f <- file(filename, "w")
	on.exit(close(f))
	cat("#1.2", "\n", file=f, append=TRUE, sep="")
	cat(rows, "\t", columns, "\n", file=f, append=TRUE, sep="")
	cat("Name", "\t", file=f, append=TRUE, sep="")
	cat("Description", file=f, append=TRUE, sep="")
	names <- colnames(gct$data)
	for(j in 1:length(names)) {
		cat("\t", names[j], file=f, append=TRUE, sep="")
	}
	cat("\n", file=f, append=TRUE, sep="")
	write.table(m, file=f, append=TRUE, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
}


write.cls <- function(factor, filename) {
	file <- file(filename, "w")
	on.exit(close(file))
 	codes <- as.character(factor)
	cat(file=file, length(codes), length(levels(factor)), "1\n")
	levels <- levels(factor)
	cat(file=file, "# ")
	num.levels <- length(levels)
    if(num.levels-1 != 0) {
	    for(i in 1:(num.levels-1)) {
		    cat(file=file, levels[i])
		    cat(file=file, " ")
	    }
	}
	cat(file=file, levels[num.levels])
	cat(file=file, "\n")
	num.samples <- length(codes)
	if(num.samples-1 != 0) {
	    for(i in 1:(num.samples-1)) {
		    cat(file=file, codes[i])
		    cat(file=file, " ")
	    }
	}
	cat(file=file, codes[num.samples])
}


export_ma_plot <- function(data, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            plotMA(data)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            plotMA(data)
            dev.off()

            cat(paste("\nExport MA-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export MA-plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_plot <- function(data, rootname, intgroup, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            pca_data <- plotPCA(data, intgroup=intgroup, returnData=TRUE)
            percentVar <- round(100 * attr(pca_data, "percentVar"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(pca_data, aes(PC1, PC2, color=group)) +
                geom_point(size=5, shape=19) +
                xlab(paste0("PC1: ",percentVar[1], "% variance")) +
                ylab(paste0("PC2: ",percentVar[2], "% variance")) + 
                geom_label_repel(aes(label=name), point.padding=0.5, box.padding=0.5, check_overlap = TRUE, show.legend = FALSE)
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ggplot(pca_data, aes(PC1, PC2, color=group)) +
                geom_point(size=5, shape=19) +
                xlab(paste0("PC1: ",percentVar[1], "% variance")) +
                ylab(paste0("PC2: ",percentVar[2], "% variance")) + 
                geom_label_repel(aes(label=name), point.padding=0.5, box.padding=0.5, check_overlap = TRUE, show.legend = FALSE)
            )
            dev.off()

            cat(paste("\nExport PCA-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA-plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_heatmap <- function(mat_data, column_data, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            pheatmap(
                mat=mat_data,
                annotation_col=column_data,
                cluster_rows=FALSE,
                show_rownames=TRUE,
                cluster_cols=FALSE
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            pheatmap(
                mat=mat_data,
                annotation_col=column_data,
                cluster_rows=FALSE,
                show_rownames=TRUE,
                cluster_cols=FALSE
            )
            dev.off()

            cat(paste("\nExport expression heatmap to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export expression heatmap to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


apply_cutoff <- function(data, cutoff, colnames){
    print(paste("Apply rpkm cutoff ", cutoff, sep=""))
    data = data[data[,colnames[1]] >= cutoff | data[,colnames[2]] >= cutoff,]
    return( data )
}


assert_args <- function(args){
    print("Check input parameters")
    if (is.null(args$ualias) | is.null(args$talias)){
        print("--ualias or --talias were not set, use default values based on the expression file names")
        for (i in 1:length(args$untreated)) {
            args$ualias = append(args$ualias, head(unlist(strsplit(basename(args$untreated[i]), ".", fixed = TRUE)), 1))
        }
        for (i in 1:length(args$treated)) {
            args$talias = append(args$talias, head(unlist(strsplit(basename(args$treated[i]), ".", fixed = TRUE)), 1))
        }
    } else {
        if ( (length(args$ualias) != length(args$untreated)) | (length(args$talias) != length(args$treated)) ){
            cat("\nNot correct number of inputs provided as -u, -t, -ua, -ut")
            quit(save = "no", status = 1, runLast = FALSE)
        }
    }

    if (length(args$treated) == 1 || length(args$untreated) == 1){
        args$batchfile <- NULL  # reset batchfile to NULL. We don't need it for DESeq even if it was provided
    }

    if (!is.null(args$batchfile)){
        batch_metadata <- read.table(
            args$batchfile,
            sep=get_file_type(args$batchfile),
            row.names=1,
            col.names=c("name", "batch"),
            header=FALSE,
            stringsAsFactors=FALSE
        )
        rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
        if (all(is.element(c(args$ualias, args$talias), rownames(batch_metadata)))){
            args$batchfile <- batch_metadata  # dataframe
        } else {
            cat("\nMissing values in batch metadata file. Skipping multi-factor analysis\n")
            args$batchfile = NULL
        }
    } 

    return (args)
}


get_args <- function(){
    parser <- ArgumentParser(description='Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups')
    parser$add_argument("-u",  "--untreated", help='Untreated CSV/TSV isoforms expression files',    type="character", required="True", nargs='+')
    parser$add_argument("-t",  "--treated",   help='Treated CSV/TSV isoforms expression files',      type="character", required="True", nargs='+')
    parser$add_argument("-ua", "--ualias",    help='Unique aliases for untreated expression files. Default: basenames of -u without extensions', type="character", nargs='*')
    parser$add_argument("-ta", "--talias",    help='Unique aliases for treated expression files. Default: basenames of -t without extensions',   type="character", nargs='*')
    parser$add_argument("-un", "--uname",     help='Name for untreated condition, use only letters and numbers', type="character", default="untreated")
    parser$add_argument("-tn", "--tname",     help='Name for treated condition, use only letters and numbers',   type="character", default="treated")
    parser$add_argument("-bf", "--batchfile", help='Metadata file for multi-factor analysis. Headerless TSV/CSV file. First column - names from --ualias and --talias, second column - batch group name. Default: None', type="character")
    parser$add_argument("-cu", "--cutoff",    help='Minimum threshold for rpkm filtering. Default: 5', type="double", default=5)
    parser$add_argument("-o",  "--output",    help='Output prefix. Default: deseq',    type="character", default="./deseq")
    parser$add_argument("-d",  "--digits",    help='Precision, number of digits to print. Default: 3', type="integer", default=3)
    parser$add_argument("-p",  "--threads",   help='Threads',            type="integer",   default=1)
    args <- assert_args(parser$parse_args(gsub("'|\"| ", "_", commandArgs(trailingOnly = TRUE))))
    return (args)
}


args <- get_args()


# Set threads
register(MulticoreParam(args$threads))


# Load isoforms/genes/tss files
raw_data <- load_isoform_set(args$treated, args$talias, READ_COL, RPKM_COL, args$tname, INTERSECT_BY, args$digits, args$batchfile, load_isoform_set(args$untreated, args$ualias, READ_COL, RPKM_COL, args$uname, INTERSECT_BY, args$digits, args$batchfile))
collected_isoforms <- apply_cutoff(raw_data$collected_isoforms, args$cutoff, raw_data$rpkm_colnames)
read_count_cols = raw_data$read_colnames
column_data = raw_data$column_data
print(paste("Number of rows common for all input files ", nrow(collected_isoforms), sep=""))
print(head(collected_isoforms))
print("DESeq categories")
print(column_data)
print("DESeq count data")
countData = collected_isoforms[read_count_cols]
print(head(countData))

# Run DESeq or DESeq2
if (length(args$treated) > 1 && length(args$untreated) > 1){
    
    suppressMessages(library(DESeq2))

    if (!is.null(args$batchfile)){
        print("Run DESeq2 multi-factor analysis")
        design=~batch+conditions  # We use simple +, because batch is not biologically interesting for us.
    } else {
        print("Run DESeq2 analysis")
        design=~conditions
    }
    
    dse <- DESeqDataSetFromMatrix(countData=countData, colData=column_data, design=design)
    dsq <- DESeq(dse)

    # for norm count file. Batch correction doens't influence it
    normCounts <- counts(dsq, normalized=TRUE)
    rownames(normCounts) <- toupper(collected_isoforms[,c("GeneId")])

    res <- results(dsq, contrast=c("conditions", args$uname, args$tname))
    export_ma_plot(res, paste(args$output, "_ma_plot", sep=""))

    # for PCA and heatmap
    vst <- varianceStabilizingTransformation(dse)

    if (!is.null(args$batchfile)){
        assay(vst) <- limma::removeBatchEffect(assay(vst), vst$batch)
        pca_intgroup <- c("conditions", "batch")
    } else {
        pca_intgroup <- c("conditions")
    }
    export_pca_plot(vst, paste(args$output, "_pca_plot", sep=""), pca_intgroup)

    vsd <- assay(vst)
    rownames(vsd) <- collected_isoforms[,c("GeneId")]
    mat <- vsd[order(rowMeans(counts(dsq, normalized=TRUE)), decreasing=TRUE)[1:30],]

    DESeqRes <- as.data.frame(res[,c(2,5,6)])
} else {
    print("Run DESeq analysis")
    suppressMessages(library(DESeq))
    cds <- newCountDataSet(countData, column_data[,"conditions"])
    cdsF <- estimateSizeFactors(cds)
    cdsD <- estimateDispersions(cdsF, method="blind", sharingMode="fit-only", fitType="local")
    normCounts <- counts(cdsD, normalized=TRUE)
    rownames(normCounts) <- toupper(collected_isoforms[,c("GeneId")])
    res <- nbinomTest(cdsD, args$uname, args$tname)
    infLFC <- is.infinite(res$log2FoldChange)
    res$log2FoldChange[infLFC] <- log2((res$baseMeanB[infLFC]+0.1)/(res$baseMeanA[infLFC]+0.1))

    export_ma_plot(res, paste(args$output, "_ma_plot", sep=""))

    vsd <- exprs(varianceStabilizingTransformation(cdsD))
    rownames(vsd) <- collected_isoforms[,c("GeneId")]
    mat <- vsd[order(rowMeans(counts(cdsD, normalized=TRUE)), decreasing=TRUE)[1:30],]

    DESeqRes <- res[,c(6,7,8)]
}


# Normalized counts table for GCT export
normCountsGct <- list(rowDescriptions=c(rep("na", times=length(row.names(normCounts)))), data=as.matrix(normCounts))


# Create phenotype table for CLS export
phenotype_labels <- gsub("\\s|\\t", "_", column_data[colnames(normCounts), "conditions"])
phenotype_data <- as.factor(phenotype_labels)
phenotype_data <- factor(phenotype_data, levels=unique(phenotype_labels))

# Expression data heatmap of the 30 most highly expressed genes
export_heatmap(mat, column_data, paste(args$output, "_expression_heatmap", sep=""))

# Filter DESeq/DESeq2 output
DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0;
DESeqRes[is.na(DESeqRes)] <- 1;
DESeqRes <- format(DESeqRes, digits=args$digits)


# Add metadata columns to the DESeq results
collected_isoforms <- data.frame(cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_count_cols], DESeqRes), check.names=F, check.rows=F)
collected_isoforms[,"'-LOG10(pval)'"] <- format(-log(as.numeric(collected_isoforms$pval), 10), digits=args$digits)
collected_isoforms[,"'-LOG10(padj)'"] <- format(-log(as.numeric(collected_isoforms$padj), 10), digits=args$digits)


# Export DESeq results to file
collected_isoforms_filename <- paste(args$output, "_report.tsv", sep="")
write.table(collected_isoforms,
            file=collected_isoforms_filename,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
print(paste("Export DESeq report to ", collected_isoforms_filename, sep=""))


# Export DESeq normalized counts to GSEA compatible file
gct_filename <- paste(args$output, "_counts.gct", sep="")
write.gct(normCountsGct, file=gct_filename)
print(paste("Export normalized counts to ", gct_filename, sep=""))


# Export phenotype data to GSEA compatible file
cls_filename <- paste(args$output, "_phenotypes.cls", sep="")
write.cls(phenotype_data, file=cls_filename)
print(paste("Export phenotype data to ", cls_filename, sep=""))