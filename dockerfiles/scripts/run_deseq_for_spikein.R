#!/usr/bin/env Rscript
options(warn=-1)
options("width"=300)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(Glimma))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(hopach))
suppressMessages(library(cmapR))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))


##########################################################################################
#
#
#
# v0.0.1
#   - base script copied from `datirium/workflows/dockerfiles/scripts/run_deseq.R` (at v1.0.0)
#   - if both groups have >1 samples, runs DESeq without calculating sizeFactors by default, they are instead set manually to 1 for all samples
#   - this is intended to run RNA-Seq data that has been normalized via spike-in sequences (e.g. ERCC ExFold)
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
# DESeq/DESeq2 always compares untreated_vs_treated groups (condition-1-vs-condition-2)
# 
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Use -tn and -un to set custom names for treated and untreated conditions
#
# Use -ta and -ua to set aliases for input expression files. Should be unique
# Exports GCT and CLS files to be used by GSEA. GCT files is always with uppercase GeneId
# 
##########################################################################################


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")
RPKM_UNTREATED_ALIAS <- "RpkmCondition1"
RPKM_TREATED_ALIAS <- "RpkmCondition2"


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


export_mds_html_plot <- function(norm_counts_data, location){
    tryCatch(
        expr = {
            htmlwidgets::saveWidget(
                glimmaMDS(
                    x=assay(norm_counts_data),
                    groups=as.data.frame(SummarizedExperiment::colData(norm_counts_data)),
                    labels=rownames(SummarizedExperiment::colData(norm_counts_data))
                ),
                file=location
            )
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}


export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    tryCatch(
        expr = {
            row_metadata <- row_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            col_metadata <- col_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            print(paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export GCT data to", location, sep=" "))
        }
    )
}


export_cls <- function(categories, location){
    base::tryCatch(
        expr = {
            output_stream <- base::file(location, "w")
            on.exit(base::close(output_stream), add=TRUE)           # can't put it in 'finally' as there is no access to output_stream variable
            base::cat(
                base::paste(
                    length(categories),                             # number of datasets
                    length(base::levels(categories)),               # number of different categories
                    "1",                                            # should be always 1
                    sep="\t"
                ),
                base::paste(
                    "#",
                    base::paste(
                        base::unique(as.character(categories)),     # preserves the order, but removes duplicates
                        collapse="\t"
                    ),
                    sep="\t"
                ),
                base::paste(
                    base::paste(
                        as.character(categories),
                        collapse="\t"
                    ),
                    sep="\t"
                ),
                file=output_stream,
                sep="\n"
            )
            base::print(base::paste("Exporting CLS data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export CLS data to ", location, "with error - ", e, sep=""))
        }
    )
}


get_clustered_data <- function(expression_data, center, dist, transpose) {

    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }
    if (!is.null(center)) {
        print(paste("Centering expression data by ", center, sep=""))
        if (center == "mean"){
            expression_data = expression_data - rowMeans(expression_data)    
        } else {
            expression_data = expression_data - rowMedians(data.matrix(expression_data))    
        }
    }
    print("Creating distance matrix")
    distance_matrix <- distancematrix(expression_data, dist)
    print("Running HOPACH")
    hopach_results <- hopach(expression_data, dmat=distance_matrix)

    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }

    print("Parsing cluster labels")
    clusters = as.data.frame(hopach_results$clustering$labels)
    colnames(clusters) = "label"
    clusters = cbind(
        clusters,
        "HCL"=outer(
            clusters$label,
            10^c((nchar(trunc(clusters$label))[1]-1):0),
            function(a, b) {
                paste0("c", a %/% b %% 10)
            }
        )
    )
    clusters = clusters[, c(-1), drop=FALSE]

    return (
        list(
            order=as.vector(hopach_results$clustering$order),
            expression=expression_data,
            clusters=clusters
        )
    )
}


load_isoform_set <- function(filenames, prefixes, read_colname, rpkm_colname, rpkm_colname_alias, conditions, intersect_by, digits, batch_metadata, collected_data=NULL) {
    for (i in 1:length(filenames)) {
        isoforms <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, stringsAsFactors=FALSE)
        new_read_colname = paste(prefixes[i], conditions, sep="_")
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
    collected_data$collected_isoforms[rpkm_colname_alias] = format(rowSums(collected_data$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames), digits=digits)
    collected_data$rpkm_colnames = c(collected_data$rpkm_colnames, rpkm_colname_alias)
    collected_data$collected_isoforms <- collected_data$collected_isoforms[, !colnames(collected_data$collected_isoforms) %in% rpkm_columns]
    return( collected_data )
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
                main="Top 30 genes from VST normalized read counts",
                annotation_col=column_data,
                cluster_rows=FALSE,
                show_rownames=TRUE,
                cluster_cols=FALSE
            )
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            pheatmap(
                mat=mat_data,
                main="Top 30 genes from VST normalized read counts",
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
    parser <- ArgumentParser(description="Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2)")
    parser$add_argument(
        "-u", "--untreated",
        help="Untreated (condition 1) CSV/TSV isoforms expression files",
        type="character",required="True",
        nargs="+"
    )
    parser$add_argument(
        "-t", "--treated",
        help="Treated (condition 2) CSV/TSV isoforms expression files",
        type="character", required="True",
        nargs="+"
    )
    parser$add_argument(
        "-ua", "--ualias",
        help="Unique aliases for untreated (condition 1) expression files. Default: basenames of -u without extensions",
        type="character",
        nargs="*"
    )
    parser$add_argument(
        "-ta", "--talias",
        help="Unique aliases for treated (condition 2) expression files. Default: basenames of -t without extensions",
        type="character",
        nargs="*"
    )
    parser$add_argument(
        "-un", "--uname",
        help="Name for untreated (condition 1), use only letters and numbers",
        type="character", default="untreated"
    )
    parser$add_argument(
        "-tn", "--tname",
        help="Name for treated (condition 2), use only letters and numbers",
        type="character", default="treated"
    )
    parser$add_argument(
        "-bf", "--batchfile",
        help="Metadata file for multi-factor analysis. Headerless TSV/CSV file. First column - names from --ualias and --talias, second column - batch group name. Default: None",
        type="character"
    )
    parser$add_argument(
        "-cu", "--cutoff",
        help="Minimum threshold for rpkm filtering. Default: 5",
        type="double", default=5
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization part of the analysis output only features",
            "with adjusted P-value not bigger than this value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on normalized read counts for the",
            "exploratory visualization part of the analysis. Default: do not run",
            "clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--rowdist",
        help=paste(
            "Distance metric for HOPACH row clustering. Ignored if --cluster is not",
            "provided. Default: cosangle"
        ),
        type="character", default="cosangle",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--columndist",
        help=paste(
            "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
            "provided. Default: euclid"
        ),
        type="character", default="euclid",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for feature expression prior to running",
            "clustering by row. Ignored when --cluster is not row or both.",
            "Default: do not centered"
        ),
        action="store_true"
    )
    parser$add_argument(
        "-o", "--output",
        help="Output prefix. Default: deseq",
        type="character", default="./deseq"
    )
    parser$add_argument(
        "-d", "--digits", 
        help="Precision, number of digits to print. Default: 3",
        type="integer", default=3
    )
    parser$add_argument(
        "-p", "--threads", 
        help="Threads", 
        type="integer", default=1
    )
    parser$add_argument( # ERCC counts for the untreated
        "-uer", "--untreated_ercc_counts", # args$untreated_ercc_counts
        help="Untreated (condition 1) TSV ERCC raw count files",
        type="character",required="True",
        nargs="+"
    )
    parser$add_argument( # ERCC counts for the treated
        "-ter", "--treated_ercc_counts", # args$treated_ercc_counts
        help="Treated (condition 2) TSV ERCC raw count files",
        type="character", required="True",
        nargs="+"
    )
    parser$add_argument(
      "-uaer", "--uerccalias",
      help="Unique aliases for untreated (condition 1) ERCC files. Default: basenames of -u without extensions",
      type="character", required="True",
      nargs="*"
    )
    parser$add_argument(
      "-taer", "--terccalias",
      help="Unique aliases for treated (condition 2) ERCC files. Default: basenames of -t without extensions",
      type="character", required="True",
      nargs="*"
    )
    args <- assert_args(parser$parse_args(gsub("'|\"| ", "_", commandArgs(trailingOnly = TRUE))))
    return (args)
}

args <- get_args()

# Set threads
register(MulticoreParam(args$threads))

# Load isoforms/genes/tss files
raw_data <- load_isoform_set(args$treated, args$talias, READ_COL, RPKM_COL, RPKM_TREATED_ALIAS, args$tname, INTERSECT_BY, args$digits, args$batchfile, load_isoform_set(args$untreated, args$ualias, READ_COL, RPKM_COL, RPKM_UNTREATED_ALIAS, args$uname, INTERSECT_BY, args$digits, args$batchfile))
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

# HERE WE GENERATE THE ERCC READ TABLE, WHICH WILL BE USED BY DESEQ2 TO CALCULATE THE SIZE FACTORS ----
# A function to collect and merge ercc count tables
merge_ercc_counts <- function(ercc_file_paths, ercc_file_names) {
  if (length(ercc_file_paths) != length(ercc_file_names)) {
    stop("Number of file paths and file names must be equal.")
  }
  
  # Read each file and rename the 'count' column to the corresponding file name.
  ercc_data_list <- lapply(seq_along(ercc_file_paths), function(i) {
    df <- read.table(ercc_file_paths[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    names(df)[names(df) == "count"] <- ercc_file_names[i]
    return(df)
  })
  
  # Merge all data frames by the 'ERCC_ID' column.
  merged_df <- Reduce(function(x, y) merge(x, y, by = "ERCC_ID", all = TRUE), ercc_data_list)
  
  # Remove rows with any NA values.
  merged_df <- merged_df[complete.cases(merged_df), ]
  
  # Set the row names to the ERCC_ID column and then remove the column.
  rownames(merged_df) <- merged_df$ERCC_ID
  merged_df$ERCC_ID <- NULL
  
  return(merged_df)
}

# Run the function and generate the table ----
uer_files <- args$untreated_ercc_counts
ter_files <- args$treated_ercc_counts
all_files <- c(uer_files, ter_files)
uer_names <- paste0(args$ualias, "_", args$uname) # c("MH_Control_8h_S1", "MH_Control_8h_S2", "MH_Control_8h_S3")
ter_names <- paste0(args$talias, "_", args$tname) # c("MH_5uMNTC8_8h_S1", "MH_5uMNTC8_8h_S2", "MH_5uMNTC8_8h_S3")
all_names <- c(uer_names, ter_names)

result_table <- merge_ercc_counts(all_files, all_names) # THIS IS THE ERCC COUNT TABLE
# Filtering the ercc data table (keep rows with at least one value >=5, as the value of --cutoff)
result_table <- result_table[apply(result_table, 1, max) >= 5, ]


countData <- rbind(countData, result_table)
grep("ERCC-", rownames(countData), value = T)
print(paste("There are ERCC control genes: ", length(grep("ERCC-", rownames(countData)))))

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
    
    
    # Set the DESeq matrix
    dse <- DESeqDataSetFromMatrix(countData=countData, colData=column_data, design=design)
    # Here we let DESeq2 estimate the size factors based on ERCC counts
    dse <- estimateSizeFactors(dse, controlGenes = min(grep("ERCC-",
                        rownames(countData))):max(grep("ERCC-", rownames(countData))))
    
    # check size/normalization factors
    print("DESeq sizeFactor prior to DGE run")
    print(dse$sizeFactor)
    
    # Disable normalization by setting size factors to 1 for all samples
    # We now use the ERCC to estimate the size factors
    #sizeFactors(dse) <- 1
    
    # Removing spike-in genes from the analysis before running DESeq2
    dse <- dse[!rownames(dse) %in% rownames(result_table), ]
    
    # run DESeq
    dsq <- DESeq(dse)
    
    # Check size/normalization factors after DGE
    print("DESeq sizeFactor (dsq)")
    print(dsq$sizeFactor)

    # for norm count file. Batch correction doesn't influence it
    normCounts <- counts(dsq, normalized=TRUE)
    rownames(normCounts) <- toupper(collected_isoforms[,c("GeneId")])

    res <- results(dsq, contrast=c("conditions", args$uname, args$tname))
    export_ma_plot(res, paste(args$output, "_ma_plot", sep=""))

    # for PCA and heatmap
    vst <- varianceStabilizingTransformation(dse)

    if (!is.null(args$batchfile)){
        assay(vst) <- limma::removeBatchEffect(
            assay(vst),
            vst$batch,
            design=stats::model.matrix(stats::as.formula("~conditions"), column_data)
        )
        pca_intgroup <- c("conditions", "batch")
    } else {
        pca_intgroup <- c("conditions")
    }
    export_pca_plot(vst, paste(args$output, "_pca_plot", sep=""), pca_intgroup)
    export_mds_html_plot(                                                          # need it only whne we have at least three samples
        norm_counts_data=vst,
        location=paste(args$output, "mds_plot.html", sep="_")
    )
    vsd <- assay(vst)
    rownames(vsd) <- collected_isoforms[,c("GeneId")]
    mat <- vsd[order(rowMeans(counts(dsq, normalized=TRUE)), decreasing=TRUE)[1:30],]
    DESeqRes <- as.data.frame(res[, c(1, 2, 5, 6)])
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
    DESeqRes <- res[, c(2, 6, 7, 8)]
    colnames(DESeqRes)[3] <- "pvalue"  # in DESeq2 it's pvalue, need to use the same here
}

# Expression data heatmap of the 30 most highly expressed genes
export_heatmap(mat, column_data, paste(args$output, "_expression_heatmap", sep=""))

# Filter DESeq/DESeq2 output
DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0;
DESeqRes$pvalue[is.na(DESeqRes$pvalue)] <- 1;
DESeqRes$padj[is.na(DESeqRes$padj)] <- 1;
# DESeqRes <- format(DESeqRes, digits=args$digits)   # when applied, can't be used as numbers anymore

# Add metadata columns to the DESeq results
collected_isoforms <- data.frame(cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_count_cols], DESeqRes), check.names=F, check.rows=F)
collected_isoforms[,"'-LOG10(pval)'"] <- format(-log(as.numeric(collected_isoforms$pvalue), 10), digits=args$digits)
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

# assemble and export all count table rows
row_metadata <- collected_isoforms %>%
                dplyr::mutate_at("GeneId", toupper) %>%
                dplyr::distinct(GeneId, .keep_all=TRUE) %>%                      # to prevent from failing when input files are not grouped by GeneId
                remove_rownames() %>%
                column_to_rownames("GeneId") %>%                                 # fails if GeneId is not unique (when run with not grouped by gene data)
                dplyr::select(log2FoldChange, pvalue, padj)  %>%                 # we are interested only in these three columns
                arrange(desc(log2FoldChange))

col_metadata <- column_data %>%
                mutate_at(colnames(.), as.vector)                                # need to convert to vector, because in our metadata everything was a factor

print("Exporting all normalized read counts to GCT format")
export_gct(
    counts_mat=normCounts,
    row_metadata=row_metadata,                                   # includes features as row names
    col_metadata=col_metadata,                                   # includes samples as row names
    location=paste(args$output, "_counts_all.gct", sep="")
)
# get size of matrix before filtering
read_count_matrix_all_size <- dim(normCounts)

print(
    paste(
        "Filtering normalized read counts matrix for morpheus heatmap, to include",
        "only differentially expressed features with padj <= ", args$padj
    )
)

row_metadata <- collected_isoforms %>%
                dplyr::filter(.$padj<=args$padj) %>%
                dplyr::mutate_at("GeneId", toupper) %>%
                dplyr::distinct(GeneId, .keep_all=TRUE) %>%                      # to prevent from failing when input files are not grouped by GeneId
                remove_rownames() %>%
                column_to_rownames("GeneId") %>%                                 # fails if GeneId is not unique (when run with not grouped by gene data)
                dplyr::select(log2FoldChange, pvalue, padj)  %>%                 # we are interested only in these three columns
                arrange(desc(log2FoldChange))

col_metadata <- column_data %>%
                mutate_at(colnames(.), as.vector)                                # need to convert to vector, because in our metadata everything was a factor

print("Size of the normalized read counts matrix before filtering")
print(read_count_matrix_all_size)
normCounts <- normCounts[as.vector(rownames(row_metadata)), ]
print("Size of the normalized read counts matrix after filtering")
print(dim(normCounts))

if (!is.null(args$cluster)){
    if (args$cluster == "column" || args$cluster == "both") {
        print("Clustering filtered read counts by columns")
        clustered_data = get_clustered_data(
            expression_data=normCounts,
            center=NULL,                                              # centering doesn't influence on the samples order
            dist=args$columndist,
            transpose=TRUE
        )
        col_metadata <- cbind(col_metadata, clustered_data$clusters)  # adding cluster labels
        col_metadata <- col_metadata[clustered_data$order, ]          # reordering samples order based on the HOPACH clustering resutls
        print("Reordered samples")
        print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
        print("Clustering filtered normalized read counts by rows")
        clustered_data = get_clustered_data(
            expression_data=normCounts,
            center=if(args$center) "mean" else NULL,                  # about centering normalized data https://www.biostars.org/p/387863/
            dist=args$rowdist,
            transpose=FALSE
        )
        normCounts <- clustered_data$expression                  # can be different because of centering by rows mean
        row_metadata <- cbind(row_metadata, clustered_data$clusters)  # adding cluster labels
        row_metadata <- row_metadata[clustered_data$order, ]          # reordering features order based on the HOPACH clustering results
        print("Reordered features")
        print(head(row_metadata))
    }
}

# we do not reorder normCounts based on the clustering order
# because when exportin to GCT we use row_metadata and col_metadata
# to force the proper order of rows and columns

print("Exporting filtered normalized read counts to GCT format")
export_gct(
    counts_mat=normCounts,
    row_metadata=row_metadata,                                        # includes features as row names
    col_metadata=col_metadata,                                        # includes samples as row names
    location=paste(args$output, "_counts_filtered.gct", sep="")
)

print("Exporting CLS phenotype data")
export_cls(
    categories=col_metadata[, "conditions"],
    paste(args$output, "_phenotypes.cls", sep="")
)
