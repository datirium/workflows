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
# v0.0.28
#
# - Added optional --batchcorrection parameter for DESeq2 (combatseq (provided design-formula) or limma)
# - Changed default adjusted p-value to 0.1 (in correspondance with alpha)
# - Added regulation parameter for DESeq2 (up, down, both) and appropriate altHypothesis
# - Added --lfcthreshold parameter for DESeq2 (changed default to 0.59, because it has more biological sense)
# - Changed center to min-max scaling before HOPACH clustering
# - Added blind = F parameter for DESeq2 vst and rlog (to take into account design-formula)
# - Removed RPKM filtering (not needed for DESeq2 as explained by developer while independence filtering provided)
#
#
# v1.0.0
# - Update run_deseq.R to output both all genes and filtered gene list by padj
# - copied dockerfile and run_deseq.R script from Barski lab to Datirium repo
#
#
# v0.0.27
# - Update run_deseq.R to export baseMean column
#   needed for MA-plot
#
# v0.0.26
#
# - Updated run_deseq.R with MDS plot and updated GCT export
# - Remove run_deseq_manual.R script
# - Need to install GlimmaV2 from GitHub as on official repo it's old
#
# v0.0.25
#
# - Add MDS plot and updated GCT export
#
# v0.0.24
#
# - Fix bug with pval in DESeq and pvalue in DESeq2. Now all it pvalue
#
# v0.0.23
#
# - Use RpkmCondition1 and RpkmCondition2 for RPKM columns in the output TSV file
#   We need hardcoded values for later filtering.
#
# v0.0.22
#
# - Column names for RPKM don't include spaces and brackets anymore
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
# DESeq/DESeq2 always compares untreated_vs_treated groups (condition-1-vs-condition-2)
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
RPKM_UNTREATED_ALIAS <- "RpkmCondition1"
RPKM_TREATED_ALIAS <- "RpkmCondition2"


get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
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
                        collapse=
                          "\t"
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

# Define the function to generate the markdown file
generate_md <- function(batchcorrection, batchfile, deseq_results, output_file) {
  # Initialize the markdown content
  md_content <- ""

  # Add warning message if applicable
  if (!is.null(batchcorrection) && (batchcorrection == "combatseq" || batchcorrection == "limmaremovebatcheffect") && is.null(batchfile)) {
    warning_message <- "# Warning!\n\n---\n\n**You provided a batch-correction method, but not a batch-file.**\n\nThe chosen parameter was ignored.\n\nPlease ensure that you provide a batch file when using the following batch correction methods:\n\n- **combatseq**\n- **limmaremovebatcheffect**\n\nIf you do not need batch correction, set the method to 'none'.\n\n---\n\n"
    md_content <- paste0(md_content, warning_message)
  }

  # Add DESeq results summary if provided
  if (!is.null(deseq_results)) {
    # Capture the summary output
    summary_output <- capture.output({
      summary(deseq_results)
    })

    # Parse the relevant information
    total_genes <- gsub(" with nonzero total read count", "", summary_output[2])
    total_genes <- gsub("out of ", "", total_genes)
    fdr_pvalue <- gsub("adjusted p-value < ", "", summary_output[3])
    lfc_up <- gsub(".*: ", "", summary_output[4])
    lfc_down <- gsub(".*: ", "", summary_output[5])
    outliers <- gsub(".*: ", "", summary_output[6])
    low_counts <- gsub(".*: ", "", summary_output[7])
    mean_count <- gsub("[^0-9]", "", summary_output[8])

    # Extract the LFC threshold from the summary output
    lfc_threshold <- gsub(".*LFC > ", "", summary_output[4])
    lfc_threshold <- gsub(" \\(up\\).*", "", lfc_threshold)

    # Create a data frame
    summary_df <- data.frame(
      Metric = c(
        "Total genes with non-zero read count",
        paste("LFC >", lfc_threshold, "(up)"),
        paste("LFC &lt;", lfc_threshold, "(down)"),
        "Outliers<sup>1</sup>",
        paste0("Low counts<sup>2</sup> (mean count &lt; ", mean_count, ")", collapse = "")
      ),
      Value = c(total_genes, lfc_up, lfc_down, outliers, low_counts),
      stringsAsFactors = FALSE
    )

    # Convert the data frame to markdown format using kable
    summary_output_md <- paste0(
      "# DESeq2 Results Summary\n\n---\n\n",
      "## Summary of DESeq2 Results\n\n",
      kableExtra::kable(summary_df, format = "html", escape = FALSE) %>%
        kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F),
      "\n\narguments of ?DESeq2::results():   \n<sup>1</sup> - see 'cooksCutoff',\n<sup>2</sup> - see 'independentFiltering'\n\n---\n\n"
    )

    md_content <- paste0(md_content, summary_output_md)
  }

  # Write the content to the output file
  writeLines(md_content, con = output_file)
}

get_clustered_data <- function(expression_data, dist, transpose) {
  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  expression_data <- apply(
    expression_data,
    1,
    FUN = function(x) {
      scale_min_max(x)
    }
  )

  print("Running HOPACH")
  hopach_results <- hopach(expression_data)

  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  print("Parsing cluster labels")
  clusters <- as.data.frame(hopach_results$clustering$labels)
  colnames(clusters) <- "label"
  clusters <- cbind(clusters, "HCL" = outer(clusters$label, 10^c((nchar(
    trunc(clusters$label)
  )[1] - 1):0), function(a, b) {
    paste0("c", a %/% b %% 10)
  }))
  clusters <- clusters[, c(-1), drop = FALSE]

  return(list(
    order = as.vector(hopach_results$clustering$order),
    expression = expression_data,
    clusters = clusters
  ))
}

scale_min_max <- function(x,
                          min_range = -2,
                          max_range = 2) {
  min_val <- min(x)
  max_val <- max(x)
  scaled_x <-
    (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}

load_isoform_set <- function(filenames,
                             prefixes,
                             read_colname,
                             rpkm_colname,
                             rpkm_colname_alias,
                             conditions,
                             intersect_by,
                             digits,
                             batch_metadata,
                             collected_data = NULL) {
  for (i in 1:length(filenames)) {
    isoforms <- read.table(
      filenames[i],
      sep = get_file_type(filenames[i]),
      header = TRUE,
      stringsAsFactors = FALSE
    )
    new_read_colname <- paste(prefixes[i], conditions, sep = "_")
    colnames(isoforms)[colnames(isoforms) == read_colname] <- new_read_colname
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(conditions, i, rpkm_colname,
      sep =
        " "
    )
    if (!is.null(batch_metadata)) {
      batch <- batch_metadata[prefixes[i], "batch"]
      print(
        paste(
          "Load ",
          nrow(isoforms),
          " rows from '",
          filenames[i],
          "' as '",
          new_read_colname,
          "', batch '",
          batch,
          "'",
          sep = ""
        )
      )
      column_data_frame <- data.frame(conditions, batch,
        row.names =
          c(new_read_colname)
      )
    } else {
      print(
        paste(
          "Load ",
          nrow(isoforms),
          " rows from '",
          filenames[i],
          "' as '",
          new_read_colname,
          "'",
          sep = ""
        )
      )
      column_data_frame <- data.frame(conditions, row.names = c(new_read_colname))
    }
    if (is.null(collected_data)) {
      collected_data <- list(
        collected_isoforms = isoforms,
        read_colnames = c(new_read_colname),
        column_data = column_data_frame
      )
    } else {
      collected_data$collected_isoforms <- merge(collected_data$collected_isoforms,
        isoforms,
        by = intersect_by,
        sort = FALSE
      )
      collected_data$read_colnames <- c(collected_data$read_colnames, new_read_colname)
      collected_data$column_data <- rbind(collected_data$column_data, column_data_frame)
    }
  }
  rpkm_columns <- grep(
    paste("^", conditions, " [0-9]+ ", rpkm_colname, sep = ""),
    colnames(collected_data$collected_isoforms),
    value = TRUE,
    ignore.case = TRUE
  )
  collected_data$collected_isoforms[rpkm_colname_alias] <- format(rowSums(collected_data$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames),
    digits = digits
  )
  collected_data$rpkm_colnames <- c(collected_data$rpkm_colnames, rpkm_colname_alias)
  collected_data$collected_isoforms <- collected_data$collected_isoforms[, !colnames(collected_data$collected_isoforms) %in% rpkm_columns]
  return(collected_data)
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



assert_args <- function(args) {
  print("Check input parameters")
  if (is.null(args$ualias) | is.null(args$talias)) {
    print(
      "--ualias or --talias were not set, use default values based on the expression file names"
    )
    for (i in 1:length(args$untreated)) {
      args$ualias <- append(args$ualias, head(unlist(
        strsplit(basename(args$untreated[i]), ".", fixed = TRUE)
      ), 1))
    }
    for (i in 1:length(args$treated)) {
      args$talias <- append(args$talias, head(unlist(
        strsplit(basename(args$treated[i]), ".", fixed = TRUE)
      ), 1))
    }
  } else {
    if ((length(args$ualias) != length(args$untreated)) |
      (length(args$talias) != length(args$treated))) {
      cat("\nNot correct number of inputs provided as -u, -t, -ua, -ut")
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  }

  if (length(args$treated) == 1 || length(args$untreated) == 1) {
    args$batchfile <- NULL # reset batchfile to NULL. We don't need it for DESeq even if it was provided
  }

  if (!is.null(args$batchfile)) {
    batch_metadata <- read.table(
      args$batchfile,
      sep = get_file_type(args$batchfile),
      row.names = 1,
      col.names = c("name", "batch"),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    if (all(is.element(c(args$ualias, args$talias), rownames(batch_metadata)))) {
      args$batchfile <- batch_metadata # dataframe
    } else {
      cat("\nMissing values in batch metadata file. Skipping multi-factor analysis\n")
      args$batchfile <- NULL
    }
  }

  return(args)
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
    parser$add_argument("-bf", "--batchfile",
    help = "Metadata file for multi-factor analysis. Headerless TSV/CSV file. First column - names from --ualias and --talias, second column - batch group name. Default: None", type =
      "character"
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
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
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.
      - 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.
      - 'limmaremovebatcheffect' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.
      - Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "limmaremovebatcheffect"),
    default = "none"
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold for greaterAbs and |β| < lfcThreshold for lessAbs, with p-values being two-tailed or maximum of the upper and lower tests, respectively); ",
      "'up' for upregulated genes (β > lfcThreshold in condition2 compared to condition1); ",
      "'down' for downregulated genes (β < -lfcThreshold in condition2 compared to condition1). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true"
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
        choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "abseuclid",
      "cor",
      "abscor"
    )
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
    args <- assert_args(parser$parse_args(gsub("'|\"| ", "_", commandArgs(trailingOnly = TRUE))))
    return (args)
}


args <- get_args()


# Set threads
register(MulticoreParam(args$threads))


# Load isoforms/genes/tss files
raw_data <- load_isoform_set(args$treated, args$talias, READ_COL, RPKM_COL, RPKM_TREATED_ALIAS, args$tname, INTERSECT_BY, args$digits, args$batchfile, load_isoform_set(args$untreated, args$ualias, READ_COL, RPKM_COL, RPKM_UNTREATED_ALIAS, args$uname, INTERSECT_BY, args$digits, args$batchfile))
collected_isoforms <- raw_data$collected_isoforms
read_count_cols <- raw_data$read_colnames
column_data <- raw_data$column_data
print(paste("Number of rows common for all input files ", nrow(collected_isoforms), sep=""))
print(head(collected_isoforms))
print("DESeq categories")
print(column_data)
print("DESeq count data")
countData <- collected_isoforms[read_count_cols]
print(head(countData))

# Run DESeq or DESeq2
if (length(args$treated) > 1 && length(args$untreated) > 1) {
  suppressMessages(library(DESeq2))

  if (!is.null(args$batchfile)) {
    if (args$batchcorrection == "combatseq") {
      design <- ~conditions
      print(paste(
        "Batch-correction of",
        args$batchfile,
        "using ComBat-Seq"
      ))
      counts_data <- ComBat_seq(
        as.matrix(counts_data),
        covar_mod = model.matrix(design, data = colData),
        batch = metadata[[args$batchfile]],
        # columns from counts data are alread ordered by rownames from metadata
        group = NULL
      )
    } else {
      if (args$batchcorrection == "limmaremovebatcheffect") {
        print(paste(
          "Including",
          args$batchfile,
          "as part of the design-formula"
        ))
        design <- ~ conditions + batch # We use simple +, because batch is not biologically interesting for us.
      }
    }
  } else {
    design <- ~conditions
  }

  print("Run DESeq2 analysis")

  dse <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = column_data,
    design = design
  )
  dsq <- DESeq(dse)

  # for norm count file. Batch correction doens't influence it
  normCounts <- counts(dsq, normalized = TRUE)
  rownames(normCounts) <- toupper(collected_isoforms[, c("GeneId")])

  # Determine altHypothesis and lfcThreshold based on the flag
  altHypothesis <- if (args$regulation == "up") {
    "greater"
  } else if (args$regulation == "down") {
    "less"
  } else {
    "greaterAbs"
  }

  lfcThreshold <- if (args$use_lfc_thresh) args$lfcthreshold else 0

  # Perform DESeq2 results analysis
  res <- results(
    dsq,
    contrast = c("conditions", args$uname, args$tname),
    alpha = args$fdr,
    lfcThreshold = lfcThreshold,
    independentFiltering = TRUE,
    altHypothesis = altHypothesis
  )

  generate_md(args$batchcorrection, args$batchfile, res, paste0(args$output, "_summary.md", collapse = ""))

  export_ma_plot(res, paste(args$output, "_ma_plot", sep = ""))

  # for PCA and heatmap

  if (!is.null(args$batchfile)) {
    if (args$batchcorrection == "limma") {
      vst <- DESeq2::rlog(dse, blind = FALSE)
      assay(vst) <- limma::removeBatchEffect(
        assay(vst),
        vst$batch,
        design = stats::model.matrix(stats::as.formula("~conditions"), column_data)
      )
      pca_intgroup <- c("conditions", "batch")
    } else {
      vst <- DESeq2::varianceStabilizingTransformation(dse, blind = FALSE)
      pca_intgroup <- c("conditions", "batch")
    }
  } else {
    vst <- DESeq2::varianceStabilizingTransformation(dse, blind = FALSE)
    pca_intgroup <- c("conditions")
  }

  export_pca_plot(vst, paste(args$output, "_pca_plot", sep = ""), pca_intgroup)
  export_mds_html_plot(
    # need it only whne we have at least three samples
    norm_counts_data = vst,
    location = paste(args$output, "mds_plot.html", sep = "_")
  )
  vsd <- assay(vst)
  rownames(vsd) <- collected_isoforms[, c("GeneId")]
  mat <- vsd[order(rowMeans(counts(dsq, normalized = TRUE)),
    decreasing =
      TRUE
  )[1:30], ]
  DESeqRes <- as.data.frame(res[, c(1, 2, 5, 6)])
} else {
  print("You are trying to run DESeq2 with only one replicate, which is not possible.")
}

# Expression data heatmap of the 30 most highly expressed genes
export_heatmap(
  mat,
  column_data,
  paste(args$output, "_expression_heatmap", sep = "")
)

# Filter DESeq/DESeq2 output
DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0

DESeqRes$pvalue[is.na(DESeqRes$pvalue)] <- 1

DESeqRes$padj[is.na(DESeqRes$padj)] <- 1
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
            quote=FALSE
)
print(paste("Export DESeq report to ", collected_isoforms_filename,
  sep =
    ""
))

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
    "Filtering normalized read counts matrix to include",
    "only differentially expressed features with padj <= ",
    args$fdr
  )
)

row_metadata <- collected_isoforms %>%
                dplyr::filter(.$padj<=args$fdr) %>%
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
        clustered_data <- get_clustered_data(
            expression_data = normCounts,
            dist = args$columndist,
            transpose = TRUE
        )
        col_metadata <- cbind(col_metadata, clustered_data$clusters)  # adding cluster labels
        col_metadata <- col_metadata[clustered_data$order, ]          # reordering samples order based on the HOPACH clustering resutls
        print("Reordered samples")
        print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
      print("Clustering filtered normalized read counts by rows")
      clustered_data <- get_clustered_data(
        expression_data = normCounts,
        dist = args$rowdist,
        transpose = FALSE
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