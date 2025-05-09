#!/usr/bin/env Rscript
#
# DESeq/DESeq2 differential expression analysis functions
#
# This file contains functions for performing differential expression analysis
# using DESeq2.
#
# Version: 0.1.0

#' Run DESeq2 analysis and generate results
#'
#' @param count_data Matrix or data frame of count data
#' @param col_data Data frame of sample metadata
#' @param design Design formula for DESeq2 analysis
#' @param batch_correction Batch correction method to use
#' @param batch_data Batch information (if using batch correction)
#' @param condition_names Named vector with condition names
#' @param args Command line arguments
#' @return List containing DESeq2 results and normalized counts
run_deseq2_analysis <- function(count_data, 
                                col_data, 
                                design,
                                batch_correction = "none",
                                batch_data = NULL,
                                condition_names = c(condition1 = "untreated", condition2 = "treated"),
                                args) {
  log_message("Starting DESeq2 analysis")
  
  # Check for required packages
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required but not installed")
  }
  
  # Apply batch correction with ComBat-Seq if requested
  if (batch_correction == "combatseq" && !is.null(batch_data)) {
    log_message("Applying batch correction using ComBat-Seq")
    
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("Package 'sva' is required for ComBat-Seq batch correction")
    }
    
    # Create model matrix for ComBat-Seq
    design_formula <- formula(design)
    mod <- model.matrix(design_formula, data = col_data)
    
    # Apply ComBat-Seq
    count_data <- sva::ComBat_seq(
      as.matrix(count_data),
      batch = batch_data,
      group = col_data$conditions,
      covar_mod = mod
    )
  }
  
  # Create DESeq2 dataset
  log_message("Creating DESeqDataSet")
  dse <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = design
  )
  
  # Run DESeq2
  log_message("Running DESeq2")
  dsq <- DESeq2::DESeq(dse)
  
  # Get normalized counts
  log_message("Extracting normalized counts")
  norm_counts <- DESeq2::counts(dsq, normalized = TRUE)
  
  # Determine altHypothesis based on regulation direction
  alt_hypothesis <- if (args$regulation == "up") {
    "greater"
  } else if (args$regulation == "down") {
    "less"
  } else {
    "greaterAbs"
  }
  
  # Determine LFC threshold for testing
  lfc_threshold <- if (args$use_lfc_thresh) args$lfcthreshold else 0
  
  # Get results
  log_message("Generating results with contrast")
  res <- DESeq2::results(
    dsq,
    contrast = c("conditions", condition_names["condition1"], condition_names["condition2"]),
    alpha = args$fdr,
    lfcThreshold = lfc_threshold,
    independentFiltering = TRUE,
    altHypothesis = alt_hypothesis
  )
  
  # Create VST or rlog transformation for visualization
  log_message("Creating variance-stabilized data for visualization")
  if (!is.null(batch_data) && batch_correction == "model") {
    # Apply limma's removeBatchEffect after VST
    log_message("Using limma for model-based batch correction of visualization data")
    
    if (!requireNamespace("limma", quietly = TRUE)) {
      log_warning("Package 'limma' is required for batch correction. Proceeding without batch correction.")
      vst <- DESeq2::varianceStabilizingTransformation(dse, blind = FALSE)
    } else {
      vst <- DESeq2::rlog(dse, blind = FALSE)
      batch_design <- stats::model.matrix(stats::as.formula("~conditions"), col_data)
      
      # Apply batch correction on the assay data
      assay_data <- DESeq2::assay(vst)
      corrected_data <- limma::removeBatchEffect(
        assay_data, 
        vst$batch,
        design = batch_design
      )
      
      # Replace the assay data with batch-corrected data
      DESeq2::assay(vst) <- corrected_data
    }
    
    pca_intgroup <- c("conditions", "batch")
  } else {
    # No batch correction or already applied via ComBat-Seq
    vst <- DESeq2::varianceStabilizingTransformation(dse, blind = FALSE)
    pca_intgroup <- c("conditions")
  }
  
  # Return results
  return(list(
    dse = dse,
    dsq = dsq,
    res = res,
    norm_counts = norm_counts,
    vst = vst,
    pca_intgroup = pca_intgroup
  ))
}

#' Process DESeq2 results and merge with expression data
#'
#' @param deseq_results Results from DESeq2 analysis
#' @param collected_isoforms Data frame with isoform expression data
#' @param read_colnames Vector of column names for read counts
#' @param digits Number of digits for numeric formatting
#' @return Processed and merged data frame
process_deseq_results <- function(deseq_results, collected_isoforms, read_colnames, digits) {
  log_message("Processing DESeq2 results")
  
  # Extract required columns from results
  res_df <- as.data.frame(deseq_results$res[, c("baseMean", "log2FoldChange", "pvalue", "padj")])
  
  # Handle NA values
  res_df$log2FoldChange[is.na(res_df$log2FoldChange)] <- 0
  res_df$pvalue[is.na(res_df$pvalue)] <- 1
  res_df$padj[is.na(res_df$padj)] <- 1
  
  # Merge with original data
  merged_data <- data.frame(
    cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_colnames], res_df),
    check.names = FALSE,
    check.rows = FALSE
  )
  
  # Add log-transformed p-values
  merged_data[, "'-LOG10(pval)'"] <- format(-log10(as.numeric(merged_data$pvalue)), digits = digits)
  merged_data[, "'-LOG10(padj)'"] <- format(-log10(as.numeric(merged_data$padj)), digits = digits)
  
  return(merged_data)
}

#' Generate expression heatmap of top genes
#'
#' @param vst_data Variance-stabilized data from DESeq2
#' @param count_data Raw count data
#' @param col_data Column metadata
#' @param n Number of top genes to include
#' @return Matrix of top gene expression values
get_top_expressed_genes <- function(vst_data, count_data, col_data, n = 30) {
  log_message(paste("Selecting top", n, "expressed genes for heatmap"))
  
  # Get assay data
  vsd <- DESeq2::assay(vst_data)
  
  # Order by mean expression and select top n
  top_genes <- order(rowMeans(count_data), decreasing = TRUE)[1:n]
  mat <- vsd[top_genes, ]
  
  return(mat)
} 