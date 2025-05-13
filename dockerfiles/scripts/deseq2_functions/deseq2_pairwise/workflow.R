#!/usr/bin/env Rscript
#
# Main workflow for DESeq/DESeq2 differential expression analysis
#
# This file orchestrates the entire workflow for DESeq analysis,
# from environment setup to results processing and export.
#
# Version: 0.1.0

#' Initialize the environment for DESeq analysis
#'
#' This function loads required libraries, sources dependency files,
#' and sets up error handling and logging.
initialize_environment <- function() {
  # Set options
  options(warn = -1)
  options(rlang_backtrace_on_error = "full")
  options("width" = 400)
  options(error = function() {
    message("An unexpected error occurred. Aborting script.")
    quit(save = "no", status = 1, runLast = FALSE)
  })
  
  # Set memory management options for large datasets
  options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max for global data
  options(expressions = 5000)  # Increase expression stack size
  
  # Configure garbage collection behavior
  gcinfo(FALSE)  # Disable GC messages by default
  options(gc.aggressiveness = 0)  # Default GC behavior
  
  # First load utilities which has source_with_fallback
  if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
    message("Loading utilities from Docker path: /usr/local/bin/functions/common/utilities.R")
    source("/usr/local/bin/functions/common/utilities.R")
  } else if (file.exists("functions/common/utilities.R")) {
    message("Loading utilities from relative path: functions/common/utilities.R")
    source("functions/common/utilities.R")
  } else {
    # Try one more location
    script_dir <- tryCatch({
      dirname(sys.frame(1)$ofile)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(script_dir)) {
      potential_path <- file.path(script_dir, "../common/utilities.R")
      if (file.exists(potential_path)) {
        message(paste("Loading utilities from script relative path:", potential_path))
        source(potential_path)
      } else {
        stop("Could not find common utilities.R file, which is required")
      }
    } else {
      stop("Could not find common utilities.R file, which is required")
    }
  }
  
  # Source common utility functions
  source_with_fallback("functions/common/constants.R", "/usr/local/bin/functions/common/constants.R")
  source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
  source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")
  
  # Source common visualization and export functions
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
  
  # Source DESeq-specific functions
  source_with_fallback("functions/deseq/cli_args.R", "/usr/local/bin/functions/deseq/cli_args.R")
  source_with_fallback("functions/deseq/data_processing.R", "/usr/local/bin/functions/deseq/data_processing.R")
  source_with_fallback("functions/deseq/deseq_analysis.R", "/usr/local/bin/functions/deseq/deseq_analysis.R")
  
  # Load required libraries
  load_required_libraries()

  # Configure R options
  configure_r_options()
  
  # Configure plot theme
  configure_plot_theme()
  
  # Log initialization
  log_message("Environment initialized for DESeq analysis")
}

#' Main workflow function
#'
#' @param args Command line arguments
#' @return Results of the analysis
run_workflow <- function(args) {
  log_message("Starting DESeq workflow")
  
  # Load isoforms/genes/tss files
  report_memory_usage("Before loading data")
  
  raw_data <- load_isoform_set(
    args$treated,
    args$talias,
    READ_COL,
    RPKM_COL,
    RPKM_TREATED_ALIAS,
    args$tname,
    INTERSECT_BY,
    args$digits,
    args$batchfile,
    load_isoform_set(
      args$untreated,
      args$ualias,
      READ_COL,
      RPKM_COL,
      RPKM_UNTREATED_ALIAS,
      args$uname,
      INTERSECT_BY,
      args$digits,
      args$batchfile
    )
  )
  
  # Extract data components
  collected_isoforms <- raw_data$collected_isoforms
  read_count_cols <- raw_data$read_colnames
  column_data <- raw_data$column_data
  
  log_message(paste("Number of rows common for all input files:", nrow(collected_isoforms)))
  
  # Apply RPKM filtering if requested
  if (!is.null(args$rpkm_cutoff)) {
    collected_isoforms <- filter_rpkm(collected_isoforms, args$rpkm_cutoff)
    log_message(paste("Expression data after RPKM filtering:", nrow(collected_isoforms), "rows"))
  }
  
  # Prepare count data for DESeq2
  count_data <- collected_isoforms[read_count_cols]
  log_message("Count data prepared for DESeq2 analysis")
  
  report_memory_usage("After loading data")
  
  # Run DESeq or DESeq2 based on sample count
  if (length(args$treated) > 1 && length(args$untreated) > 1) {
    log_message("Running DESeq2 analysis (multiple replicates available)")
    
    # Define design formula
    if (!is.null(args$batchfile) && args$batchcorrection == "model") {
      design <- ~conditions + batch
      log_message("Using design formula with batch effect: ~conditions + batch")
      batch_data <- args$batchfile$batch
    } else {
      design <- ~conditions
      log_message("Using standard design formula: ~conditions")
      batch_data <- NULL
    }
    
    # Run DESeq2 analysis
    deseq_results <- run_deseq2_analysis(
      count_data = count_data,
      col_data = column_data,
      design = design,
      batch_correction = args$batchcorrection,
      batch_data = batch_data,
      condition_names = c(condition1 = args$uname, condition2 = args$tname),
      args = args
    )
    
    # Generate summary markdown
    generate_md(
      args$batchcorrection, 
      args$batchfile, 
      deseq_results$res, 
      paste0(args$output, "_summary.md")
    )
    
    # Create summary plots
    create_summary_plots(
      deseq_results$dds,
      deseq_results$res,
      args$output,
      args$vst,
      args$pval,
      args$lfc
    )
    
    # Export data in various formats
    export_data(
      deseq_results$res,
      deseq_results$norm_counts,
      collected_isoforms,
      args
    )
    
    # Generate additional visualizations
    generate_visualizations(deseq_results, args)
    
    report_memory_usage("After DESeq2 analysis")
    
    # Return results
    return(deseq_results)
    
  } else {
    log_message("Running DESeq analysis with EdgeR (single replicate mode)")
    
    # Run EdgeR-based analysis for single replicate mode
    deseq_results <- run_edger_analysis(
      count_data = count_data,
      col_data = column_data,
      condition_names = c(condition1 = args$uname, condition2 = args$tname),
      args = args
    )
    
    # Export data in various formats
    export_data(
      deseq_results$res,
      deseq_results$norm_counts,
      collected_isoforms,
      args
    )
    
    # Return results
    return(deseq_results)
  }
}

#' Generate and save visualizations from DESeq2 results
#'
#' @param deseq_results Results from DESeq2 analysis
#' @param args Command line arguments
#' @return None
generate_visualizations <- function(deseq_results, args) {
  log_message("Generating visualizations")
  
  # Extract components from results
  dds <- deseq_results$dds
  res <- deseq_results$res
  
  # Create and save MA plot
  ma_plot <- function() {
    DESeq2::plotMA(res, main = "MA Plot", ylim = c(-5, 5))
  }
  save_plot(ma_plot, args$output, "ma_plot")
  
  # Create dispersion plot
  dispersion_plot <- function() {
    DESeq2::plotDispEsts(dds, main = "Dispersion Estimates")
  }
  save_plot(dispersion_plot, args$output, "dispersion_plot")
  
  # Create PCA plot if transformed data is available
  if (!is.null(deseq_results$vst_data)) {
    vst_data <- deseq_results$vst_data
    
    # Create PCA plot
    pca_plot <- DESeq2::plotPCA(vst_data, intgroup = "conditions") +
      ggplot2::ggtitle("PCA Plot") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Save PCA plot
    save_plot(pca_plot, args$output, "pca_plot")
    
    # Create interactive MDS plot
    generate_mds_plot_html(
      assay(vst_data),
      get_output_filename(args$output, "mds_plot", "html"),
      title = "MDS Plot of Samples"
    )
  }
  
  log_message("Visualizations saved successfully")
}

#' Create standard summary plots for DESeq2 analysis
#'
#' @param dds DESeq2 dds object
#' @param res DESeq2 results object
#' @param output_prefix Output file prefix
#' @param vst_transform Whether to use VST transformation
#' @param pval_threshold P-value threshold for significance
#' @param lfc_threshold Log2 fold change threshold
#' @return List of created file paths
create_summary_plots <- function(dds, res, output_prefix, vst_transform = TRUE, 
                                pval_threshold = 0.1, lfc_threshold = 1) {
  log_message("Creating summary plots for DESeq2 analysis")
  
  # Create MA plot (log fold change vs mean expression)
  log_message("Creating MA plot")
  ma_plot <- function() {
    DESeq2::plotMA(res, main = "MA Plot", ylim = c(-5, 5), alpha = pval_threshold)
  }
  
  # Save plot using consolidated function
  save_plot(ma_plot, output_prefix, "ma_plot")
  
  # Create PCA plot if data is available
  if (!is.null(dds)) {
    log_message("Creating PCA plot")
    
    # Transform data
    if (vst_transform) {
      transformed_data <- DESeq2::vst(dds, blind = FALSE)
    } else {
      transformed_data <- DESeq2::rlog(dds, blind = FALSE)
    }
    
    # Create PCA plot
    pca_plot <- DESeq2::plotPCA(transformed_data, intgroup = "conditions") +
      ggplot2::ggtitle("PCA Plot") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    # Save plot using consolidated function
    save_plot(pca_plot, output_prefix, "pca_plot")
    
    # Create expression heatmap of top differentially expressed genes
    log_message("Creating expression heatmap")
    
    # Get normalized counts
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    
    # Select top genes by adjusted p-value
    sig_genes <- res[which(res$padj < pval_threshold & abs(res$log2FoldChange) > lfc_threshold), ]
    top_genes <- head(order(sig_genes$padj), 50)
    
    if (length(top_genes) > 1) {
      # Extract count data for heatmap
      heatmap_data <- norm_counts[top_genes, ]
      
      # Scale data for better visualization
      heatmap_data_scaled <- t(scale(t(log2(heatmap_data + 1))))
      
      # Create heatmap
      heatmap_plot <- function() {
        pheatmap::pheatmap(
          heatmap_data_scaled,
          main = "Top Differentially Expressed Genes",
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          show_rownames = TRUE,
          show_colnames = TRUE,
          fontsize_row = 8,
          fontsize_col = 10
        )
      }
      
      # Save plot using consolidated function
      save_plot(heatmap_plot, output_prefix, "expression_heatmap")
    } else {
      log_message("Not enough significant genes for heatmap", "WARNING")
    }
  }
  
  # Verify expected output files with our consolidated function
  verify_outputs(output_prefix, "deseq", fail_on_missing = FALSE)
  
  log_message("Summary plots created successfully")
}

#' Export analysis results in various formats
#'
#' @param deseq_results DESeq2 results object 
#' @param norm_counts Normalized counts matrix
#' @param expression_data Full expression data 
#' @param args Command line arguments
#' @return Paths to exported files
#' @export
export_data <- function(deseq_results, norm_counts, expression_data, args) {
  log_message("Exporting data in various formats", "STEP")
  
  # Export DESeq2 results to TSV
  results_file <- export_deseq_results(deseq_results, args$output, output_name = "report")
  
  # Generate summary markdown
  summary_file <- generate_deseq_summary(
    deseq_results, 
    get_output_filename(args$output, "summary", "md"),
    parameters = list(
      "Condition 1" = args$uname,
      "Condition 2" = args$tname,
      "Batch correction" = args$batchcorrection
    )
  )
  
  # Export normalized counts in GCT format
  count_files <- export_normalized_counts(
    norm_counts, 
    args$output, 
    threshold = args$rpkm_cutoff
  )
  
  # Create CLS file for GSEA
  sample_classes <- rep(c(args$uname, args$tname), c(length(args$untreated), length(args$treated)))
  cls_file <- write_cls_file(sample_classes, get_output_filename(args$output, "phenotypes", "cls"))
  
  # Create plots with the export_visualizations function
  viz_files <- export_visualizations(
    NULL,  # We don't need the dds object here
    deseq_results, 
    args$output,
    metadata = data.frame(
      condition = sample_classes,
      row.names = colnames(norm_counts)
    )
  )
  
  # Verify all outputs were created
  verify_outputs(args$output, "deseq", fail_on_missing = FALSE)
  
  log_message("Data exported successfully", "SUCCESS")
  
  return(list(
    results = results_file,
    summary = summary_file,
    counts = count_files,
    cls = cls_file,
    viz = viz_files
  ))
}

#' Wrapper function with memory management
#'
#' @return Result of workflow execution
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message(paste("DESeq analysis started at", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Get command line arguments
  args <- get_args()
  
  # Configure parallel processing if requested
  if (!is.null(args$threads) && args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"))
    register(MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode")
  }
  
  # Run the main workflow
  results <- run_workflow(args)
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(paste("DESeq analysis completed at", format(end_time, "%Y-%m-%d %H:%M:%S")))
  log_message(paste("Total execution time:", round(as.numeric(duration), 2), "minutes"))
  
  # Clean up large objects to free memory
  rm(results)
  invisible(gc())
  
  # Final memory report
  report_memory_usage("Final")
  
  log_message("DESeq analysis completed successfully")
} 