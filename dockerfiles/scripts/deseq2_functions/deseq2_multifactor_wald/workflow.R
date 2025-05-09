#!/usr/bin/env Rscript
#
# Workflow functions for DESeq2 LRT Step 2
#

#' Main workflow function with memory management
#' 
#' Wraps the run_workflow function with memory tracking
#' 
#' @return Results from the analysis workflow
#' @export
main_with_memory_management <- function() {
  # Report initial memory usage
  report_memory_usage("Initial")
  
  # Get command-line arguments
  args <- get_args()
  
  # Run the main workflow
  results <- tryCatch({
    run_workflow(args)
  }, error = function(e) {
    log_message(paste("Error in workflow execution:", e$message), "ERROR")
    # Force garbage collection
    gc(verbose = FALSE)
    # Re-raise the error
    stop(e)
  })
  
  # Final memory usage report
  report_memory_usage("Final")
  
  return(results)
}

#' Main workflow function for DESeq2 LRT Step 2
#'
#' @param args Command-line arguments
#' @return Final results list containing normalized counts, expression data, and contrast results
#' @export

initialize_environment <- function() {
  # First, make sure we have the utilities module
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
        stop("Could not find utilities.R file")
      }
    } else {
      stop("Could not find utilities.R file")
    }
  }
  
  # Now we have access to source_with_fallback and other utilities
  # Source common functions
  source_with_fallback("functions/common/constants.R", "/usr/local/bin/functions/common/constants.R")
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
  source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
  source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")

  # Source DESeq2 LRT Step 2 specific functions
  source_with_fallback("functions/deseq2_lrt_step_2/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_2/cli_args.R")
  source_with_fallback("functions/deseq2_lrt_step_2/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_2/data_processing.R")
  source_with_fallback("functions/deseq2_lrt_step_2/contrast_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_2/contrast_analysis.R")
  
  # Load required libraries
  load_required_libraries()

  # Configure R options
  configure_r_options()
  
  # Configure plot theme
  configure_plot_theme()
  
  log_message("Environment initialized for DESeq2 LRT Step 2 analysis")
}


run_workflow <- function(args) {
  log_message("Starting DESeq2 LRT Step 2 workflow")
  
  # Load input data from Step 1
  log_message(paste("Loading DESeq2 object data from", args$dsq_obj_data))
  step1_data <- readRDS(args$dsq_obj_data)
  
  # Load contrast data
  log_message(paste("Loading contrast data from", args$contrast_df))
  contrast_df <- read.delim(args$contrast_df, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Get contrast indices to process
  contrast_indices <- as.numeric(unlist(strsplit(args$contrast_indices, ",")))
  log_message(paste("Processing contrast indices:", paste(contrast_indices, collapse = ", ")))
  
  # Check if contrast indices are valid
  if (any(contrast_indices > nrow(contrast_df))) {
    stop("Invalid contrast index. Maximum index is ", nrow(contrast_df))
  }
  
  # Subset contrast data to process only the specified contrasts
  contrast_df <- contrast_df[contrast_indices, , drop = FALSE]
  
  # Initialize variables
  dds <- step1_data$dds
  expDataDf <- step1_data$expDataDf
  
  # Check for required columns in the DESeq dataset
  if (is.null(dds)) {
    stop("DESeq2 dataset is missing from step 1 data")
  }
  
  # Extract metadata from DESeq2 object
  col_data <- as.data.frame(colData(dds))
  design_formula <- design(dds)
  
  # Verify metadata and consistency between steps
  col_data <- validate_metadata(col_data, args$batchcorrection, design_formula)
  verify_step_consistency(col_data, design_formula, args)
  
  # Identify all factors in the DESeq dataset for reference
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  log_message(paste("Available factors in DESeq2 dataset:", paste(factor_names, collapse = ", ")))
  
  # Rebuild dds with reference levels if needed
  reference_levels <- list()
  if (exists("factors", where = step1_data)) {
    for (factor_name in step1_data$factors) {
      if (factor_name %in% factor_names) {
        reference_levels[[factor_name]] <- levels(dds[[factor_name]])[1]
      } else {
        log_warning(paste("Factor", factor_name, "from step 1 not found in DESeq2 dataset"))
      }
    }
    log_message(paste("Using reference levels:", paste(names(reference_levels), reference_levels, sep = "=", collapse = ", ")))
  } else {
    log_warning("No factors list found in step 1 data - using factor information from DESeq2 dataset")
    for (factor_name in factor_names) {
      if (factor_name != "batch") {  # Skip batch factor
        reference_levels[[factor_name]] <- levels(dds[[factor_name]])[1]
      }
    }
  }
  
  # Get normalized counts from DESeq2
  normCounts <- counts(dds, normalized = TRUE)
  log_message(paste("Extracted normalized counts with dimensions:", nrow(normCounts), "x", ncol(normCounts)))
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    log_message(paste("Applying", args$batchcorrection, "batch correction to normalized counts"))
    normCounts <- apply_batch_correction(
      count_data = normCounts,
      metadata_df = col_data,
      batch_method = args$batchcorrection,
      normalized = TRUE  # These are normalized counts
    )
    log_message("Batch correction completed")
  } else {
    log_message("No batch correction applied to normalized counts")
  }
  
  # Create MDS plot for sample visualization
  log_message("Creating MDS plot for sample visualization")
  with_error_handling({
    create_mds_plot(
      normCounts = normCounts,
      col_metadata = col_data,
      output_file = "mds_plot.html",
      interactive = TRUE
    )
  })
  
  # Process each contrast
  results_list <- list()
  
  for (i in 1:nrow(contrast_df)) {
    contrast_row <- contrast_df[i, ]
    contrast_name <- contrast_row$contrast_name
    log_message(paste("Processing contrast", i, "of", nrow(contrast_df), ":", contrast_name))
    
    # Get contrast results
    contrast_res <- with_error_handling({
      get_contrast_res(step1_data, contrast_row, args)
    })
    
    if (is.null(contrast_res)) {
      log_warning(paste("Failed to process contrast:", contrast_name))
      next
    }
    
    # Add results to expression data
    updated_expDataDf <- with_error_handling({
      add_metadata_to_results(expDataDf, contrast_res, contrast_name)
    })
    
    if (!is.null(updated_expDataDf)) {
      expDataDf <- updated_expDataDf
      log_message(paste("Added results for contrast", contrast_name, "to expression data"))
    }
    
    # Store results
    results_list[[contrast_name]] <- contrast_res
    
    # Force garbage collection after each contrast to manage memory
    gc(verbose = FALSE)
  }
  
  # Prepare final results
  final_results <- list(
    dds = dds,
    normCounts = normCounts,
    expDataDf = expDataDf,
    results = results_list,
    contrasts = contrast_df
  )
  
  # Save results
  if (!is.null(args$output)) {
    output_file <- paste0(args$output, "_deseq2_step2_results.rds")
    log_message(paste("Saving results to", output_file))
    saveRDS(final_results, file = output_file)
  }
  
  # Export reports
  log_message("Exporting reports and visualizations", "STEP")
  
  # Export gene expression table
  gene_exp_file <- export_deseq_results(expDataDf, args$output, output_name="gene_exp_table")
  log_message(paste("Exported gene expression table to", gene_exp_file), "INFO")
  
  # Export normalized counts in GCT format
  count_files <- export_normalized_counts(normCounts, args$output, threshold=args$rpkm_cutoff)
  log_message(paste("Exported normalized counts to", count_files$all_counts), "INFO")
  
  # Create interactive MDS plot
  mds_file <- get_output_filename(args$output, "mds_plot", "html")
  generate_mds_plot_html(
    normCounts, 
    mds_file,
    metadata=col_data,
    color_by=factor_names[1],  # Use first factor for coloring
    title="MDS Plot of Samples"
  )
  log_message(paste("Created interactive MDS plot at", mds_file), "INFO")
  
  # Create visualizations for each contrast
  for (contrast_name in names(results_list)) {
    contrast_results <- results_list[[contrast_name]]
    contrast_prefix <- paste0(args$output, "_", make.names(contrast_name))
    
    # Create MA plot for this contrast
    ma_plot <- function() {
      DESeq2::plotMA(contrast_results, main=paste("MA Plot -", contrast_name))
    }
    ma_files <- save_plot(ma_plot, contrast_prefix, "ma_plot")
    log_message(paste("Created MA plot for", contrast_name), "INFO")
    
    # Create volcano plot for this contrast
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      results_df <- as.data.frame(contrast_results)
      results_df$significance <- ifelse(
        results_df$padj < args$fdr & abs(results_df$log2FoldChange) > args$lfcthreshold,
        "Significant", 
        "Not Significant"
      )
      
      volcano_plot <- ggplot2::ggplot(
        results_df, 
        ggplot2::aes(
          x = log2FoldChange, 
          y = -log10(pvalue),
          color = significance
        )
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
        ggplot2::labs(
          title = paste("Volcano Plot -", contrast_name),
          x = "log2 Fold Change",
          y = "-log10(p-value)"
        ) +
        ggplot2::theme_minimal()
      
      volcano_files <- save_plot(volcano_plot, contrast_prefix, "volcano_plot")
      log_message(paste("Created volcano plot for", contrast_name), "INFO")
    }
  }
  
  # Verify output files
  verify_outputs(args$output, "lrt_step2", fail_on_missing=FALSE)
  
  log_message("DESeq2 LRT Step 2 workflow completed successfully", "SUCCESS")
  
  return(final_results)
} 