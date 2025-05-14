#!/usr/bin/env Rscript

# --- DESeq2 analysis functions ---

# Main function to run DESeq2 analysis
run_deseq2 <- function(count_data, metadata, design_formula, reduced_formula, args) {
  logger::info("Running DESeq2 analysis with LRT")
  
  # Create DESeqDataSet object
  dds <- try(DESeqDataSetFromMatrix(
    countData = count_data,
    colData = metadata,
    design = design_formula
  ))
  
  if (inherits(dds, "try-error")) {
    stop("Failed to create DESeqDataSet from matrix")
  }
  
  # Run DESeq with LRT
  logger::info(paste("Running DESeq with full model:", deparse(design_formula)))
  logger::info(paste("Reduced model:", deparse(reduced_formula)))
  
  dds <- try(DESeq(
    dds,
    test = "LRT",
    reduced = reduced_formula,
    parallel = (args$threads > 1),
    BPPARAM = MulticoreParam(args$threads)
  ))
  
  if (inherits(dds, "try-error")) {
    stop("DESeq2 analysis failed")
  }
  
  # Extract results with FDR threshold
  results <- try(results(
    dds,
    alpha = args$fdr,
    lfcThreshold = args$lfcthreshold
  ))
  
  if (inherits(results, "try-error")) {
    stop("Failed to extract results from DESeq2 analysis")
  }
  
  # Sort results by p-value
  results <- results[order(results$pvalue), ]
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Get VST transformed data
  vst_data <- try(vst(dds, blind = FALSE))
  
  if (inherits(vst_data, "try-error")) {
    logger::warn("VST transformation failed, will attempt to use rlog instead")
    vst_data <- try(rlog(dds, blind = FALSE))
    
    if (inherits(vst_data, "try-error")) {
      logger::error("Both VST and rlog transformations failed")
      vst_data <- NULL
    }
  }
  
  return(list(
    dds = dds,
    results = results,
    normalized_counts = normalized_counts,
    vst_data = vst_data
  ))
}

# Function to export results from DESeq2 analysis
export_results <- function(deseq_results, expression_data, metadata, args, batch_warning = NULL, rpkm_filtered_count = NULL) {
  log_message("Exporting DESeq2 LRT Step 1 results", "STEP")
  
  # Ensure output directory exists
  output_dir <- dirname(args$output_prefix)
  if (!is.null(output_dir) && output_dir != "" && output_dir != ".") {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # Check if output_prefix is missing, use output as fallback
  if (is.null(args$output_prefix)) {
    log_message("output_prefix is NULL, using output instead", "WARNING")
    args$output_prefix <- args$output
  }
  
  # Check if deseq_results is NULL
  if (is.null(deseq_results)) {
    log_message("DESeq2 results object is NULL, unable to export results", "ERROR")
    return(NULL)
  }
  
  # Safely extract result components
  dds <- if (is.list(deseq_results) && !is.null(deseq_results$dds)) deseq_results$dds else NULL
  lrt_res <- if (is.list(deseq_results)) {
    if (!is.null(deseq_results$lrt_res)) {
      deseq_results$lrt_res
    } else if (!is.null(deseq_results$results)) {
      deseq_results$results  # Backward compatibility
    } else {
      NULL
    }
  } else {
    NULL
  }
  
  # Safely extract other components
  norm_counts <- if (is.list(deseq_results)) {
    if (!is.null(deseq_results$normCounts)) {
      deseq_results$normCounts
    } else if (!is.null(deseq_results$normalized_counts)) {
      deseq_results$normalized_counts  # Backward compatibility
    } else if (!is.null(dds)) {
      counts(dds, normalized = TRUE)
    } else {
      NULL
    }
  } else {
    NULL
  }
  
  vst_data <- if (is.list(deseq_results) && !is.null(deseq_results$vst_data)) deseq_results$vst_data else NULL
  contrasts <- if (is.list(deseq_results) && !is.null(deseq_results$contrasts)) deseq_results$contrasts else NULL
  
  design_formula <- if (is.list(deseq_results) && !is.null(deseq_results$design_formula)) {
    deseq_results$design_formula
  } else if (!is.null(dds)) {
    design(dds)
  } else {
    NULL
  }
  
  reduced_formula <- if (is.list(deseq_results) && !is.null(deseq_results$reduced_formula)) {
    deseq_results$reduced_formula
  } else {
    NULL
  }
  
  # Save DESeq2 dataset and results as RDS files if available
  if (!is.null(dds)) {
    export_deseq_object(dds, args$output_prefix, "dds")
  }
  
  if (!is.null(deseq_results)) {
    export_deseq_object(deseq_results, args$output_prefix, "contrasts")
  }
  
  # Export contrasts table if available
  if (!is.null(contrasts)) {
    export_deseq_results(contrasts, args$output_prefix, output_name="contrasts_table")
  }
  
  # Export LRT results if available
  if (!is.null(lrt_res)) {
    export_deseq_results(lrt_res, args$output_prefix, output_name="gene_exp_table")
  }
  
  # Export normalized counts in GCT format if available
  if (!is.null(dds)) {
    export_normalized_counts(dds, args$output_prefix, threshold=args$rpkm_cutoff)
  } else if (!is.null(norm_counts)) {
    # If dds is not available but we have normalized counts
    write_gct_file(norm_counts, get_output_filename(args$output_prefix, "counts_all", "gct"))
    
    # Apply filtering if threshold is provided
    if (!is.null(args$rpkm_cutoff)) {
      # Filter by row means
      row_means <- rowMeans(norm_counts)
      filtered_counts <- norm_counts[row_means >= args$rpkm_cutoff, , drop = FALSE]
      write_gct_file(filtered_counts, get_output_filename(args$output_prefix, "counts_filtered", "gct"))
    }
  }
  
  # Generate MDS plot HTML if vst_data is available
  if (!is.null(vst_data) && is.object(vst_data)) {
    # Check if we can extract the assay data
    vst_matrix <- try(assay(vst_data), silent = TRUE)
    if (!inherits(vst_matrix, "try-error") && !is.null(vst_matrix)) {
      generate_mds_plot_html(
        vst_matrix, 
        get_output_filename(args$output_prefix, "mds_plot", "html"),
        metadata = metadata, 
        color_by = if ("condition" %in% colnames(metadata)) "condition" else names(metadata)[1]
      )
    }
  }
  
  # If batch correction was applied, create batch-corrected MDS plot
  if (!is.null(args$batchcorrection) && args$batchcorrection != "none" && 
      !is.null(metadata) && is.data.frame(metadata) && "batch" %in% colnames(metadata)) {
    log_message("Generating batch-corrected MDS plot", "INFO")
    # Use the corrected VST data if available
    if (is.list(deseq_results) && !is.null(deseq_results$vst_data_corrected)) {
      vst_corrected_matrix <- try(assay(deseq_results$vst_data_corrected), silent = TRUE)
      if (!inherits(vst_corrected_matrix, "try-error") && !is.null(vst_corrected_matrix)) {
        generate_mds_plot_html(
          vst_corrected_matrix,
          get_output_filename(args$output_prefix, "mds_plot_corrected", "html"),
          metadata = metadata,
          color_by = if ("condition" %in% colnames(metadata)) "condition" else names(metadata)[1]
        )
      }
    }
  }
  
  # Generate summary markdown if lrt_res is available
  if (!is.null(lrt_res)) {
    design_formula_text <- if (!is.null(design_formula)) {
      if (is.language(design_formula)) deparse(design_formula) else as.character(design_formula)
    } else {
      "Not available"
    }
    
    reduced_formula_text <- if (!is.null(reduced_formula)) {
      if (is.language(reduced_formula)) deparse(reduced_formula) else as.character(reduced_formula)
    } else {
      "Not available"
    }
    
    generate_deseq_summary(
      lrt_res, 
      get_output_filename(args$output_prefix, "lrt_result", "md"),
      title = "DESeq2 LRT Analysis Summary",
      parameters = list(
        "Design formula" = design_formula_text,
        "Reduced formula" = reduced_formula_text,
        "Batch correction" = args$batchcorrection,
        "RPKM cutoff" = args$rpkm_cutoff,
        "Filtered genes" = rpkm_filtered_count
      )
    )
  }
  
  # Create visualizations if dds and lrt_res are available
  if (!is.null(dds) && !is.null(lrt_res)) {
    export_visualizations(
      dds, 
      lrt_res, 
      args$output_prefix,
      transformed_counts = vst_data,
      metadata = metadata,
      color_by = if (!is.null(metadata) && "condition" %in% colnames(metadata)) 
        "condition" else NULL
    )
  }
  
  # Create alignment stats barchart if dds is available
  if (!is.null(dds) && requireNamespace("ggplot2", quietly = TRUE)) {
    # Generate alignment stats from metadata and counts
    alignment_stats <- data.frame(
      sample = colnames(counts(dds)),
      read_count = colSums(counts(dds))
    )
    
    # Create the plot
    alignment_plot <- ggplot2::ggplot(
      alignment_stats, 
      ggplot2::aes(x = sample, y = read_count)
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = "Alignment Statistics",
        x = "Sample",
        y = "Read Count"
      )
    
    # Save plot
    save_plot(alignment_plot, "", "alignment_stats_barchart", formats = "png")
  }
  
  log_message("Results exported successfully", "SUCCESS")
  return(TRUE)
}

# Function to run DESeq2 LRT analysis
run_deseq2_lrt <- function(count_data, sample_metadata, design_formula, reduced_formula, 
                           batchcorrection = "none", threads = 1) {
  require(DESeq2)
  
  # Create DESeqDataSet
  message("Creating DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_metadata,
    design = design_formula
  )
  
  # Apply batch correction if requested
  if (batchcorrection == "combatseq") {
    message("Applying ComBat-Seq batch correction...")
    require(sva)
    
    # Check if 'batch' column exists in sample metadata
    if (!("batch" %in% colnames(sample_metadata))) {
      stop("Batch correction requested but 'batch' column not found in metadata")
    }
    
    # Get batch information
    batch <- sample_metadata$batch
    
    # Get model matrix from full design
    group <- sample_metadata[, sapply(sample_metadata, is.factor)]
    mod <- model.matrix(design_formula, data = group)
    
    # Apply ComBat-Seq
    corrected_counts <- ComBat_seq(
      counts = count_data,
      batch = batch,
      group = group
    )
    
    # Update DESeqDataSet with corrected counts
    dds <- DESeqDataSetFromMatrix(
      countData = corrected_counts,
      colData = sample_metadata,
      design = design_formula
    )
  }
  
  # Set up parallel processing
  if (threads > 1) {
    message(paste("Using", threads, "threads for parallel processing"))
    BiocParallel::register(BiocParallel::MulticoreParam(threads))
  }
  
  # Run DESeq2 with LRT
  message("Running DESeq2 LRT analysis...")
  dds <- DESeq(
    dds,
    test = "LRT",
    reduced = reduced_formula,
    parallel = (threads > 1)
  )
  
  return(dds)
}

# Function to extract results from DESeq2 analysis
extract_deseq2_results <- function(dds, fdr = 0.1, lfcthreshold = 0.59, use_lfc_thresh = FALSE) {
  require(DESeq2)
  
  # Determine if we're dealing with LRT or Wald test
  is_lrt <- mcols(dds)$modelMatrixType == "LRT"
  
  if (is_lrt) {
    message(paste("Extracting DESeq2 LRT results with FDR =", fdr))
    
    # For LRT test, don't use lfcThreshold parameter
    res <- results(
      dds,
      alpha = fdr,
      altHypothesis = "greaterAbs",
      pAdjustMethod = "BH",
      independentFiltering = TRUE,
      cooksCutoff = TRUE,
      parallel = FALSE
    )
  } else {
    # Wald test logic
    if (use_lfc_thresh) {
      message(paste("Extracting DESeq2 Wald results with FDR =", fdr, 
                    "and LFC threshold =", lfcthreshold, "as hypothesis threshold"))
      
      # If use_lfc_thresh is TRUE, use lfcthreshold directly in hypothesis test
      res <- results(
        dds,
        alpha = fdr,
        lfcThreshold = lfcthreshold,
        altHypothesis = "greaterAbs",
        pAdjustMethod = "BH",
        independentFiltering = TRUE,
        cooksCutoff = TRUE,
        parallel = FALSE
      )
    } else {
      message(paste("Extracting DESeq2 Wald results with FDR =", fdr, 
                    "and post-filtering with LFC threshold =", lfcthreshold))
      
      # If use_lfc_thresh is FALSE, use lfcThreshold = 0 and filter later
      res <- results(
        dds,
        alpha = fdr,
        lfcThreshold = 0,
        altHypothesis = "greaterAbs",
        pAdjustMethod = "BH",
        independentFiltering = TRUE,
        cooksCutoff = TRUE,
        parallel = FALSE
      )
    }
  }
  
  # Convert to data frame for easier manipulation
  res_df <- as.data.frame(res)
  
  # Add gene names if available
  if (!is.null(mcols(dds)$symbol)) {
    res_df$symbol <- mcols(dds)$symbol
  }
  
  # Apply post-filtering for LFC if use_lfc_thresh is FALSE and it's not an LRT test
  if (!is_lrt && !use_lfc_thresh && lfcthreshold > 0) {
    message(paste("Post-filtering results with |log2FoldChange| >=", lfcthreshold))
    
    # Mark genes not meeting threshold as non-significant
    significant_genes <- !is.na(res_df$padj) & res_df$padj < fdr & abs(res_df$log2FoldChange) >= lfcthreshold
    
    # Create a new column for LFC-filtered p-values
    res_df$padj_lfc_filtered <- res_df$padj
    res_df$padj_lfc_filtered[!significant_genes] <- NA
    
    # Create a filtering status column
    res_df$lfc_significant <- significant_genes
    
    message(paste("After LFC filtering:", sum(significant_genes, na.rm = TRUE), "of", 
                  sum(!is.na(res_df$padj) & res_df$padj < fdr), "significant genes remain"))
  }
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  return(res_df)
}

# Function to apply limma batch correction on the results
apply_limma_batch_correction <- function(dds, sample_metadata) {
  require(limma)
  
  message("Applying limma batch correction to normalized counts...")
  
  # Check if 'batch' column exists in sample metadata
  if (!("batch" %in% colnames(sample_metadata))) {
    stop("Batch correction requested but 'batch' column not found in metadata")
  }
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Create design matrix
  model <- model.matrix(~ 0 + group, data = sample_metadata)
  colnames(model) <- levels(sample_metadata$group)
  
  # Create batch variable
  batch <- sample_metadata$batch
  
  # Apply removeBatchEffect
  corrected_counts <- removeBatchEffect(
    normalized_counts,
    batch = batch,
    design = model
  )
  
  return(corrected_counts)
}

# Function to validate essential arguments
validate_essential_args <- function(args) {
  # Check if args is a valid list
  if (!is.list(args)) {
    stop("Arguments must be provided as a list")
  }
  
  # Check essential arguments
  essential_args <- c("meta", "design", "reduced", "input", "name")
  missing_args <- essential_args[!essential_args %in% names(args)]
  
  if (length(missing_args) > 0) {
    stop(paste("Missing essential arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Check if any essential arguments are NULL or empty
  empty_args <- c()
  for (arg in essential_args) {
    if (is.null(args[[arg]]) || (is.character(args[[arg]]) && length(args[[arg]]) == 0)) {
      empty_args <- c(empty_args, arg)
    }
  }
  
  if (length(empty_args) > 0) {
    stop(paste("These essential arguments are empty:", paste(empty_args, collapse=", ")))
  }
  
  # Check that input and name have equal lengths
  if (length(args$input) != length(args$name)) {
    logger::warn(paste("Mismatch between number of input files (", length(args$input), 
                        ") and sample names (", length(args$name), ")"))
    
    # If more samples than inputs, trim sample names
    if (length(args$name) > length(args$input)) {
      logger::warn("Trimming extra sample names to match input files")
      args$name <- args$name[1:length(args$input)]
    }
    
    # If more inputs than samples, trim input files
    if (length(args$input) > length(args$name)) {
      logger::warn("Trimming extra input files to match sample names")
      args$input <- args$input[1:length(args$name)]
    }
  }
  
  logger::info("Essential arguments validated successfully")
  return(args)
}

# Main function to perform DESeq2 LRT analysis
run_deseq2_lrt_analysis <- function(args) {
  # Log the start of the analysis
  logger::info("Starting DESeq2 analysis with LRT")
  logger::debug("Arguments:", args)
  
  # Validate essential arguments
  validate_essential_args(args)
  
  # Load the metadata file
  metadata <- load_metadata(args$meta)
  
  # Load expression data
  if (is.null(args$input) || length(args$input) == 0) {
    stop("No input expression files provided")
  }
  
  if (is.null(args$name) || length(args$name) == 0) {
    stop("No sample names provided")
  }
  
  # Make sure name doesn't contain non-sample names
  if (length(args$name) > length(args$input)) {
    logger::warn("More sample names than input files provided. Trimming extra names.")
    args$name <- args$name[1:length(args$input)]
  }
  
  # Now proceed with loading the data
  expression_data <- load_expression_data(args$input, args$name)
  
  # Extract read counts
  counts_cols <- grep(READ_COL, colnames(expression_data), ignore.case = TRUE, value = TRUE)
  if (length(counts_cols) == 0) {
    stop("No count columns found in expression data")
  }
  
  count_data <- expression_data[, c(INTERSECT_BY, counts_cols)]
  
  # Set gene IDs as row names
  row.names(count_data) <- count_data[[INTERSECT_BY]]
  count_data[[INTERSECT_BY]] <- NULL
  
  # Check if there are any zero rows
  zero_rows <- rowSums(count_data) == 0
  if (any(zero_rows)) {
    logger::info(paste("Removing", sum(zero_rows), "rows with zero counts"))
    count_data <- count_data[!zero_rows, , drop = FALSE]
  }
  
  # Set up DESeq2 design from formula string
  if (is.null(args$design) || args$design == "") {
    stop("No valid design formula provided")
  }
  design_formula <- as.formula(args$design)
  
  # Set up reduced design for LRT
  if (is.null(args$reduced) || args$reduced == "") {
    stop("No valid reduced design formula provided")
  }
  reduced_formula <- as.formula(args$reduced)
  
  # Prepare metadata and ensure it matches count data samples
  metadata <- prepare_metadata_for_deseq(metadata, colnames(count_data))
  
  # Run DESeq2 analysis
  deseq_results <- run_deseq2(count_data, metadata, design_formula, reduced_formula, args)
  
  # Export results
  export_results(deseq_results, expression_data, metadata, args)
} 