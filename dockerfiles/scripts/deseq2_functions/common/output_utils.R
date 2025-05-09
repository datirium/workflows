#!/usr/bin/env Rscript
#
# Streamlined output utility functions for DESeq2 workflows
# This file contains consolidated output functions to reduce redundancy
#

#' Create standardized filename for output
#'
#' @param prefix Output prefix specified by the user
#' @param stem Type of output (e.g., "ma_plot", "gene_exp_table") 
#' @param extension File extension without the dot (e.g., "png", "pdf")
#' @return Standardized filename string
#' @export
get_output_filename <- function(prefix, stem, extension) {
  # Clean prefix (remove trailing slashes and dots)
  prefix <- gsub("[/\\.]$", "", prefix)
  
  # Special cases based on CWL expected outputs
  if (stem == "alignment_stats_barchart" && extension == "png") {
    return("alignment_stats_barchart.png")
  } 
  
  # Normal case: prefix_stem.extension
  if (stem == "") {
    return(paste0(prefix, ".", extension))
  } else {
    return(paste0(prefix, "_", stem, ".", extension))
  }
}

#' Save a plot in multiple formats
#'
#' @param plot A ggplot2 object or function that creates a plot
#' @param prefix Output prefix
#' @param stem Middle part of the filename
#' @param formats Vector of formats to save (default: c("png", "pdf"))
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param ... Additional parameters passed to plot device
#' @return List of created filenames
#' @export
save_plot <- function(plot, prefix, stem, formats = c("png", "pdf"), 
                      width = 8, height = 6, ...) {
  results <- list()
  
  for (format in formats) {
    file_path <- get_output_filename(prefix, stem, format)
    
    # For ggplot2 plots
    if (inherits(plot, "ggplot")) {
      ggplot2::ggsave(
        filename = file_path,
        plot = plot,
        width = width,
        height = height,
        dpi = 300,
        units = "in",
        device = if(format == "pdf") cairo_pdf else NULL,
        ...
      )
    } else {
      # For base R plots or function calls
      if (format == "png") {
        png(file_path, width = width, height = height, units = "in", res = 300)
      } else if (format == "pdf") {
        pdf(file_path, width = width, height = height)
      }
      
      if (is.function(plot)) plot() else print(plot)
      dev.off()
    }
    
    results[[format]] <- file_path
    message(paste("Saved", format, "plot to", file_path))
  }
  
  return(results)
}

#' Write data to a GCT file for GSEA
#'
#' @param data Expression data matrix or data frame (genes in rows, samples in columns)
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_gct_file <- function(data, file_path) {
  # Make sure it's a matrix or data frame
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  
  # Extract dimensions
  n_rows <- nrow(data)
  n_cols <- ncol(data)
  
  # Ensure we have row names
  if (is.null(rownames(data))) {
    rownames(data) <- paste0("GENE_", 1:n_rows)
  }
  
  # Open connection to file
  con <- file(file_path, "w")
  
  # Write GCT header
  writeLines("#1.2", con)
  writeLines(paste(n_rows, n_cols, sep = "\t"), con)
  
  # Write column headers
  col_headers <- c("NAME", "Description", colnames(data))
  writeLines(paste(col_headers, collapse = "\t"), con)
  
  # Write data
  for (i in 1:n_rows) {
    # Create description (use "na" as a placeholder)
    description <- "na"
    
    # Get data for the row
    row_data <- data[i, ]
    row_name <- rownames(data)[i]
    
    # Format the output line
    row_str <- paste(c(row_name, description, as.character(row_data)), collapse = "\t")
    writeLines(row_str, con)
  }
  
  close(con)
  message(paste("Created GCT file:", file_path))
  
  return(file_path)
}

#' Create a CLS file for GSEA
#'
#' @param sample_classes Vector of sample class labels 
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_cls_file <- function(sample_classes, file_path) {
  # Get unique classes
  unique_classes <- unique(sample_classes)
  n_classes <- length(unique_classes)
  n_samples <- length(sample_classes)
  
  # Create header line
  header <- paste(n_samples, n_classes, 1, sep = " ")
  
  # Create class names line
  class_names <- paste("#", paste(unique_classes, collapse = " "))
  
  # Create sample assignments line (0-based indices)
  class_indices <- sapply(sample_classes, function(x) which(unique_classes == x) - 1)
  assignments <- paste(class_indices, collapse = " ")
  
  # Write to file
  writeLines(c(header, class_names, assignments), file_path)
  message(paste("Created CLS file:", file_path))
  
  return(file_path)
}

#' Write DESeq2 results to a TSV file
#'
#' @param results DESeq2 results object or data frame
#' @param file_path Output file path
#' @param add_rank Whether to add a rank column (log2FC * -log10(pvalue))
#' @return The path to the created file
#' @export
write_deseq_results <- function(results, file_path, add_rank = TRUE) {
  # Convert to data frame if needed
  if (!is.data.frame(results)) {
    results_df <- as.data.frame(results)
  } else {
    results_df <- results
  }
  
  # Add rank column for GSEA if requested
  if (add_rank && "log2FoldChange" %in% colnames(results_df) && "pvalue" %in% colnames(results_df)) {
    results_df$rank <- sign(results_df$log2FoldChange) * -log10(results_df$pvalue)
  }
  
  # Sort by adjusted p-value if available
  if ("padj" %in% colnames(results_df)) {
    results_df <- results_df[order(results_df$padj), ]
  }
  
  # Write to TSV file
  write.table(results_df, file = file_path, sep = "\t", quote = FALSE, row.names = TRUE)
  message(paste("Created results TSV file:", file_path))
  
  return(file_path)
}

#' Create a simple markdown summary file
#'
#' @param content Character vector of content lines
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_markdown_summary <- function(content, file_path) {
  writeLines(content, file_path)
  message(paste("Created markdown summary:", file_path))
  
  return(file_path)
}

#' Generate a standard summary markdown for DESeq2 results
#'
#' @param deseq_results DESeq2 results object
#' @param file_path Output file path 
#' @param title Title for the summary
#' @param parameters List of parameters to include
#' @param fdr_cutoff FDR cutoff for counting significant genes
#' @return Path to the created markdown file
#' @export
generate_deseq_summary <- function(deseq_results, file_path, title = "DESeq2 Analysis Summary",
                                 parameters = list(), fdr_cutoff = 0.1) {
  # Initialize content
  lines <- c(
    paste("#", title),
    "",
    paste("Analysis performed on", Sys.Date()),
    "",
    "## Parameters",
    ""
  )
  
  # Add parameters
  for (name in names(parameters)) {
    lines <- c(lines, paste("*", name, ":", parameters[[name]]))
  }
  
  # Get summary statistics
  if (is.data.frame(deseq_results)) {
    res_df <- deseq_results
  } else {
    res_df <- as.data.frame(deseq_results)
  }
  
  total_genes <- nrow(res_df)
  sig_genes <- sum(res_df$padj < fdr_cutoff, na.rm = TRUE)
  
  # Add up/down regulation if log2FoldChange exists
  if ("log2FoldChange" %in% colnames(res_df)) {
    up_genes <- sum(res_df$padj < fdr_cutoff & res_df$log2FoldChange > 0, na.rm = TRUE)
    down_genes <- sum(res_df$padj < fdr_cutoff & res_df$log2FoldChange < 0, na.rm = TRUE)
    
    # Add result summary
    lines <- c(lines,
      "",
      "## Results",
      "",
      paste("* Total genes analyzed:", total_genes),
      paste("* Significant genes (FDR <", fdr_cutoff, "):", sig_genes),
      paste("* Up-regulated genes:", up_genes),
      paste("* Down-regulated genes:", down_genes)
    )
  } else {
    # LRT or other test without fold change
    lines <- c(lines,
      "",
      "## Results",
      "",
      paste("* Total genes analyzed:", total_genes),
      paste("* Significant genes (FDR <", fdr_cutoff, "):", sig_genes)
    )
  }
  
  # Write to file
  write_markdown_summary(lines, file_path)
  
  return(file_path)
}

#' Generate interactive MDS plot HTML
#'
#' @param transformed_data Transformed count data (e.g., vst or rlog)
#' @param file_path Output file path
#' @param metadata Sample metadata (optional)
#' @param color_by Column name in metadata to use for coloring (optional)
#' @param title Plot title
#' @return Path to created HTML file
#' @export
generate_mds_plot_html <- function(transformed_data, file_path, metadata = NULL,
                                color_by = NULL, title = "MDS Plot") {
  # Check required packages
  if (!requireNamespace("plotly", quietly = TRUE) || 
      !requireNamespace("htmlwidgets", quietly = TRUE) ||
      !requireNamespace("limma", quietly = TRUE)) {
    warning("Cannot create MDS plot: missing required packages (plotly, htmlwidgets, or limma)")
    return(NULL)
  }
  
  # Calculate MDS
  mds_data <- limma::plotMDS(transformed_data, plot = FALSE)
  
  # Create data frame for plotting
  mds_df <- data.frame(
    x = mds_data$x,
    y = mds_data$y,
    sample = colnames(transformed_data)
  )
  
  # Add metadata if provided
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    # Match samples to metadata
    if (all(mds_df$sample %in% rownames(metadata))) {
      idx <- match(mds_df$sample, rownames(metadata))
      mds_df$color <- metadata[idx, color_by]
    } else {
      warning("Sample names in transformed data don't match metadata rownames")
    }
  }
  
  # Create plotly figure
  p <- plotly::plot_ly(mds_df, x = ~x, y = ~y, text = ~sample,
                     mode = "markers", type = "scatter",
                     marker = list(size = 10))
  
  if ("color" %in% colnames(mds_df)) {
    p <- plotly::add_trace(p, color = ~color)
  }
  
  p <- plotly::layout(p, title = title,
                    xaxis = list(title = "Leading logFC dim 1"),
                    yaxis = list(title = "Leading logFC dim 2"))
  
  # Save to HTML file
  htmlwidgets::saveWidget(p, file_path, selfcontained = TRUE)
  message(paste("Generated MDS plot HTML at", file_path))
  
  return(file_path)
}

#' Verify all outputs for a workflow exist
#'
#' @param output_prefix Output prefix used for files
#' @param workflow_type Type of workflow: "deseq", "lrt_step1", or "lrt_step2"
#' @param fail_on_missing Whether to stop if files are missing
#' @return Logical indicating success
#' @export
verify_outputs <- function(output_prefix, workflow_type, fail_on_missing = FALSE) {
  # Clean prefix
  prefix <- gsub("[/\\.]$", "", output_prefix)
  
  # Define expected files based on workflow type
  expected_files <- switch(workflow_type,
    "deseq" = c(
      paste0(prefix, "_report.tsv"),
      paste0(prefix, "_summary.md"),
      paste0(prefix, "_counts_all.gct"),
      paste0(prefix, "_counts_filtered.gct"),
      paste0(prefix, "_phenotypes.cls"),
      paste0(prefix, "_ma_plot.png"),
      paste0(prefix, "_expression_heatmap.png"),
      paste0(prefix, "_pca_plot.png"),
      paste0(prefix, "_mds_plot.html")
    ),
    "lrt_step1" = c(
      paste0(prefix, "_contrasts_table.tsv"),
      paste0(prefix, "_contrasts.rds"),
      paste0(prefix, "_gene_exp_table.tsv"),
      paste0(prefix, "_mds_plot.html"),
      paste0(prefix, "_counts_all.gct"),
      paste0(prefix, "_counts_filtered.gct"),
      paste0(prefix, "_lrt_result.md"),
      "alignment_stats_barchart.png"
    ),
    "lrt_step2" = c(
      paste0(prefix, "_gene_exp_table.tsv"),
      paste0(prefix, "_mds_plot.html"),
      paste0(prefix, "_counts_all.gct"),
      paste0(prefix, "_counts_filtered.gct")
    ),
    stop("Unknown workflow type: ", workflow_type)
  )
  
  # Check for missing files
  missing_files <- expected_files[!file.exists(expected_files)]
  
  # Report results
  if (length(missing_files) > 0) {
    msg <- paste("Missing expected output files:", 
                paste(missing_files, collapse = ", "))
    
    if (fail_on_missing) {
      stop(msg)
    } else {
      warning(msg)
      return(FALSE)
    }
  }
  
  message("All expected output files were created successfully")
  return(TRUE) 