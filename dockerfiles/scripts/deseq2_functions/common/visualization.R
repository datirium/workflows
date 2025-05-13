#!/usr/bin/env Rscript
#
# Common visualization functions for DESeq2 analysis
#

#' Create an MDS plot for sample visualization
#'
#' @param normCounts Normalized count matrix
#' @param col_metadata Column metadata with sample information
#' @param output_file Output HTML file name
#' @param interactive Whether to create an interactive plot
#' @param color_by Column name to color points by (optional)
#' @param label_by Column name to label points by (optional)
#' @param title Plot title (optional)
#' @return MDS plot object
#' @export
create_mds_plot <- function(normCounts, col_metadata, 
                           output_file = "mds_plot.html", 
                           interactive = TRUE,
                           color_by = NULL, 
                           label_by = NULL, 
                           title = "Multi-dimensional scaling plot") {
  log_message("Creating MDS plot...", "STEP")
  
  # Check dimensions
  if (ncol(normCounts) < 2) {
    log_message("Cannot create MDS plot with fewer than 2 samples", "WARNING")
    return(NULL)
  }
  
  if (nrow(normCounts) < 10) {
    log_message("MDS plot may be unreliable with fewer than 10 genes", "WARNING")
  }
  
  # Calculate distances between samples using correlation
  log_message("Calculating sample distances...")
  sample_dists <- as.dist(1 - cor(normCounts, method = "spearman"))
  
  # Perform MDS
  log_message("Performing MDS...")
  mds_result <- cmdscale(sample_dists, k = 2)
  
  # Create a data frame for plotting
  mds_df <- as.data.frame(mds_result)
  colnames(mds_df) <- c("MDS1", "MDS2")
  
  # Add sample names
  mds_df$Sample <- rownames(mds_df)
  
  # Add metadata if available
  if (!is.null(col_metadata) && nrow(col_metadata) > 0) {
    # Ensure rownames match
    common_samples <- intersect(rownames(mds_df), rownames(col_metadata))
    if (length(common_samples) < nrow(mds_df)) {
      log_message(paste("Only", length(common_samples), "of", nrow(mds_df), "samples found in metadata"), "WARNING")
    }
    
    if (length(common_samples) > 0) {
      metadata_subset <- col_metadata[common_samples, , drop = FALSE]
      mds_df <- cbind(mds_df, metadata_subset)
    }
  }
  
  # Set color_by if not specified but available
  if (is.null(color_by) && ncol(mds_df) > 3) {
    # Find first factor column
    factor_cols <- sapply(mds_df, is.factor)
    if (any(factor_cols)) {
      color_by <- names(factor_cols)[which(factor_cols)[1]]
      log_message(paste("Automatically coloring by", color_by), "INFO")
    } else {
      color_by <- "Sample"
    }
  } else if (is.null(color_by)) {
    color_by <- "Sample"
  }
  
  # Set label_by if not specified
  if (is.null(label_by)) {
    label_by <- "Sample"
  }
  
  # Create static plot with ggplot2
  log_message("Creating static plot...")
  p <- ggplot2::ggplot(mds_df, ggplot2::aes_string(x = "MDS1", y = "MDS2", color = color_by, label = label_by)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Samples colored by", color_by),
      x = "MDS Dimension 1",
      y = "MDS Dimension 2"
    ) +
    ggplot2::theme_bw()
  
  # Create interactive plot if requested
  if (interactive) {
    log_message("Creating interactive plot...")
    
    # Create interactive Plotly version
    p_interactive <- plotly::ggplotly(p, tooltip = c("label", "MDS1", "MDS2", color_by))
    
    # Save the plot
    htmlwidgets::saveWidget(p_interactive, output_file, selfcontained = TRUE)
    log_message(paste("Interactive MDS plot saved to", output_file), "SUCCESS")
  } else {
    # Save static plot
    ggplot2::ggsave(output_file, p, width = 8, height = 6)
    log_message(paste("Static MDS plot saved to", output_file), "SUCCESS")
  }
  
  return(p)
}

#' Export GCT file format for visualization and analysis tools
#'
#' @param expr_data Expression data matrix
#' @param row_metadata Row metadata
#' @param col_metadata Column metadata
#' @param output_prefix Output file prefix
#' @param filename Specific filename (overrides prefix if provided)
#' @return Path to the saved GCT file
#' @export
export_gct_data <- function(expr_data, row_metadata, col_metadata, 
                           output_prefix = NULL, filename = NULL) {
  log_message("Exporting data in GCT format...", "STEP")
  
  # Set output filename
  if (is.null(filename)) {
    filename <- paste0(output_prefix, "_counts_all.gct")
  }
  
  # Check dimensions
  if (nrow(expr_data) == 0 || ncol(expr_data) == 0) {
    report_error(
      "Cannot export empty expression matrix",
      details = paste("Expression matrix has dimensions", nrow(expr_data), "x", ncol(expr_data)),
      recommendations = "Check that filtering has not removed all rows/columns"
    )
  }
  
  # Ensure row metadata has correct dimensions
  if (nrow(row_metadata) != nrow(expr_data)) {
    log_message(paste("Row metadata has", nrow(row_metadata), "rows, but expression data has", 
                     nrow(expr_data), "rows"), "WARNING")
    
    # Try to match by rownames
    common_rows <- intersect(rownames(expr_data), rownames(row_metadata))
    if (length(common_rows) == 0) {
      report_error(
        "No common rows between expression data and row metadata",
        details = "Row metadata and expression data have no common row names",
        recommendations = c(
          "Check that row names are consistent between datasets",
          "Ensure that expression data and metadata are from the same source"
        )
      )
    }
    
    # Subset to common rows
    log_message(paste("Using", length(common_rows), "common rows between expression data and row metadata"), "INFO")
    expr_data <- expr_data[common_rows, , drop = FALSE]
    row_metadata <- row_metadata[common_rows, , drop = FALSE]
  }
  
  # Ensure column metadata has correct dimensions
  if (nrow(col_metadata) != ncol(expr_data)) {
    log_message(paste("Column metadata has", nrow(col_metadata), "rows, but expression data has", 
                     ncol(expr_data), "columns"), "WARNING")
    
    # Try to match by rownames/colnames
    common_cols <- intersect(colnames(expr_data), rownames(col_metadata))
    if (length(common_cols) == 0) {
      report_error(
        "No common columns between expression data and column metadata",
        details = "Column metadata and expression data have no common column/row names",
        recommendations = c(
          "Check that column names and row names are consistent between datasets",
          "Ensure that expression data and metadata are from the same source"
        )
      )
    }
    
    # Subset to common columns
    log_message(paste("Using", length(common_cols), "common columns between expression data and column metadata"), "INFO")
    expr_data <- expr_data[, common_cols, drop = FALSE]
    col_metadata <- col_metadata[common_cols, , drop = FALSE]
  }
  
  # Create GCT struct
  log_message("Creating GCT structure...")
  
  # Clean up metadata - replace NA values with "na"
  row_metadata[] <- lapply(row_metadata, function(x) {
    ifelse(is.na(x), "na", x)
  })
  
  col_metadata[] <- lapply(col_metadata, function(x) {
    ifelse(is.na(x), "na", x)
  })
  
  # Create dataset
  gct_obj <- cmapR::GCT(
    mat = as.matrix(expr_data),
    rdesc = as.data.frame(row_metadata),
    cdesc = as.data.frame(col_metadata)
  )
  
  # Write to file
  log_message(paste("Writing GCT file to", filename), "INFO")
  cmapR::write_gct(gct_obj, filename)
  
  log_message(paste("GCT file exported successfully to", filename), "SUCCESS")
  return(filename)
}

#' Create a filtered GCT file based on FDR cutoff
#'
#' @param expr_data Expression data matrix
#' @param row_metadata Row metadata with FDR values
#' @param col_metadata Column metadata
#' @param fdr_column Name of the FDR column in row_metadata
#' @param fdr_cutoff FDR cutoff value
#' @param output_prefix Output file prefix
#' @param filename Specific filename (overrides prefix if provided)
#' @return Path to the saved filtered GCT file
#' @export
export_filtered_gct <- function(expr_data, row_metadata, col_metadata, 
                               fdr_column, fdr_cutoff = 0.1,
                               output_prefix = NULL, filename = NULL) {
  log_message(paste("Creating filtered GCT file with FDR cutoff", fdr_cutoff), "STEP")
  
  # Set output filename
  if (is.null(filename)) {
    filename <- paste0(output_prefix, "_counts_filtered.gct")
  }
  
  # Check if FDR column exists
  if (!(fdr_column %in% colnames(row_metadata))) {
    log_message(paste("FDR column", fdr_column, "not found in row metadata"), "WARNING")
    return(NULL)
  }
  
  # Apply FDR filter
  fdr_values <- row_metadata[[fdr_column]]
  significant_rows <- which(!is.na(fdr_values) & fdr_values <= fdr_cutoff)
  
  if (length(significant_rows) == 0) {
    log_message(paste("No rows pass the FDR cutoff of", fdr_cutoff), "WARNING")
    return(NULL)
  }
  
  log_message(paste(length(significant_rows), "rows pass the FDR cutoff of", fdr_cutoff), "INFO")
  
  # Subset data
  filtered_expr <- expr_data[significant_rows, , drop = FALSE]
  filtered_row_metadata <- row_metadata[significant_rows, , drop = FALSE]
  
  # Export filtered GCT
  filtered_file <- export_gct_data(
    filtered_expr, 
    filtered_row_metadata, 
    col_metadata, 
    output_prefix = NULL,
    filename = filename
  )
  
  return(filtered_file)
} 