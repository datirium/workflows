#!/usr/bin/env Rscript
#
# Export functions for DESeq2 workflows
# Focuses on specialized export functionality not covered by output_utils.R
#

#' Export normalized count data in formats suitable for downstream analysis
#'
#' @param dds DESeq2 object with normalized counts
#' @param output_prefix Output prefix for files
#' @param threshold Optional RPKM threshold for filtering
#' @return List of paths to exported files
#' @export
export_normalized_counts <- function(dds, output_prefix, threshold = NULL) {
  # Ensure we have the DESeq2 package
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package is required for this function")
  }
  
  # Get normalized counts
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  # File paths for all counts
  all_counts_path <- get_output_filename(output_prefix, "counts_all", "gct")
  
  # Write all counts to GCT
  write_gct_file(norm_counts, all_counts_path)
  
  result_paths <- list(all_counts = all_counts_path)
  
  # Apply threshold filtering if specified
  if (!is.null(threshold)) {
    # Calculate row means
    row_means <- rowMeans(norm_counts)
    
    # Apply threshold
    filtered_counts <- norm_counts[row_means >= threshold, ]
    
    # File path for filtered counts
    filtered_counts_path <- get_output_filename(output_prefix, "counts_filtered", "gct")
    
    # Write filtered counts to GCT
    write_gct_file(filtered_counts, filtered_counts_path)
    
    result_paths$filtered_counts <- filtered_counts_path
    message(paste("Filtered counts from", nrow(norm_counts), "to", nrow(filtered_counts), "genes"))
  }
  
  return(result_paths)
}

#' Export DESeq2 results to various formats
#'
#' @param deseq_results DESeq2 results object
#' @param output_prefix Output prefix for files
#' @param gene_info Optional gene annotation data frame
#' @param output_name Type of output (e.g., "report", "gene_exp_table")
#' @return Path to exported results file
#' @export
export_deseq_results <- function(deseq_results, output_prefix, 
                                gene_info = NULL, output_name = "report") {
  # Convert to data frame if needed
  if (!is.data.frame(deseq_results)) {
    results_df <- as.data.frame(deseq_results)
  } else {
    results_df <- deseq_results
  }
  
  # Add gene info if provided
  if (!is.null(gene_info) && is.data.frame(gene_info)) {
    # Check if rownames match
    common_genes <- intersect(rownames(results_df), rownames(gene_info))
    
    if (length(common_genes) > 0) {
      # Merge based on common genes
      results_df <- cbind(
        results_df[common_genes, , drop = FALSE],
        gene_info[common_genes, , drop = FALSE]
      )
      
      message(paste("Added gene annotations for", length(common_genes), "genes"))
    } else {
      message("No matching genes found between results and gene info")
    }
  }
  
  # Create output filename
  if (output_name == "report") {
    output_file <- get_output_filename(output_prefix, "report", "tsv")
  } else if (output_name == "gene_exp_table") {
    output_file <- get_output_filename(output_prefix, "gene_exp_table", "tsv")
  } else if (output_name == "contrasts_table") {
    output_file <- get_output_filename(output_prefix, "contrasts_table", "tsv")
  } else {
    output_file <- get_output_filename(output_prefix, output_name, "tsv")
  }
  
  # Export to TSV
  write_deseq_results(results_df, output_file)
  
  return(output_file)
}

#' Export DESeq2 results to files specifically for GSEA analysis
#'
#' @param deseq_results DESeq2 results object
#' @param norm_counts Normalized count matrix
#' @param metadata Sample metadata
#' @param condition_col Column in metadata containing condition information
#' @param output_prefix Output prefix for files
#' @return List of paths to exported files
#' @export
export_for_gsea <- function(deseq_results, norm_counts, metadata, condition_col, output_prefix) {
  result_paths <- list()
  
  # 1. Export ranked gene list with metric
  ranked_results <- as.data.frame(deseq_results)
  
  # Calculate ranking metric (signed -log10 p-value)
  if ("log2FoldChange" %in% colnames(ranked_results) && "pvalue" %in% colnames(ranked_results)) {
    ranked_results$metric <- sign(ranked_results$log2FoldChange) * -log10(ranked_results$pvalue)
    
    # Sort by metric for GSEA
    ranked_results <- ranked_results[order(ranked_results$metric, decreasing = TRUE), ]
    
    # Export to TSV
    ranked_file <- get_output_filename(output_prefix, "ranked_genes", "rnk")
    write.table(
      data.frame(
        gene = rownames(ranked_results),
        metric = ranked_results$metric
      ),
      file = ranked_file,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    
    result_paths$ranked_file <- ranked_file
    message(paste("Exported ranked gene list for GSEA to", ranked_file))
  }
  
  # 2. Export expression data in GCT format
  gct_file <- get_output_filename(output_prefix, "expression", "gct")
  write_gct_file(norm_counts, gct_file)
  result_paths$expression_file <- gct_file
  
  # 3. Export phenotype data in CLS format
  if (!is.null(metadata) && condition_col %in% colnames(metadata)) {
    # Get conditions for samples in the same order as columns in norm_counts
    sample_conditions <- metadata[colnames(norm_counts), condition_col]
    
    cls_file <- get_output_filename(output_prefix, "phenotypes", "cls")
    write_cls_file(sample_conditions, cls_file)
    result_paths$phenotype_file <- cls_file
    message(paste("Exported phenotype classes for GSEA to", cls_file))
  }
  
  return(result_paths)
}

#' Export DESeq2 object to RDS file
#'
#' @param dds DESeq2 object
#' @param output_prefix Output prefix
#' @param type Type of DESeq2 object (e.g., "dds", "contrasts")
#' @return Path to exported RDS file
#' @export
export_deseq_object <- function(dds, output_prefix, type = "contrasts") {
  output_file <- get_output_filename(output_prefix, type, "rds")
  saveRDS(dds, output_file)
  message(paste("Saved DESeq2 object to", output_file))
  return(output_file)
}

#' Create a standard set of visualization outputs
#'
#' @param deseq_object DESeq2 object
#' @param deseq_results DESeq2 results object
#' @param output_prefix Output prefix for files
#' @param transformed_counts Optional pre-transformed counts (VST or rlog)
#' @param metadata Sample metadata
#' @param color_by Column in metadata to use for coloring
#' @param top_n_genes Number of top genes to include in heatmap
#' @return List of paths to exported visualization files
#' @export
export_visualizations <- function(deseq_object, deseq_results, output_prefix,
                                transformed_counts = NULL, metadata = NULL,
                                color_by = NULL, top_n_genes = 50) {
  result_paths <- list()
  
  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for visualizations")
    return(result_paths)
  }
  
  # 1. MA Plot (if we have log2FoldChange)
  if ("log2FoldChange" %in% colnames(deseq_results) && 
      "baseMean" %in% colnames(deseq_results)) {
      
    # Convert results to data frame if needed
    if (!is.data.frame(deseq_results)) {
      results_df <- as.data.frame(deseq_results)
    } else {
      results_df <- deseq_results
    }
    
    # Create MA plot
    ma_plot <- ggplot2::ggplot(
      results_df, 
      ggplot2::aes(x = log10(baseMean), y = log2FoldChange, color = padj < 0.1)
    ) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::scale_color_manual(values = c("grey", "red"), name = "Significant") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::xlab("log10(Mean Expression)") +
      ggplot2::ylab("log2(Fold Change)") +
      ggplot2::ggtitle("MA Plot: Fold Change vs Mean Expression")
    
    # Save plot
    ma_files <- save_plot(ma_plot, output_prefix, "ma_plot")
    result_paths$ma_plot <- ma_files
  }
  
  # 2. MDS Plot as HTML interactive visualization
  # Get transformed data if not provided
  if (is.null(transformed_counts) && !is.null(deseq_object)) {
    # Try to perform VST
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      transformed_counts <- tryCatch({
        DESeq2::vst(deseq_object, blind = TRUE)
      }, error = function(e) {
        message("Could not perform VST: ", e$message)
        NULL
      })
    }
  }
  
  if (!is.null(transformed_counts)) {
    # Create interactive MDS plot
    mds_file <- get_output_filename(output_prefix, "mds_plot", "html")
    mds_path <- generate_mds_plot_html(
      assay(transformed_counts), 
      mds_file,
      metadata = metadata,
      color_by = color_by,
      title = "MDS Plot of Samples"
    )
    
    if (!is.null(mds_path)) {
      result_paths$mds_plot <- mds_path
    }
    
    # 3. Heatmap of top genes
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      # Get top differentially expressed genes
      if (!is.null(deseq_results) && "padj" %in% colnames(deseq_results)) {
        top_genes <- head(rownames(deseq_results)[order(deseq_results$padj)], top_n_genes)
        
        if (length(top_genes) > 0) {
          # Extract data for these genes
          heatmap_data <- assay(transformed_counts)[top_genes, ]
          
          # Create annotation if metadata is available
          if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
            # Create annotation data frame
            anno_df <- data.frame(
              Condition = metadata[colnames(heatmap_data), color_by],
              row.names = colnames(heatmap_data)
            )
            
            # Create heatmap with annotation
            heatmap_plot <- pheatmap::pheatmap(
              heatmap_data,
              annotation_col = anno_df,
              scale = "row",
              clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              show_rownames = (nrow(heatmap_data) <= 50),
              main = paste("Top", length(top_genes), "Differentially Expressed Genes"),
              silent = TRUE
            )
          } else {
            # Create heatmap without annotation
            heatmap_plot <- pheatmap::pheatmap(
              heatmap_data,
              scale = "row",
              clustering_distance_rows = "correlation",
              clustering_distance_cols = "correlation",
              show_rownames = (nrow(heatmap_data) <= 50),
              main = paste("Top", length(top_genes), "Differentially Expressed Genes"),
              silent = TRUE
            )
          }
          
          # Save heatmap
          heatmap_files <- save_plot(heatmap_plot, output_prefix, "expression_heatmap")
          result_paths$heatmap <- heatmap_files
        }
      }
    }
  }
  
  return(result_paths)
} 