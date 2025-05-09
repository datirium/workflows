#!/usr/bin/env Rscript
#
# Data processing functions for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for loading, preprocessing, and filtering
# expression data for the DESeq analysis workflow.
#
# Version: 0.1.0

# Define constants for column names
READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")
RPKM_UNTREATED_ALIAS <- "RpkmCondition1"
RPKM_TREATED_ALIAS <- "RpkmCondition2"

#' Determine file type from filename extension
#'
#' @param filename The path to the file
#' @return The appropriate separator character for the file type
get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
}

#' Filter expression data by RPKM threshold
#'
#' @param expression_df Data frame containing expression data
#' @param n RPKM threshold value
#' @return Filtered data frame
filter_rpkm <- function(expression_df, n) {
  log_message(paste("Filtering expression data with RPKM cutoff:", n))
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for RPKM filtering")
  }
  
  filtered_df <- expression_df %>%
    dplyr::filter(dplyr::if_any(dplyr::contains("Rpkm"), ~. > n))
  
  log_message(paste("Filtered", nrow(expression_df) - nrow(filtered_df), "rows by RPKM cutoff"))
  
  return(filtered_df)
}

#' Load expression data from isoform files and merge them
#'
#' @param input_files Vector of file paths
#' @param sample_aliases Vector of sample aliases
#' @param read_column Column name for read counts
#' @param rpkm_column Column name for RPKM values
#' @param rpkm_alias Alias for combined RPKM values
#' @param condition_name Condition name (e.g., "treated" or "untreated")
#' @param intersect_by Vector of column names to use for merging
#' @param digits Number of digits for rounding
#' @param batch_data Batch metadata for multi-factor analysis
#' @param prev_data Previously collected data (optional)
#' @return List containing merged expression data and metadata
load_isoform_set <- function(input_files, sample_aliases, read_column, rpkm_column, rpkm_alias, condition_name, intersect_by, digits = NULL, batch_data = NULL, prev_data = NULL) {
  log_message(paste("Loading isoform set for condition:", condition_name))
  
  # Check and validate inputs
  if (is.null(input_files) || length(input_files) == 0) {
    stop("No input files provided")
  }
  
  # Validate sample aliases
  if (is.null(sample_aliases) || length(sample_aliases) == 0) {
    log_warning("No sample aliases provided. Using input file basenames")
    sample_aliases <- basename(input_files)
    # Remove extensions from filenames
    sample_aliases <- gsub("\\.(csv|tsv)$", "", sample_aliases, ignore.case = TRUE)
  }
  
  # Check for mismatched file and alias counts
  if (length(input_files) != length(sample_aliases)) {
    log_warning(paste("Mismatch between input files count (", length(input_files), ") and sample aliases count (", length(sample_aliases), ")"))
    
    # Handle two specific cases:
    # 1. More files than aliases: generate aliases from file basenames
    # 2. More aliases than files: truncate aliases to match files
    
    if (length(input_files) > length(sample_aliases)) {
      log_message("More input files than sample aliases. Using file basenames for missing aliases.")
      
      # Generate aliases from file paths for the missing slots
      missing_aliases_count <- length(input_files) - length(sample_aliases)
      file_basenames <- basename(input_files[(length(sample_aliases)+1):length(input_files)])
      file_basenames <- gsub("\\.(csv|tsv)$", "", file_basenames, ignore.case = TRUE)
      
      # Append the generated aliases
      sample_aliases <- c(sample_aliases, file_basenames)
      log_message(sprintf("Added %d aliases from file basenames. Now have %d aliases.", 
                      missing_aliases_count, length(sample_aliases)))
    } else if (length(sample_aliases) > length(input_files)) {
      log_message("More sample aliases than input files. Truncating aliases list.")
      sample_aliases <- sample_aliases[1:length(input_files)]
    }
  }
  
  # Clean sample aliases - replace spaces, quotes, etc. with underscores
  sample_aliases <- gsub("'|\"| ", "_", sample_aliases)
  
  # Initialize the result container
  result <- list()
  result$collected_isoforms <- NULL
  
  for (i in 1:length(input_files)) {
    # Log processing
    log_message(paste("Loading file", i, "of", length(input_files), ":", input_files[i]))
    
    # Get file type based on extension
    file_type <- get_file_type(input_files[i])
    
    # Read file with appropriate delimiter
    current_data <- tryCatch({
      read.table(
        input_files[i],
        sep = file_type,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      log_error(paste("Error reading file:", input_files[i], "-", e$message))
      return(NULL)
    })
    
    # Skip if file couldn't be read
    if (is.null(current_data)) {
      log_warning(paste("Skipping file due to read error:", input_files[i]))
      next
    }
    
    # Rename columns
    new_read_colname <- paste(sample_aliases[i], condition_name, sep = "_")
    colnames(current_data)[colnames(current_data) == read_column] <- new_read_colname
    colnames(current_data)[colnames(current_data) == rpkm_column] <- paste(condition_name, i, rpkm_column, sep = " ")
    
    # Create column data
    if (!is.null(batch_data)) {
      batch <- batch_data[sample_aliases[i], "batch"]
      log_message(
        paste0(
          "Loaded ", nrow(current_data),
          " rows from '", input_files[i],
          "' as '", new_read_colname,
          "', batch '", batch, "'"
        )
      )
      column_data_frame <- data.frame(condition_name, batch, row.names = c(new_read_colname))
    } else {
      log_message(
        paste0(
          "Loaded ", nrow(current_data),
          " rows from '", input_files[i],
          "' as '", new_read_colname, "'"
        )
      )
      column_data_frame <- data.frame(condition_name, row.names = c(new_read_colname))
    }
    
    # Initialize or update collected data
    if (is.null(prev_data)) {
      result$collected_isoforms <- current_data
      result$read_colnames <- c(new_read_colname)
      result$column_data <- column_data_frame
      result$rpkm_colnames <- character(0)
    } else {
      result$collected_isoforms <- merge(
        prev_data$collected_isoforms,
        current_data,
        by = intersect_by,
        sort = FALSE
      )
      result$read_colnames <- c(prev_data$read_colnames, new_read_colname)
      result$column_data <- rbind(prev_data$column_data, column_data_frame)
    }
  }
  
  # Calculate average RPKM for the condition
  rpkm_columns <- grep(
    paste("^", condition_name, " [0-9]+ ", rpkm_column, sep = ""),
    colnames(result$collected_isoforms),
    value = TRUE,
    ignore.case = TRUE
  )
  
  result$collected_isoforms[rpkm_alias] <- format(
    rowSums(result$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(input_files),
    digits = digits
  )
  
  # Update RPKM column names and remove individual RPKM columns
  result$rpkm_colnames <- c(result$rpkm_colnames, rpkm_alias)
  result$collected_isoforms <- result$collected_isoforms[, !colnames(result$collected_isoforms) %in% rpkm_columns]
  
  return(result)
}

#' Generate an analysis summary markdown file
#'
#' @param batchcorrection The batch correction method used
#' @param batchfile The batch file provided (if any)
#' @param deseq_results Results from DESeq/DESeq2 analysis
#' @param output_file Path to the output markdown file
#' @return None
generate_md <- function(batchcorrection, batchfile, deseq_results, output_file) {
  # Initialize the markdown content
  md_lines <- c("# DESeq2 Summary\n")

  # Add warning message if applicable
  if (!is.null(batchcorrection) && (batchcorrection == "combatseq" || batchcorrection == "model") && is.null(batchfile)) {
    warning_lines <- c(
      "# Warning!\n",
      "---\n",
      "**You provided a batch-correction method, but not a batch-file.**\n",
      "The chosen parameter was ignored.\n",
      "Please ensure that you provide a batch file when using the following batch correction methods:\n",
      "- **combatseq**", 
      "- **model**\n",
      "If you do not need batch correction, set the method to 'none'.\n",
      "---\n"
    )
    md_lines <- c(md_lines, warning_lines)
  }

  # Add DESeq results summary if provided
  if (!is.null(deseq_results)) {
    # Extract total genes with non-zero read count
    total_genes <- sum(!is.na(deseq_results$baseMean) & deseq_results$baseMean > 0)

    # Get thresholds from deseq_results metadata
    metadata_res <- metadata(deseq_results)
    lfc_threshold <- ifelse(!is.null(metadata_res$lfcThreshold), metadata_res$lfcThreshold, 0)
    p_adj_threshold <- ifelse(!is.null(metadata_res$alpha), metadata_res$alpha, 0.1)

    # Extract counts
    res_nonzero <- deseq_results[!is.na(deseq_results$baseMean) & deseq_results$baseMean > 0,]

    # Calculate stats
    lfc_up <- sum(
      res_nonzero$log2FoldChange > lfc_threshold &
        res_nonzero$padj < p_adj_threshold,
      na.rm = TRUE
    )
    lfc_down <- sum(
      res_nonzero$log2FoldChange < -lfc_threshold &
        res_nonzero$padj < p_adj_threshold,
      na.rm = TRUE
    )
    outliers <- sum(
      is.na(res_nonzero$pvalue) & !is.na(res_nonzero$baseMean),
      na.rm = TRUE
    )
    low_counts <- sum(
      is.na(res_nonzero$padj) & is.na(res_nonzero$pvalue),
      na.rm = TRUE
    )
    mean_count <- if (!is.null(metadata_res$filterThreshold)) {
      round(metadata_res$filterThreshold, 2)
    } else {
      "-"
    }

    # Add results summary
    md_lines <- c(md_lines, 
      "\n## Analysis Results\n",
      paste("* Total genes with non-zero counts:", total_genes),
      paste("* Significant genes (FDR <", p_adj_threshold, "):", lfc_up + lfc_down),
      paste("* Up-regulated genes (LFC >", lfc_threshold, "):", lfc_up, 
            sprintf("(%.1f%%)", ifelse(total_genes > 0, (lfc_up / total_genes) * 100, 0))),
      paste("* Down-regulated genes (LFC < -", lfc_threshold, "):", lfc_down,
            sprintf("(%.1f%%)", ifelse(total_genes > 0, (lfc_down / total_genes) * 100, 0))),
      paste("* Outliers with high Cook's distance:", outliers),
      paste("* Genes filtered due to low counts:", low_counts),
      if (mean_count != "-") {
        paste("* Mean count threshold:", mean_count)
      }
    )
    
    # Add footnotes
    md_lines <- c(md_lines,
      "\n## Notes\n",
      "1. Outliers are genes with high Cook's distance (see [cooksCutoff](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier)).",
      "2. Low counts are genes filtered out due to low mean counts (see [independentFiltering](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#independent-filtering-of-results))."
    )
  }

  # Write the content to the output file using our consolidated function
  write_markdown_summary(md_lines, output_file)
} 