#!/usr/bin/env Rscript
#
# Data processing functions for DESeq2 LRT Step 2
#

#' Rebuild DESeq2 object with specific factor reference levels
#'
#' @param dds DESeq2 data object
#' @param factor_ref_list Named list of factors and their reference levels
#' @return Updated DESeq2 data object
#' @export
rebuild_dds <- function(dds, factor_ref_list) {
  # Convert colData to a plain data.frame so we can manipulate it
  colData_df <- as.data.frame(colData(dds))

  # Re-level factors as specified in factor_ref_list
  # factor_ref_list might be list(grouped_ref = c("grouped" = "0H_GFP_N"))
  for (f_name in names(factor_ref_list)) {
    # factor_ref_list[[f_name]] is the new reference level for factor f_name
    if (f_name %in% colnames(colData_df)) {
      ref_level <- factor_ref_list[[f_name]]
      
      # Check if the specified level exists in the factor
      if (is.factor(colData_df[[f_name]])) {
        current_levels <- levels(colData_df[[f_name]])
        if (ref_level %in% current_levels) {
          colData_df[[f_name]] <- relevel(as.factor(colData_df[[f_name]]), ref = ref_level)
          log_message(paste("Re-leveled factor", f_name, "with reference level", ref_level))
        } else {
          log_warning(paste("Reference level", ref_level, "not found in factor", f_name,
                            "- available levels:", paste(current_levels, collapse = ", ")))
        }
      } else {
        log_warning(paste("Column", f_name, "is not a factor - cannot re-level"))
      }
    } else {
      log_warning(paste("Factor", f_name, "not found in colData - cannot re-level"))
    }
  }

  # Build a new DESeqDataSet
  dds_new <- DESeqDataSetFromMatrix(
    countData = counts(dds),
    colData = colData_df,
    design = design(dds)
  )
  return(dds_new)
}

#' Harmonize sample names between data frames
#'
#' @param expression_data_df Expression data frame
#' @param read_col Column name pattern for read counts
#' @return List with read counts data frame and cleaned column names
#' @export
harmonize_sample_names <- function(expression_data_df, read_col = "TotalReads") {
  # Extract read counts columns
  read_counts_columns <- grep(
    paste(read_col, sep = ""),
    colnames(expression_data_df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(read_counts_columns) == 0) {
    log_warning(paste("No columns matching pattern", read_col, "found in expression data"))
    read_counts_columns <- grep("counts|reads", colnames(expression_data_df), 
                               value = TRUE, ignore.case = TRUE)
    
    if (length(read_counts_columns) == 0) {
      log_warning("No read count columns identified, using all non-metadata columns")
      # Guess which columns might be read counts (exclude common metadata columns)
      metadata_patterns <- c("gene", "transcript", "id", "name", "symbol", "length", 
                            "chr", "start", "end", "strand", "description")
      read_counts_columns <- colnames(expression_data_df)[
        !grepl(paste(metadata_patterns, collapse = "|"), 
              colnames(expression_data_df), ignore.case = TRUE)
      ]
    }
  }
  
  read_counts_data_df <- expression_data_df[read_counts_columns]
  log_message("Read Counts Data Extracted")
  
  # Clean and standardize column names
  original_colnames <- colnames(read_counts_data_df)
  log_message("Cleaning column names")
  
  cleaned_colnames <- sapply(original_colnames, function(s) {
    # Remove the last word if separated by space
    words <- unlist(strsplit(s, " ", fixed = TRUE))
    if (length(words) > 1) {
      paste(head(words, -1), collapse = "_")
    } else {
      s
    }
  })
  
  # Convert to lowercase and trim whitespace
  cleaned_colnames <- tolower(cleaned_colnames)
  cleaned_colnames <- trimws(cleaned_colnames)
  cleaned_colnames <- gsub("[^[:alnum:]_]", "", cleaned_colnames)
  
  # Assign cleaned column names back to the dataframe
  colnames(read_counts_data_df) <- cleaned_colnames
  
  return(list(
    read_counts_data_df = read_counts_data_df,
    read_counts_columns = read_counts_columns,
    cleaned_colnames = cleaned_colnames
  ))
}

#' Apply batch correction to normalized counts
#'
#' @param dds DESeq2 data object
#' @param col_data Metadata data frame with sample information
#' @param batch_correction_method Batch correction method ("model", "combatseq", or "none")
#' @return Normalized counts matrix with batch effects removed
#' @export
apply_batch_correction <- function(dds, col_data, batch_correction_method) {
  # Validate the batch correction method
  valid_methods <- c("none", "model", "combatseq")
  if (!batch_correction_method %in% valid_methods) {
    log_warning(paste("Invalid batch correction method:", batch_correction_method,
                     "- must be one of:", paste(valid_methods, collapse = ", ")))
    log_warning("Defaulting to 'none' for batch correction")
    batch_correction_method <- "none"
  }
  
  # Check if a batch column exists when using batch correction
  if (batch_correction_method != "none" && !"batch" %in% colnames(col_data)) {
    log_warning("Batch correction requested but no 'batch' column found in metadata")
    log_warning("Will proceed without batch correction")
    batch_correction_method <- "none"
  }
  
  if (batch_correction_method == "model" && "batch" %in% colnames(col_data)) {
    # Case: After DESeq2, apply rlog transformation and remove batch effects using limma
    log_message("Applying rlog transformation and limma batch effect removal")
    rlog_transformed <- rlog(dds, blind = FALSE)
    rlog_counts <- assay(rlog_transformed)

    # Get the design formula without batch effect
    design_formula_text <- as.character(design(dds))[2]
    # Remove the batch term if it exists in the formula
    if (grepl("batch", design_formula_text)) {
      design_formula_text <- stringr::str_remove(design_formula_text, "\\s*\\+\\s*batch")
      log_message(paste("Removed batch term from design formula:", design_formula_text))
    }
    design_formula <- as.formula(paste0("~", design_formula_text))

    # Create design matrix for removeBatchEffect
    log_message("Creating design matrix for batch effect removal")
    design_matrix <- model.matrix(design_formula, data = col_data)

    # Apply removeBatchEffect on rlog normalized counts
    log_message("Removing batch effects using limma")
    norm_counts <- limma::removeBatchEffect(rlog_counts, 
                                         batch = col_data$batch, 
                                         design = design_matrix)
    
    return(norm_counts)
  } else if (batch_correction_method == "combatseq" && "batch" %in% colnames(col_data)) {
    # Case: Apply ComBat-Seq for batch correction
    log_message("Applying ComBat-Seq batch correction on raw counts")
    
    # ComBat-Seq works on raw counts
    raw_counts <- counts(dds)
    
    # Create a model matrix for biological variables to preserve
    design_formula_text <- as.character(design(dds))[2]
    if (grepl("batch", design_formula_text)) {
      design_formula_text <- stringr::str_remove(design_formula_text, "\\s*\\+\\s*batch")
    }
    design_formula <- as.formula(paste0("~", design_formula_text))
    
    log_message("Creating model matrix for biological variables")
    bio_model <- model.matrix(design_formula, data = col_data)
    
    # Apply ComBat-Seq
    log_message("Running ComBat-Seq")
    corrected_counts <- sva::ComBat_seq(
      counts = raw_counts,
      batch = col_data$batch,
      group = NULL,
      covar_mod = bio_model
    )
    
    # Apply rlog transformation to the batch-corrected counts
    log_message("Applying rlog transformation to batch-corrected counts")
    # Create a temporary DESeqDataSet with the corrected counts
    dds_corrected <- DESeqDataSetFromMatrix(
      countData = corrected_counts,
      colData = col_data,
      design = design(dds)
    )
    
    norm_counts <- assay(rlog(dds_corrected, blind = FALSE))
    return(norm_counts)
  } else {
    # Default case: No batch correction, just apply rlog
    log_message("Applying rlog transformation without batch correction")
    norm_counts <- assay(rlog(dds, blind = FALSE))
    return(norm_counts)
  }
}

#' Prepare final data frame with ordered columns
#'
#' @param annotated_expression_combined_df Combined expression data with annotations
#' @param dds DESeq2 data object
#' @return Data frame with reordered columns
#' @export
prepare_final_data <- function(annotated_expression_combined_df, dds) {
  # Remove duplicate columns
  annotated_expression_combined_df <- annotated_expression_combined_df[, !duplicated(names(annotated_expression_combined_df))]
  
  # Specify the desired order of some key columns
  annotation_cols <- c("GeneId", "gene_id", "RefseqId", "gene_name", "Chrom", "TxStart", "TxEnd", "Strand", "baseMean")
  
  # Keep only annotation columns that actually exist
  annotation_cols <- annotation_cols[annotation_cols %in% colnames(annotated_expression_combined_df)]
  
  if (length(annotation_cols) == 0) {
    log_warning("No standard annotation columns found in data")
    # Try to identify any potential ID columns
    id_cols <- grep("id$|name$|symbol$", colnames(annotated_expression_combined_df), 
                   value = TRUE, ignore.case = TRUE)
    if (length(id_cols) > 0) {
      annotation_cols <- id_cols
      log_message(paste("Using identified ID columns:", paste(id_cols, collapse = ", ")))
    }
  }
  
  # Get the current sample order from the DESeq object
  sample_order <- colnames(dds)
  
  # Extract contrast columns: any column name that contains "vs" or standard result column suffixes
  contrast_cols <- grep("_vs_|_LFC$|_pvalue$|_FDR$|_padj$|log2FoldChange|pvalue|padj", 
                      colnames(annotated_expression_combined_df), value = TRUE)
  
  # Create the final order:
  # 1. Annotation columns (in the specified order, if they exist)
  # 2. The sample columns (in the exact order they appear in the DESeq object)
  # 3. All contrast-related columns
  final_cols <- c(
    intersect(annotation_cols, colnames(annotated_expression_combined_df)),
    intersect(sample_order, colnames(annotated_expression_combined_df)),
    intersect(contrast_cols, colnames(annotated_expression_combined_df))
  )
  
  # Add any remaining columns that weren't categorized
  remaining_cols <- base::setdiff(colnames(annotated_expression_combined_df), final_cols)
  final_cols <- c(final_cols, remaining_cols)
  
  # Reorder the columns of the data frame
  annotated_expression_combined_df <- annotated_expression_combined_df[, final_cols, drop = FALSE]
  
  return(annotated_expression_combined_df)
}

#' Verify consistency between Step 1 and Step 2 analysis
#'
#' @param col_data Column data from DESeq2 dataset
#' @param design_formula Design formula from DESeq2 dataset
#' @param args Command-line arguments
#' @return TRUE if consistent, otherwise stops with an error
#' @export
verify_step_consistency <- function(col_data, design_formula, args) {
  log_message("Verifying consistency between Step 1 and Step 2", "INFO")
  
  # Extract design formula as text for comparison
  design_text <- as.character(design_formula)[2]
  
  # Check if design_text contains the required LRT factor
  if (!is.null(args$lrt_factor)) {
    if (!grepl(args$lrt_factor, design_text, fixed = TRUE)) {
      stop(paste("LRT factor", args$lrt_factor, "not found in design formula from Step 1"))
    }
    log_message(paste("Verified LRT factor", args$lrt_factor, "is in design formula"))
  }
  
  # Check batch correction consistency
  if (args$batchcorrection != "none") {
    # Check if batch column exists in the column data
    if (!"batch" %in% colnames(col_data)) {
      log_warning("Batch correction requested but 'batch' column not found in DESeq2 data")
      log_warning("Batch correction will be disabled")
      args$batchcorrection <- "none"
    } else {
      log_message("Verified batch column exists for batch correction")
    }
    
    # For model-based batch correction, check if batch is in the design formula
    if (args$batchcorrection == "model" && !grepl("batch", design_text, fixed = TRUE)) {
      log_warning("Model-based batch correction requested but 'batch' not in design formula")
      log_warning("Setting batch correction to 'limma' instead of 'model'")
      args$batchcorrection <- "limma"
    }
  }
  
  # Check if the factors referenced in contrast_df exist in the column data
  if (!is.null(args$contrast_df) && file.exists(args$contrast_df)) {
    contrast_df <- read.delim(args$contrast_df, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # Extract factor names from contrast data
    factor_names <- unique(c(
      contrast_df$factor,
      contrast_df$specificity_group
    ))
    
    # Remove NA or empty values
    factor_names <- factor_names[!is.na(factor_names) & factor_names != ""]
    
    # Check if each factor exists in column data
    missing_factors <- factor_names[!factor_names %in% colnames(col_data)]
    
    if (length(missing_factors) > 0) {
      log_warning(paste("The following factors from contrast file are not in column data:", 
                       paste(missing_factors, collapse=", ")))
    }
  }
  
  # All checks passed
  log_message("Step 1 and Step 2 consistency verified", "SUCCESS")
  return(TRUE)
} 