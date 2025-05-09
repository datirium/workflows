#!/usr/bin/env Rscript

# --- Data processing functions ---

# Constants for column name patterns
READ_COL <- "Read"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- "GeneId"

# Function to check the delimiter used in a file
check_file_delimiter <- function(file_path) {
  # Default to tab delimiter
  delimiter <- "\t"
  
  # Check file extension
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    delimiter <- ","
  } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    delimiter <- "\t"
  } else {
    # Try to detect delimiter from file content
    first_line <- readLines(file_path, n = 1)
    
    if (grepl(",", first_line)) {
      delimiter <- ","
    } else if (grepl("\t", first_line)) {
      delimiter <- "\t"
    } else if (grepl(";", first_line)) {
      delimiter <- ";"
    }
  }
  
  message(paste("Detected delimiter for file", file_path, ":", delimiter))
  return(delimiter)
}

# Function to load expression data from multiple files
load_expression_data <- function(args) {
  message("Loading expression data from files...")
  
  # Validate input files exist
  for (file in args$input) {
    if (!file.exists(file)) {
      stop(paste("Input file does not exist:", file))
    }
  }
  
  # Determine file formats and load data
  expr_data_list <- list()
  for (i in seq_along(args$input)) {
    file_path <- args$input[i]
    file_name <- if (!is.null(args$name) && length(args$name) >= i) args$name[i] else basename(file_path)
    
    message(paste("Loading file", i, "of", length(args$input), ":", file_path))
    
    # Determine file format
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      delimiter <- ","
    } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
      delimiter <- "\t"
    } else {
      # Default to CSV
      message(paste("Unknown file extension for", file_path, "- assuming CSV format"))
      delimiter <- ","
    }
    
    # Load file with better error handling
    tryCatch({
      data <- read.table(file_path, sep = delimiter, header = TRUE, stringsAsFactors = FALSE)
      message(paste("  Loaded", nrow(data), "rows and", ncol(data), "columns"))
      
      # Validate required columns
      required_cols <- c("GeneId", "TotalReads", "Rpkm")
      missing_cols <- required_cols[!required_cols %in% colnames(data)]
      
      if (length(missing_cols) > 0) {
        warning(paste("File", file_path, "is missing required columns:", paste(missing_cols, collapse = ", ")))
        
        # Try to fix common column name variations
        col_mappings <- list(
          "GeneId" = c("gene_id", "geneid", "gene", "Gene_ID", "Gene"),
          "TotalReads" = c("total_reads", "reads", "count", "counts", "read_count"),
          "Rpkm" = c("rpkm", "FPKM", "fpkm")
        )
        
        for (req_col in required_cols) {
          if (!req_col %in% colnames(data) && any(col_mappings[[req_col]] %in% colnames(data))) {
            # Find the first matching alternative column name
            alt_col <- col_mappings[[req_col]][col_mappings[[req_col]] %in% colnames(data)][1]
            message(paste("  Using column", alt_col, "as", req_col))
            colnames(data)[colnames(data) == alt_col] <- req_col
          }
        }
        
        # Check again after attempted fixes
        missing_cols <- required_cols[!required_cols %in% colnames(data)]
        if (length(missing_cols) > 0) {
          stop(paste("File", file_path, "is still missing required columns after fixes:", paste(missing_cols, collapse = ", ")))
        }
      }
      
      # Add sample name to column names for clarity
      read_col <- paste(file_name, "TotalReads", sep = " ")
      rpkm_col <- paste(file_name, "Rpkm", sep = " ")
      
      colnames(data)[colnames(data) == "TotalReads"] <- read_col
      colnames(data)[colnames(data) == "Rpkm"] <- rpkm_col
      
      expr_data_list[[i]] <- data
      
    }, error = function(e) {
      stop(paste("Error loading file", file_path, ":", e$message))
    })
  }
  
  # Merge all data frames
  message("Merging data from all input files...")
  if (length(expr_data_list) == 0) {
    stop("No valid expression data loaded")
  } else if (length(expr_data_list) == 1) {
    merged_data <- expr_data_list[[1]]
  } else {
    # Use reduce with merge to join all data frames by GeneId
    merged_data <- Reduce(function(x, y) merge(x, y, by = "GeneId", all = TRUE), expr_data_list)
  }
  
  # Check for NA values and handle them
  na_count <- sum(is.na(merged_data))
  if (na_count > 0) {
    message(paste("Found", na_count, "NA values in merged data - replacing with zeros"))
    merged_data[is.na(merged_data)] <- 0
  }
  
  message(paste("Final merged data has", nrow(merged_data), "rows and", ncol(merged_data), "columns"))
  return(merged_data)
}

# Function to filter expression data based on RPKM values
filter_rpkm <- function(expression_data, rpkm_cutoff) {
  if (is.null(rpkm_cutoff) || rpkm_cutoff <= 0) {
    message("No RPKM filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Get RPKM columns
  rpkm_cols <- grep(RPKM_COL, colnames(expression_data), ignore.case = TRUE, value = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No RPKM columns found, skipping RPKM filtering")
    return(expression_data)
  }
  
  # Keep rows where any RPKM value exceeds cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  filtered_data <- expression_data[keep_rows, , drop = FALSE]
  
  message(paste("RPKM filtering: removed", sum(!keep_rows), "rows, kept", sum(keep_rows), "rows"))
  
  return(filtered_data)
}

# Function to process count data and apply batch correction
process_count_data <- function(args, counts, metadata, design_formula) {
  message("Processing count data...")
  
  # Set up variables for batch correction result
  batch_warning <- NULL
  countData <- counts
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    message(paste("Applying batch correction method:", args$batchcorrection))
    
    if (!"batch" %in% colnames(metadata)) {
      batch_warning <- "Batch correction requested but no 'batch' column found in metadata. Batch correction skipped."
      warning(batch_warning)
    } else {
      # Check batch correction method
      if (args$batchcorrection == "combatseq") {
        # Apply ComBat-Seq batch correction
        countData <- apply_combatseq_correction(counts, metadata, design_formula)
      } else if (args$batchcorrection == "model") {
        # Include batch in model (no data modification needed)
        design_formula <- update(design_formula, ~ . + batch)
        message("Added batch to design formula: ", deparse(design_formula))
      } else {
        batch_warning <- paste("Unknown batch correction method:", args$batchcorrection, ". Batch correction skipped.")
        warning(batch_warning)
      }
    }
  }
  
  # Return processed data
  return(list(
    countData = countData,
    design_formula = design_formula,
    batch_warning = batch_warning
  ))
}

# Helper function for ComBat-Seq batch correction
apply_combatseq_correction <- function(counts, metadata, design_formula) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    warning("Package 'sva' is required for ComBat-Seq batch correction but not available. Skipping batch correction.")
    return(counts)
  }
  
  message("Applying ComBat-Seq batch correction...")
  
  # Extract batch information
  batch <- metadata$batch
  
  # Get model matrix from design formula
  mod <- model.matrix(design_formula, data = metadata)
  
  # Apply ComBat-Seq
  corrected_counts <- sva::ComBat_seq(
    counts = as.matrix(counts),
    batch = batch,
    group = mod
  )
  
  # Convert back to data frame
  corrected_counts <- as.data.frame(corrected_counts)
  
  # Ensure row and column names are preserved
  rownames(corrected_counts) <- rownames(counts)
  colnames(corrected_counts) <- colnames(counts)
  
  message("Batch correction applied successfully")
  
  return(corrected_counts)
}

# Function to scale expression data based on scaling_type
scale_expression_data <- function(expression_data, scaling_type) {
  # Ensure scaling_type is a character
  scaling_type <- as.character(scaling_type)
  
  if (scaling_type == "none" || is.null(scaling_type)) {
    message("No scaling applied to expression data")
    return(expression_data)
  }
  
  message(paste("Applying", scaling_type, "scaling to expression data"))
  
  # Apply scaling
  if (scaling_type == "minmax") {
    # Min-max scaling to range [-2, 2]
    scaled_data <- apply(expression_data, 1, function(x) {
      min_val <- min(x, na.rm = TRUE)
      max_val <- max(x, na.rm = TRUE)
      if (max_val == min_val) {
        return(rep(0, length(x)))
      }
      return(4 * (x - min_val) / (max_val - min_val) - 2)
    })
    # Transpose back to original orientation
    scaled_data <- t(scaled_data)
    
  } else if (scaling_type == "zscore") {
    # Z-score scaling (mean=0, sd=1)
    scaled_data <- t(apply(expression_data, 1, function(x) {
      mean_val <- mean(x, na.rm = TRUE)
      sd_val <- sd(x, na.rm = TRUE)
      if (sd_val == 0) {
        return(rep(0, length(x)))
      }
      return((x - mean_val) / sd_val)
    }))
  } else {
    stop(paste("Unsupported scaling type:", scaling_type))
  }
  
  return(scaled_data)
}

# Function to perform clustering based on cluster method
perform_clustering <- function(expression_data, cluster_method, row_distance, column_distance, k, kmax) {
  if (cluster_method == "none" || is.null(cluster_method)) {
    message("No clustering performed")
    return(list(expression_data = expression_data, row_order = 1:nrow(expression_data), col_order = 1:ncol(expression_data)))
  }
  
  message(paste("Performing", cluster_method, "clustering"))
  
  # Load required library
  require(hopach)
  
  # Initialize clustering results
  row_clustering <- NULL
  col_clustering <- NULL
  
  # Check what parameters are accepted by hopach in this version
  hopach_args <- names(formals(hopach::hopach))
  message("Available hopach parameters: ", paste(hopach_args, collapse=", "))
  
  # Perform row clustering if requested
  if (cluster_method %in% c("row", "both")) {
    message(paste("Row clustering with distance:", row_distance))
    
    tryCatch({
      # Try modern parameter names first
      if ("distmethod" %in% hopach_args) {
        row_clustering <- hopach::hopach(
          data = expression_data,
          diss = NULL,
          distmethod = as.character(row_distance),
          K = as.numeric(k),
          kmax = as.numeric(kmax)
        )
      } 
      # Fall back to old parameter names
      else if ("dist.method" %in% hopach_args) {
        row_clustering <- hopach::hopach(
          data = expression_data,
          dmat = NULL,
          dist.method = as.character(row_distance),
          k = as.numeric(k),
          kmax = as.numeric(kmax),
          khigh = NULL
        )
      }
      # If neither works, use minimal parameters
      else {
        row_clustering <- hopach::hopach(
          data = expression_data
        )
      }
    }, error = function(e) {
      # Attempt with default parameters if specific ones fail
      message(paste("Error in row clustering:", e$message))
      message("Attempting row clustering with default parameters")
      tryCatch({
        row_clustering <<- hopach::hopach(data = expression_data)
      }, error = function(e2) {
        message(paste("Error in basic row clustering:", e2$message))
        message("Skipping row clustering")
      })
    })
  }
  
  # Perform column clustering if requested
  if (cluster_method %in% c("column", "both")) {
    message(paste("Column clustering with distance:", column_distance))
    
    tryCatch({
      # Try modern parameter names first
      if ("distmethod" %in% hopach_args) {
        col_clustering <- hopach::hopach(
          data = t(expression_data),  # Transpose for column clustering
          diss = NULL,
          distmethod = as.character(column_distance),
          K = as.numeric(k),
          kmax = as.numeric(kmax)
        )
      } 
      # Fall back to old parameter names
      else if ("dist.method" %in% hopach_args) {
        col_clustering <- hopach::hopach(
          data = t(expression_data),  # Transpose for column clustering
          dmat = NULL,
          dist.method = as.character(column_distance),
          k = as.numeric(k),
          kmax = as.numeric(kmax),
          khigh = NULL
        )
      }
      # If neither works, use minimal parameters
      else {
        col_clustering <- hopach::hopach(
          data = t(expression_data)  # Transpose for column clustering
        )
      }
    }, error = function(e) {
      # Attempt with default parameters if specific ones fail
      message(paste("Error in column clustering:", e$message))
      message("Attempting column clustering with default parameters")
      tryCatch({
        col_clustering <<- hopach::hopach(data = t(expression_data))
      }, error = function(e2) {
        message(paste("Error in basic column clustering:", e2$message))
        message("Skipping column clustering")
      })
    })
  }
  
  # Handle the case where clustering fails completely
  if (is.null(row_clustering) && cluster_method %in% c("row", "both")) {
    message("Row clustering failed, using original row order")
    row_order <- 1:nrow(expression_data)
  } else if (!is.null(row_clustering)) {
    row_order <- row_clustering$final$order
  } else {
    row_order <- 1:nrow(expression_data)
  }
  
  if (is.null(col_clustering) && cluster_method %in% c("column", "both")) {
    message("Column clustering failed, using original column order")
    col_order <- 1:ncol(expression_data)
  } else if (!is.null(col_clustering)) {
    col_order <- col_clustering$final$order
  } else {
    col_order <- 1:ncol(expression_data)
  }
  
  # Reorder data based on clustering
  ordered_data <- expression_data[row_order, col_order]
  
  return(list(
    expression_data = ordered_data,
    row_order = row_order,
    col_order = col_order,
    row_clustering = row_clustering,
    col_clustering = col_clustering
  ))
}

# Function to filter expression data based on RPKM cutoff
filter_by_rpkm <- function(expression_data, sample_metadata, rpkm_cutoff) {
  if (is.null(rpkm_cutoff)) {
    message("No RPKM-based filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Find columns containing "Rpkm" in their names
  rpkm_cols <- grep("Rpkm", colnames(expression_data), ignore.case = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No columns with 'Rpkm' in their names found, no filtering applied")
    return(expression_data)
  }
  
  # Keep rows where any RPKM column exceeds the cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  message(paste("Filtered out", sum(!keep_rows), "of", nrow(expression_data), "rows"))
  
  return(expression_data[keep_rows, ])
}

# Function to load and validate metadata file
load_metadata <- function(args) {
  message(paste("Loading metadata from", args$meta))
  
  # Validate metadata file exists
  if (!file.exists(args$meta)) {
    stop(paste("Metadata file does not exist:", args$meta))
  }
  
  # Determine file format
  if (grepl("\\.csv$", args$meta, ignore.case = TRUE)) {
    delimiter <- ","
  } else if (grepl("\\.tsv$", args$meta, ignore.case = TRUE)) {
    delimiter <- "\t"
  } else {
    # Default to CSV
    message(paste("Unknown metadata file extension -", args$meta, "- assuming CSV format"))
    delimiter <- ","
  }
  
  # Load metadata with error handling
  tryCatch({
    metadata <- read.table(args$meta, sep = delimiter, header = TRUE, row.names = 1)
    message(paste("Loaded metadata with", nrow(metadata), "samples and", ncol(metadata), "variables"))
    
    # Validate metadata structure
    if (nrow(metadata) == 0) {
      stop("Metadata file contains no samples")
    }
    
    # Check if metadata has batch column when batch correction is requested
    if (args$batchcorrection != "none") {
      if (!"batch" %in% colnames(metadata)) {
        warning("Batch correction requested but 'batch' column not found in metadata")
        message("Switching batch correction method to 'none'")
        args$batchcorrection <- "none"
      } else {
        # Ensure batch column is numeric
        if (!is.numeric(metadata$batch)) {
          message("Converting batch column to numeric")
          metadata$batch <- as.numeric(as.factor(metadata$batch))
        }
      }
    }
    
    # Check if all input sample names are in metadata
    if (!is.null(args$name)) {
      missing_samples <- args$name[!args$name %in% rownames(metadata)]
      if (length(missing_samples) > 0) {
        warning(paste("Some sample names are missing in metadata:", paste(missing_samples, collapse = ", ")))
        
        # Try to find matches with different case
        for (sample in missing_samples) {
          case_insensitive_match <- grep(paste0("^", sample, "$"), rownames(metadata), ignore.case = TRUE)
          if (length(case_insensitive_match) > 0) {
            message(paste("Found case-insensitive match for", sample, ":", rownames(metadata)[case_insensitive_match[1]]))
          }
        }
        
        message("Continuing with available samples...")
      }
    }
    
    return(metadata)
    
  }, error = function(e) {
    stop(paste("Error loading metadata:", e$message))
  })
}

# Function to prepare metadata for DESeq2
prepare_metadata_for_deseq <- function(metadata, count_colnames) {
  # Extract sample names from count data column names
  sample_pattern <- "^([^\\s]+)\\s.*$"
  sample_names <- gsub(sample_pattern, "\\1", count_colnames)
  
  # Ensure all samples in count data have metadata
  missing_samples <- sample_names[!sample_names %in% metadata$sample]
  
  if (length(missing_samples) > 0) {
    stop(paste("Metadata missing for samples:", paste(missing_samples, collapse=", ")))
  }
  
  # Subset metadata to only include samples in count data
  metadata_subset <- metadata[metadata$sample %in% sample_names, , drop = FALSE]
  
  # Ensure order matches count data
  metadata_subset <- metadata_subset[match(sample_names, metadata_subset$sample), , drop = FALSE]
  
  # Ensure column order matches count data
  if (!identical(metadata_subset$sample, sample_names)) {
    logger::warn("Metadata sample order doesn't match count data. Reordering.")
    metadata_subset <- metadata_subset[match(sample_names, metadata_subset$sample), , drop = FALSE]
  }
  
  # Convert character columns to factors for DESeq2
  for (col in colnames(metadata_subset)) {
    if (is.character(metadata_subset[[col]])) {
      metadata_subset[[col]] <- as.factor(metadata_subset[[col]])
    }
  }
  
  return(metadata_subset)
}

#' Perform quality control and filtering on expression data
#' @param data Expression data frame
#' @param args Command line arguments
#' @return Filtered expression data frame
#' @export
perform_qc_filtering <- function(data, args) {
  message("Performing quality control and filtering...")
  
  # Initial count
  initial_count <- nrow(data)
  message(paste("Initial gene count:", initial_count))
  
  # Apply RPKM filtering if specified
  if (!is.null(args$rpkm_cutoff)) {
    message(paste("Applying RPKM cutoff of", args$rpkm_cutoff))
    
    # Find RPKM columns
    rpkm_cols <- grep("Rpkm$", colnames(data), value = TRUE)
    
    if (length(rpkm_cols) == 0) {
      warning("No RPKM columns found for filtering")
    } else {
      # Keep genes where at least one sample has RPKM >= cutoff
      rpkm_filter <- apply(data[, rpkm_cols, drop = FALSE], 1, function(x) any(x >= args$rpkm_cutoff))
      data <- data[rpkm_filter, ]
      message(paste("After RPKM filtering:", nrow(data), "genes remain"))
      message(paste(initial_count - nrow(data), "genes removed by RPKM filtering"))
    }
  }
  
  # Handle test mode if enabled
  if (args$test_mode) {
    message("Test mode enabled - sampling 500 genes")
    if (nrow(data) > 500) {
      # Use stratified sampling to keep representative genes
      set.seed(42) # for reproducibility
      data <- data[sample(nrow(data), 500), ]
      message("Sampled 500 genes for test mode")
    }
  }
  
  # Check for duplicate gene IDs
  if (any(duplicated(data$GeneId))) {
    dup_count <- sum(duplicated(data$GeneId))
    warning(paste("Found", dup_count, "duplicate gene IDs"))
    
    # Make gene IDs unique
    message("Making gene IDs unique")
    data$GeneId <- make.unique(as.character(data$GeneId))
  }
  
  message(paste("Final filtered data has", nrow(data), "genes"))
  return(data)
} 