#!/usr/bin/env Rscript
#
# Common clustering functions for DESeq2 analysis
#

#' Generate clustered data using HOPACH with memory optimization
#'
#' @param expression_data Expression data matrix
#' @param by Clustering dimension ("row" or "col")
#' @param k Number of levels for clustering
#' @param kmax Maximum number of clusters at each level
#' @param dist Distance metric
#' @param scaling_type Type of scaling to apply
#' @param max_features Maximum number of features to use for clustering (for memory optimization)
#' @param use_subset Whether to use a subset of high-variance features for large datasets
#' @return List with clustered data results
#' @export
get_clustered_data <- function(expression_data, by = "row", k = 3, kmax = 5, dist = "cosangle", 
                              scaling_type = "zscore", max_features = 5000, use_subset = TRUE) {
  start_time <- proc.time()
  log_message(paste("Starting HOPACH clustering with", by, "clustering, using", scaling_type, "scaling"))

  if (!(by %in% c("row", "col"))) {
    stop("Invalid value for 'by'. Choose either 'row' or 'col'.")
  }

  # Get dimensions for logging
  n_rows <- nrow(expression_data)
  n_cols <- ncol(expression_data)
  log_message(paste("Input data dimensions:", n_rows, "rows x", n_cols, "columns"))
  
  # Large dataset handling
  orig_expression_data <- expression_data  # Keep original for later
  subset_used <- FALSE
  row_subset <- NULL
  
  # For row clustering on large datasets, use a subset of high-variance features
  if (by == "row" && n_rows > max_features && use_subset) {
    log_message(paste("Large dataset detected. Using top", max_features, "features by variance for clustering"))
    
    # Calculate row variances
    row_vars <- apply(expression_data, 1, var)
    
    # Select high-variance rows
    row_subset <- order(row_vars, decreasing = TRUE)[1:max_features]
    expression_data <- expression_data[row_subset, , drop = FALSE]
    
    log_message(paste("Subset created with dimensions:", nrow(expression_data), "rows x", ncol(expression_data), "columns"))
    subset_used <- TRUE
  }

  # If clustering by columns, transpose so that columns become rows for scaling.
  orientation_transposed <- FALSE
  if (by == "col") {
    log_message("Transposing expression data for column clustering")
    expression_data <- t(expression_data)
    orientation_transposed <- TRUE
  }

  # Apply scaling - more memory efficient implementation
  log_message(paste("Applying", scaling_type, "scaling"))
  expression_data <- switch(scaling_type,
    "none" = expression_data,  # No scaling
    "minmax" = {
      # Min-max scaling with more efficient implementation
      log_message("Applying min-max scaling")
      t(apply(expression_data, 1, function(x) {
        min_val <- min(x, na.rm = TRUE)
        max_val <- max(x, na.rm = TRUE)
        if (max_val == min_val) return(rep(0, length(x)))
        return((x - min_val) / (max_val - min_val) * 4 - 2)  # Scale to [-2, 2]
      }))
    },
    "zscore" = {
      # Z-score scaling with NA handling
      log_message("Applying z-score scaling")
      scaled_data <- t(scale(t(expression_data), center = TRUE, scale = TRUE))
      # Handle zero-variance rows
      scaled_data[is.na(scaled_data)] <- 0
      scaled_data
    },
    stop("Invalid scaling type. Choose 'none', 'minmax', or 'zscore'.")
  )

  # If data was transposed for column scaling, transpose back to original orientation.
  if (orientation_transposed) {
    log_message("Transposing scaled data back to original orientation")
    expression_data <- t(expression_data)
  }

  # Run garbage collection before HOPACH to free memory
  gc()
  log_message(paste0("Running HOPACH for ", nrow(expression_data), " features"))
  
  # Run HOPACH with error handling
  hopach_results <- tryCatch({
    hopach::hopach(expression_data,
                  verbose = TRUE,
                  K       = k,
                  kmax    = kmax,
                  khigh   = kmax,
                  d = dist)
  }, error = function(e) {
    log_message(paste("HOPACH clustering failed:", e$message), "ERROR")
    
    # Try with different parameters if it failed
    if (dist == "cosangle") {
      log_message("Retrying with 'euclid' distance", "WARNING")
      hopach::hopach(expression_data,
                    verbose = TRUE,
                    K       = k,
                    kmax    = kmax,
                    khigh   = kmax,
                    d = "euclid")
    } else {
      # Rethrow the error if we already tried with a different distance
      stop("HOPACH clustering failed with multiple distance metrics")
    }
  })

  log_message("Parsing cluster labels")
  # Final labels (no prefix 'c' in hopach by default)
  final_labels      <- hopach_results$clustering$labels[hopach_results$clustering$order]
  # Convert to character
  final_labels_char <- as.character(final_labels)

  # Each digit in the label corresponds to a level.
  # Determine max number of levels
  max_levels <- max(nchar(final_labels_char))

  # Create a data frame for clusters
  clusters           <- data.frame(Label = final_labels_char, stringsAsFactors = FALSE)
  
  # If we used a subset for row clustering, need to map back to full dataset
  if (subset_used && by == "row") {
    log_message("Mapping cluster results back to full dataset")
    
    # Get the ordered row indices from the subset
    ordered_subset_indices <- row_subset[hopach_results$clustering$order]
    
    # Create new order vector for all rows, putting non-clustered rows at the end
    all_indices <- seq_len(n_rows)
    non_subset_indices <- setdiff(all_indices, row_subset)
    full_order <- c(ordered_subset_indices, non_subset_indices)
    
    # Set rownames for clusters from ordered subset
    rownames(clusters) <- rownames(orig_expression_data)[ordered_subset_indices]
    
    # Create NA entries for non-clustered rows
    non_clustered <- data.frame(
      Label = rep(NA, length(non_subset_indices)), 
      stringsAsFactors = FALSE
    )
    rownames(non_clustered) <- rownames(orig_expression_data)[non_subset_indices]
    
    # Combine clustered and non-clustered
    clusters <- rbind(clusters, non_clustered)
    
    # Use the original data, but in the new order
    expression_data <- orig_expression_data[full_order, , drop = FALSE]
    
    # Update the order to include all indices
    hopach_results$clustering$order <- full_order
  } else {
    # Normal case - use the order from HOPACH directly
    rownames(clusters) <- rownames(expression_data)[hopach_results$clustering$order]
  }

  # Split labels into levels efficiently
  # For each label, we split into characters and assign to new columns
  level_columns <- vector("list", max_levels)
  for (i in 1:max_levels) {
    # Extract the ith character from each label, or NA if too short
    level_columns[[i]] <- substr(clusters$Label, i, i)
    level_columns[[i]][nchar(clusters$Label) < i] <- NA
  }
  
  # Convert to data frame
  level_data <- as.data.frame(level_columns, stringsAsFactors = FALSE)
  colnames(level_data) <- paste0("Cluster_Level_", seq_len(max_levels))

  # Combine into clusters
  clusters <- cbind(clusters, level_data)

  # Rename the Label column
  clusters$HCL   <- paste0("c", clusters$Label)
  clusters$Label <- NULL

  end_time <- proc.time() - start_time
  log_message(paste("HOPACH clustering completed in", round(end_time["elapsed"], 2), "seconds"))

  return(list(
    order      = hopach_results$clustering$order,
    expression = expression_data,
    clusters   = clusters,
    subset_used = subset_used
  ))
}

#' Cluster and reorder data based on clustering results
#'
#' @param normCounts Normalized count matrix
#' @param col_metadata Column metadata
#' @param row_metadata Row metadata
#' @param args Command-line arguments
#' @return List with clustered and reordered data
#' @export
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {
  start_time <- proc.time()

  if (args$cluster != "none") {
    # Determine clustering parameters based on dataset size
    if (args$test_mode) {
      k    <- 2
      kmax <- 2
      max_features <- 100
    } else {
      k    <- args$k
      kmax <- args$kmax
      
      # Adjust max_features based on dataset size
      n_rows <- nrow(normCounts)
      if (n_rows > 20000) {
        max_features <- 5000
        log_message(paste("Large dataset detected with", n_rows, "rows. Limiting clustering to top 5000 features."), "WARNING")
      } else if (n_rows > 10000) {
        max_features <- 8000
        log_message(paste("Medium-large dataset detected with", n_rows, "rows. Limiting clustering to top 8000 features."), "INFO")
      } else {
        max_features <- n_rows
      }
    }
    
    # Column clustering if requested
    if (args$cluster == "column" || args$cluster == "both") {
      log_message("Performing column clustering")
      clustered_data_cols <- get_clustered_data(
        normCounts, 
        by = "col", 
        k = k, 
        kmax = kmax, 
        scaling_type = args$scaling_type, 
        dist = args$columndist,
        max_features = max_features,
        use_subset = FALSE  # Don't use subset for column clustering - usually manageable
      )
      normCounts          <- normCounts[, clustered_data_cols$order, drop = FALSE]
      col_metadata        <- col_metadata[clustered_data_cols$order, , drop = FALSE]
      # After reordering, cbind cluster info
      col_metadata        <- cbind(col_metadata, clustered_data_cols$clusters)
    }
    
    # Row clustering if requested
    if (args$cluster == "row" || args$cluster == "both") {
      log_message("Performing row clustering")
      clustered_data_rows <- get_clustered_data(
        normCounts, 
        by = "row", 
        k = k, 
        kmax = kmax, 
        scaling_type = args$scaling_type, 
        dist = args$rowdist,
        max_features = max_features,
        use_subset = TRUE  # Use subset for row clustering if dataset is large
      )
      normCounts          <- clustered_data_rows$expression
      
      # Handle row metadata based on whether we used a subset
      if (clustered_data_rows$subset_used) {
        log_message("Reordering row metadata to match clustered subset")
        # Use the order from clustering to reorder row metadata
        row_metadata <- row_metadata[clustered_data_rows$order, , drop = FALSE]
      } else {
        # Standard reordering
        row_metadata <- row_metadata[clustered_data_rows$order, , drop = FALSE]
      }
      
      # Add cluster annotations
      row_metadata <- cbind(row_metadata, clustered_data_rows$clusters)
    }
  } else {
    log_message("No clustering requested")
  }

  end_time <- proc.time() - start_time
  log_message(paste("Clustering and reordering completed in", round(end_time["elapsed"], 2), "seconds"))

  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
} 