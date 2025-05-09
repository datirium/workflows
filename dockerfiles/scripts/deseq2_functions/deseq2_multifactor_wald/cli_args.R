#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq2 LRT Step 2
#

#' Define and parse command line arguments
#'
#' @return Parsed arguments list
#' @export
get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 analysis using contrasts from previous LRT step")
  parser$add_argument(
    "--dsq_obj_data",
    help = "RDS file containing contrasts and expression data from step 1",
    required = TRUE,
    type = "character"
  )
  parser$add_argument(
    "--contrast_df",
    required = TRUE,
    help = "TSV file containing contrasts data",
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = "Batch correction method to use: 'none', 'combatseq', or 'model'",
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  parser$add_argument(
    "--contrast_indices",
    help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "FDR cutoff for significance filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Use lfcthreshold as the null hypothesis value in the results function call.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison: 'both', 'up', or 'down'.",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts.",
      "Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Type of scaling for expression data: 'minmax' or 'zscore'.",
      "Default: zscore"
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = paste(
      "Distance metric for HOPACH row clustering.",
      "Default: cosangle"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "cosangle"
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering.",
      "Default: euclid"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "euclid"
  )
  parser$add_argument(
    "--k",
    help = "Number of levels for Hopach clustering (1-15). Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help = "Maximum number of clusters at each level (2-9). Default: 5.",
    type = "integer",
    default = 5
  )
  parser$add_argument(
    "--output",
    help = "Output prefix. Default: deseq-lrt-step-2",
    type = "character",
    default = "deseq-lrt-step-2"
  )
  parser$add_argument(
    "--threads",
    help = "Number of threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run in test mode (first 500 rows only)",
    action = "store_true",
    default = FALSE
  )
  
  # Parse arguments with better error handling
  tryCatch({
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    # Get all command line arguments
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Initialize an empty list for our parsed arguments
    args <- list()
    
    # Make sure to explicitly extract required arguments first
    required_args <- c("dsq_obj_data", "contrast_df", "contrast_indices")
    for (req_arg in required_args) {
      req_flag <- paste0("--", req_arg)
      arg_idx <- which(all_args == req_flag)
      if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
        args[[req_arg]] <- all_args[arg_idx[1] + 1]
        message(paste("Directly extracted required argument:", req_arg, "=", args[[req_arg]]))
      }
    }
    
    # First pass: process all flags with values
    i <- 1
    while (i <= length(all_args)) {
      current_arg <- all_args[i]
      
      # Check if this is a flag argument (starts with --)
      if (grepl("^--", current_arg)) {
        arg_name <- sub("^--", "", current_arg)
        
        # Check if the next item exists and is not a flag
        if (i < length(all_args) && !grepl("^--", all_args[i + 1])) {
          arg_value <- all_args[i + 1]
          
          # Set the value
          args[[arg_name]] <- arg_value
          
          i <- i + 2  # Skip the value
        } else {
          # This is a boolean flag
          args[[arg_name]] <- TRUE
          i <- i + 1
        }
      } else {
        # This is a positional argument, move on
        i <- i + 1
      }
    }
    
    # Show what we parsed
    message("Manually parsed arguments:")
    for (arg_name in names(args)) {
      message(paste0("  ", arg_name, ": ", args[[arg_name]]))
    }
    
    return(args)
  })
  
  # Validate required arguments
  required_args <- c("dsq_obj_data", "contrast_df", "contrast_indices")
  missing_args <- required_args[!required_args %in% names(args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Process contrast indices
  if (!is.null(args$contrast_indices)) {
    args$contrast_indices <- as.numeric(unlist(strsplit(args$contrast_indices, ",")))
  }
  
  # Convert boolean string values if needed
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      if (is.character(args[[arg_name]])) {
        args[[arg_name]] <- toupper(args[[arg_name]]) %in% c("TRUE", "T", "YES", "Y", "1")
      }
    }
  }
  
  # Convert numeric values
  for (arg_name in c("fdr", "lfcthreshold", "k", "kmax", "threads")) {
    if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
      if (grepl("^[0-9.]+$", args[[arg_name]])) {
        args[[arg_name]] <- as.numeric(args[[arg_name]])
      }
    }
  }
  
  return(args)
} 