#!/usr/bin/env Rscript

# --- Main workflow functions ---

#' Initialize the environment for DESeq2 LRT Step 1 analysis
#'
#' Loads required libraries, sources dependency files, and configures environment
initialize_environment <- function() {
  # Display startup message
  message("Starting DESeq2 LRT Step 1 Analysis")
  message("Working directory:", getwd())
  
  # Print command line arguments for debugging purposes
  args <- commandArgs(trailingOnly = TRUE)
  message("Command line arguments received:")
  message(paste(args, collapse = " "))
  
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

  # Source DESeq2 LRT Step 1 specific functions
  source_with_fallback("functions/deseq2_lrt_step_1/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
  source_with_fallback("functions/deseq2_lrt_step_1/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
  source_with_fallback("functions/deseq2_lrt_step_1/deseq2_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
  source_with_fallback("functions/deseq2_lrt_step_1/contrast_generation.R", "/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")
  
  # Load required libraries with clear error messages
  tryCatch({
    message("Loading required libraries...")
    suppressPackageStartupMessages({
      required_packages <- c(
        "DESeq2",
        "BiocParallel",
        "data.table",
        "ggplot2",
        "plotly",
        "sva",
        "hopach",
        "stringr"
      )
      
      for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(paste("Required package not found:", pkg))
        }
        message(paste("Loading package:", pkg))
        library(pkg, character.only = TRUE)
      }
    })
  }, error = function(e) {
    stop(paste("Error loading libraries:", e$message))
  })

  # Configure R options
  configure_r_options()
  
  # Configure plot theme
  configure_plot_theme()
  
  log_message("Environment initialized for DESeq2 LRT Step 1 analysis")
}

# Load and validate metadata
load_and_validate_metadata <- function(args) {
  message("Loading metadata...")
  
  # Get the file delimiter
  delimiter <- check_file_delimiter(args$meta)
  
  # Load metadata
  metadata_df <- read.table(
    args$meta,
    sep = delimiter,
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  # Clean metadata column and row names
  colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
  rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  
  message(glue::glue("Loaded metadata for {nrow(metadata_df)} samples with {ncol(metadata_df)} covariates"))
  
  # Check design formulas
  design_formula <- as.formula(args$design)
  
  # Apply comprehensive metadata validation using the common utility function
  metadata_df <- validate_metadata(metadata_df, args$batchcorrection, design_formula)
  
  # Add formulas to metadata for convenience
  attr(metadata_df, "design_formula") <- design_formula
  attr(metadata_df, "reduced_formula") <- as.formula(args$reduced)
  
  return(metadata_df)
}

# Load and validate expression data
load_and_validate_expression_data <- function(args, metadata_df) {
  message("Loading expression data...")
  
  # Clean sample names for consistency
  clean_names <- clean_sample_names(args$name)
  
  # Trim any trailing whitespace which can cause issues
  clean_names <- trimws(clean_names)
  
  # Check if we have the correct number of sample names
  if (length(clean_names) == 0) {
    stop("No sample names provided. Please provide sample names with --name parameter.")
  }
  
  # Validate we have the right number of files and names
  if (length(args$input) != length(clean_names)) {
    warning(sprintf("Mismatch between number of input files (%d) and sample names (%d).",
                    length(args$input), length(clean_names)))
    
    # Handle two specific cases:
    # 1. More files than names: use file basename as names
    # 2. More names than files: truncate names to match files
    
    if (length(args$input) > length(clean_names)) {
      message("More input files than sample names. Using file basenames for missing names.")
      
      # Generate names from file paths for the missing slots
      missing_names_count <- length(args$input) - length(clean_names)
      file_basenames <- basename(args$input[(length(clean_names)+1):length(args$input)])
      file_basenames <- gsub("\\.(tsv|csv)$", "", file_basenames)
      
      # Append the generated names
      clean_names <- c(clean_names, file_basenames)
      message(sprintf("Added %d names from file basenames. Now have %d names.", 
                      missing_names_count, length(clean_names)))
    } else if (length(clean_names) > length(args$input)) {
      message("More sample names than input files. Truncating sample names list.")
      clean_names <- clean_names[1:length(args$input)]
    }
  }
  
  # Load expression data
  message(sprintf("Loading expression data for %d files with %d sample names", 
                  length(args$input), length(clean_names)))
  expression_data_df <- load_expression_data(args$input, clean_names, READ_COL, RPKM_COL, INTERSECT_BY)
  message(glue::glue("Loaded expression data for {nrow(expression_data_df)} genes from {length(args$input)} files"))
  
  # Apply RPKM filtering if specified
  rpkm_filtered_count <- NULL
  if (!is.null(args$rpkm_cutoff)) {
    message(glue::glue("Applying RPKM cutoff of {args$rpkm_cutoff}..."))
    initial_gene_count <- nrow(expression_data_df)
    
    expression_data_df <- filter_rpkm(expression_data_df, args$rpkm_cutoff)
    
    # Ensure at least some genes remain
    if (nrow(expression_data_df) == 0) {
      report_error(
        "No genes remaining after RPKM filtering!",
        details = "All genes were removed due to the RPKM cutoff being too high.",
        recommendations = c("Reduce the RPKM threshold to retain more genes.")
      )
    }
    
    # Calculate how many genes were removed
    rpkm_filtered_count <- initial_gene_count - nrow(expression_data_df)
    
    message(glue::glue("{rpkm_filtered_count} genes removed, {nrow(expression_data_df)} genes retained"))
  }
  
  # Process count data
  read_counts_columns <- grep(
    paste(READ_COL, sep = ""),
    colnames(expression_data_df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  message(glue::glue("Found {length(read_counts_columns)} read count columns"))
  
  # Check if GeneId column exists
  if (!INTERSECT_BY %in% colnames(expression_data_df)) {
    stop(paste("Required column", INTERSECT_BY, "not found in expression data"))
  }
  
  # Check for duplicate GeneId values
  if (any(duplicated(expression_data_df[[INTERSECT_BY]]))) {
    dup_genes <- expression_data_df[[INTERSECT_BY]][duplicated(expression_data_df[[INTERSECT_BY]])]
    message(glue::glue("Warning: Found {length(dup_genes)} duplicate gene identifiers"))
    message(glue::glue("First few duplicates: {paste(head(dup_genes), collapse=', ')}"))
    
    # Make gene IDs unique
    expression_data_df[[INTERSECT_BY]] <- make.unique(as.character(expression_data_df[[INTERSECT_BY]]))
    message("Made gene identifiers unique")
  }
  
  # Extract read count data with safer approach
  tryCatch({
    # First extract only needed columns
    read_counts_subset <- expression_data_df[, c(INTERSECT_BY, read_counts_columns)]
    
    # Check for duplicate column names
    if (any(duplicated(colnames(read_counts_subset)))) {
      dup_cols <- colnames(read_counts_subset)[duplicated(colnames(read_counts_subset))]
      stop(paste("Duplicate column names in read count data:", paste(dup_cols, collapse=", ")))
    }
    
    # Make column names clean for downstream processing
    temp_colnames <- colnames(read_counts_subset)
    temp_colnames[-1] <- lapply(temp_colnames[-1], function(s) {
      # Extract the sample name (before the space)
      parts <- unlist(strsplit(s, " ", fixed = TRUE))
      if (length(parts) > 1) {
        return(parts[1])
      } else {
        return(s)
      }
    })
    colnames(read_counts_subset) <- temp_colnames
    
    # Now convert GeneId to row names
    # Clone the data frame
    read_counts_data_df <- read_counts_subset
    
    # Convert GeneId column to uppercase
    read_counts_data_df[[INTERSECT_BY]] <- toupper(read_counts_data_df[[INTERSECT_BY]])
    
    # Keep only first instance of each GeneId
    read_counts_data_df <- read_counts_data_df[!duplicated(read_counts_data_df[[INTERSECT_BY]]), ]
    
    # Check for duplicated gene IDs again after toupper conversion
    if (any(duplicated(read_counts_data_df[[INTERSECT_BY]]))) {
      dup_genes <- read_counts_data_df[[INTERSECT_BY]][duplicated(read_counts_data_df[[INTERSECT_BY]])]
      stop(paste("Duplicate gene IDs after uppercase conversion:", paste(head(dup_genes), collapse=", ")))
    }
    
    # Set row names safely
    rownames(read_counts_data_df) <- read_counts_data_df[[INTERSECT_BY]]
    read_counts_data_df <- read_counts_data_df[, -1, drop = FALSE] # Remove GeneId column
    
  }, error = function(e) {
    stop(paste("Error processing count data:", e$message))
  })
  
  # Clean count data column names with improved method
  cleaned_colnames <- clean_sample_names(colnames(read_counts_data_df))
  
  # Check for duplicates in cleaned column names
  if (any(duplicated(cleaned_colnames))) {
    dup_cols <- cleaned_colnames[duplicated(cleaned_colnames)]
    message(paste("Warning: Duplicate column names after cleaning:", paste(dup_cols, collapse=", ")))
    
    # Make column names unique
    cleaned_colnames <- make.unique(cleaned_colnames)
    message("Made column names unique")
  }
  
  # Apply cleaned column names
  colnames(read_counts_data_df) <- cleaned_colnames
  
  # Verify sample name consistency between metadata and counts
  validate_sample_consistency(metadata_df, read_counts_data_df)
  
  # Reorder count data columns to match metadata
  read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]
  
  # Return all results
  return(list(
    expression_df = expression_data_df,
    counts = read_counts_data_df,
    design_formula = attr(metadata_df, "design_formula"),
    reduced_formula = attr(metadata_df, "reduced_formula"),
    rpkm_filtered_count = rpkm_filtered_count
  ))
}

# Main execution function to organize workflow
run_deseq_analysis <- function(args) {
  # Start time tracking for the entire analysis
  total_start_time <- proc.time()
  
  # Setup and validation
  message("=== DESeq2 Analysis Pipeline ===")
  message(glue::glue("Analysis started at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
  
  # 1. Load and validate metadata
  metadata_df <- load_and_validate_metadata(args)
  
  # 2. Load and validate expression data
  expression_data <- load_and_validate_expression_data(args, metadata_df)
  
  # 3. Process count data and apply batch correction
  count_data_results <- process_count_data(args, expression_data$counts, metadata_df, 
                                         expression_data$design_formula)
  
  # 4. Run DESeq2 analysis
  deseq_results <- run_deseq2(count_data_results$countData, 
                            metadata_df, 
                            count_data_results$design_formula,
                            args$reduced,
                            args)
  
  # Add required fields to deseq_results that might be missing
  if (is.list(deseq_results)) {
    # Add design and reduced formulas if not already present
    if (is.null(deseq_results$design_formula)) {
      deseq_results$design_formula <- count_data_results$design_formula
    }
    
    if (is.null(deseq_results$reduced_formula)) {
      deseq_results$reduced_formula <- as.formula(args$reduced)
    }
    
    # Rename results to lrt_res for compatibility with export_results
    if (!is.null(deseq_results$results) && is.null(deseq_results$lrt_res)) {
      deseq_results$lrt_res <- deseq_results$results
    }
    
    # Add normCounts if only normalized_counts exists
    if (!is.null(deseq_results$normalized_counts) && is.null(deseq_results$normCounts)) {
      deseq_results$normCounts <- deseq_results$normalized_counts
    }
  }
  
  # Verify args has output_prefix
  if (is.null(args$output_prefix) && !is.null(args$output)) {
    args$output_prefix <- args$output
  }
  
  # 5. Process and export results with proper error handling
  tryCatch({
    export_results(deseq_results, expression_data$expression_df, 
                 metadata_df, args, count_data_results$batch_warning,
                 expression_data$rpkm_filtered_count)
  }, error = function(e) {
    log_message(paste("Error in export_results:", e$message), "ERROR")
    log_message("Continuing with analysis despite export error", "WARNING")
  })
  
  # Report total elapsed time
  total_elapsed <- proc.time() - total_start_time
  message(glue::glue("=== Analysis completed in {round(total_elapsed['elapsed']/60, 1)} minutes ==="))
  message(glue::glue("Finished at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
}

# Harmonize parameter names for compatibility
harmonize_parameters <- function(args) {
  # This function is retained for backward compatibility with any code
  # that might still expect the old parameter names
  # No action needed since we've handled parameter mapping in ArgumentParser
  return(args)
}

# Check required parameters for analysis
validate_analysis_params <- function(args) {
  # Check critical parameters
  critical_params <- c("input", "name", "meta", "design", "reduced")
  missing_params <- critical_params[!critical_params %in% names(args) | sapply(args[critical_params], is.null)]
  
  if (length(missing_params) > 0) {
    stop(paste("Missing critical parameters:", paste(missing_params, collapse=", ")))
  }
  
  # Validate batch correction parameter
  if (args$batchcorrection != "none") {
    # Check if metadata file exists
    if (!file.exists(args$meta)) {
      stop("Metadata file does not exist")
    }
    
    # Get the file delimiter
    delimiter <- check_file_delimiter(args$meta)
    
    # Read metadata to check for batch column
    metadata <- read.table(args$meta, sep=delimiter, header=TRUE)
    
    if (!"batch" %in% colnames(metadata)) {
      warning("Batch correction requested but 'batch' column not found in metadata. Batch correction will be disabled.")
      args$batchcorrection <- "none"
    }
  }
  
  return(args)
}

# Run the main process - adapter between main_with_memory_management and run_deseq_analysis
run_main_process <- function(args) {
  # Start timing
  start_time <- Sys.time()
  message("Starting DESeq2 LRT analysis process...")

  # Configure parallel processing based on thread count
  if (args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"), "CONFIG")
    BiocParallel::register(BiocParallel::MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode", "CONFIG")
  }

  # Run the main analysis workflow
  run_deseq_analysis(args)

  # Report completion and elapsed time
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  message(sprintf("DESeq2 LRT analysis completed in %.2f minutes", as.numeric(elapsed)))
}

# Main script execution with enhanced error handling
main <- function(args = NULL) {
  # Set up error handling
  tryCatch({
    # Parse arguments if not provided
    if (is.null(args)) {
      args <- get_args()
    }
    
    # Harmonize parameter names for compatibility
    args <- harmonize_parameters(args)
    
    # Validate analysis parameters
    args <- validate_analysis_params(args)
    
    # Setup debugging and logging
    if ("debug" %in% names(args) && args$debug) {
      enable_debug_mode(args)
    }
    
    # Configure BiocParallel
    BiocParallel::register(BiocParallel::MulticoreParam(args$threads))
    log_message(glue::glue("Using {args$threads} CPU threads for parallel processing"), "INFO")
    
    # Run the analysis
    run_deseq_analysis(args)
    
    # Report successful completion
    log_message("DESeq2 analysis completed successfully! ✓", "SUCCESS")
    
  }, error = function(e) {
    # Get detailed error information
    error_msg <- paste("Analysis failed:", e$message)
    
    # Set up detailed traceback capture
    options(rlang_backtrace_on_error = "full")
    
    # Get call stack information
    call_stack <- sys.calls()
    call_stack_formatted <- paste(format(call_stack), collapse = "\n")
    
    # Get traceback
    error_trace <- paste(capture.output(traceback()), collapse="\n")
    
    # Get session information
    session_info <- paste(capture.output(sessionInfo()), collapse="\n")
    
    # Log detailed error
    log_message(error_msg, "ERROR")
    log_message("Error occurred in call:", "ERROR")
    log_message(deparse(e$call), "ERROR")
    
    # Write detailed error report
    error_report_file <- file.path(getwd(), "deseq_lrt_error_report.txt")
    writeLines(
      paste0(
        "# DESeq2 LRT Analysis Error Report\n\n",
        "## Error Details\n\n",
        "**Timestamp:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
        "**Error Message:** ", e$message, "\n\n",
        "**Error Call:** ", deparse(e$call), "\n\n",
        
        "## Stack Trace\n\n```r\n", call_stack_formatted, "\n```\n\n",
        
        "## Error Traceback\n\n```\n", error_trace, "\n```\n\n",
        
        "## Session Information\n\n```\n", session_info, "\n```\n"
      ),
      con = error_report_file
    )
    
    log_message(paste("Detailed error report written to:", error_report_file), "ERROR")
    
    # Output simple error message to stderr for script integration
    message(error_msg)
    
    # Exit with error code
    quit(status = 1)
  })
}

# Wrapper function with memory management
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  
  # Make sure we have the log_message function available
  if (!exists("log_message", mode="function")) {
    # Define a basic version if not available
    log_message <- function(message, level="INFO") {
      cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), level, ":", message, "\n")
    }
  }
  
  # Report start time in a safe way
  log_message(paste("DESeq2 LRT Step 1 analysis started at", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Run main function in a safe way
  tryCatch({
    main()
  }, error = function(e) {
    message("Error in main function: ", e$message)
    quit(status = 1)
  })
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(paste("DESeq2 LRT Step 1 analysis completed at", format(end_time, "%Y-%m-%d %H:%M:%S")))
  log_message(paste("Total execution time:", round(as.numeric(duration), 2), "minutes"))
  
  # Clean up large objects to free memory
  gc(verbose = FALSE)
  
  # Define safe memory report function if it doesn't exist
  if (!exists("report_memory_usage", mode="function")) {
    report_memory_usage <- function(label) {
      mem_used <- round(pryr::mem_used() / 1024^2, 1)
      log_message(paste(label, "memory usage:", mem_used, "MB"))
    }
  }
  
  # Final memory report
  tryCatch({
    report_memory_usage("Final")
  }, error = function(e) {
    log_message("Could not report memory usage")
  })
}

# Main workflow module for DESeq2 LRT Step 1
#
# This module coordinates the entire DESeq2 LRT analysis workflow

# Source required modules with informative messages
source_with_message <- function(file_path, fallback_path = NULL) {
  if (file.exists(file_path)) {
    cat(paste("Sourcing:", file_path, "\n"))
    source(file_path)
    return(TRUE)
  } else if (!is.null(fallback_path) && file.exists(fallback_path)) {
    cat(paste("Using fallback path:", fallback_path, "\n"))
    source(fallback_path)
    return(TRUE)
  } else {
    cat(paste("WARNING: Could not find", file_path, "\n"))
    return(FALSE)
  }
}

#' Get the number of threads to use for parallel processing
#' @param args Command line arguments
#' @return Number of threads to use
#' @export
get_threads <- function(args = NULL) {
  # First check if args has threads parameter
  if (!is.null(args) && !is.null(args$threads) && is.numeric(args$threads)) {
    threads <- args$threads
  } else {
    # Try to get from command line
    cmd_args <- commandArgs(trailingOnly = TRUE)
    thread_index <- which(cmd_args == "--threads")
    
    if (length(thread_index) > 0 && thread_index < length(cmd_args)) {
      threads <- as.numeric(cmd_args[thread_index + 1])
    } else {
      # Default to 1
      threads <- 1
    }
  }
  
  # Ensure it's at least 1
  if (is.na(threads) || threads < 1) {
    threads <- 1
  }
  
  return(threads)
}

# Source common modules
source_with_message("/usr/local/bin/functions/common/utilities.R")
source_with_message("/usr/local/bin/functions/common/logging.R")
source_with_message("/usr/local/bin/functions/common/constants.R")
source_with_message("/usr/local/bin/functions/common/visualization.R")
source_with_message("/usr/local/bin/functions/common/clustering.R")
source_with_message("/usr/local/bin/functions/common/export_functions.R")

# Source DESeq2 LRT specific modules
source_with_message("/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
source_with_message("/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
source_with_message("/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
source_with_message("/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")

#' Initialize the environment for DESeQ2 LRT analysis
#' @export
initialize_environment <- function() {
  # Load required packages with clear error messages
  required_packages <- c(
    "DESeq2",
    "BiocParallel",
    "argparse",
    "pheatmap", 
    "dplyr",
    "tidyr",
    "data.table",
    "ggplot2",
    "plotly",
    "sva",
    "limma",
    "hopach",
    "rlang",
    "cmapR"
  )
  
  # Define aliases to avoid namespace conflicts
  if (requireNamespace("dplyr", quietly = TRUE)) {
    assign("mutate", dplyr::mutate, envir = .GlobalEnv)
    assign("filter", dplyr::filter, envir = .GlobalEnv)
    assign("select", dplyr::select, envir = .GlobalEnv)
  }
  
  # Load packages with verbose output for debugging
  cat("Loading required packages:\n")
  for (pkg in required_packages) {
    cat(paste("  - Loading", pkg, "... "))
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      cat("done\n")
    }, error = function(e) {
      cat("FAILED\n")
      cat(paste("    Error:", conditionMessage(e), "\n"))
      # Continue with other packages even if one fails
    })
  }
  
  # Configure parallel processing
  threads <- get_threads()
  cat(paste("Setting up parallel processing with", threads, "threads\n"))
  register(MulticoreParam(threads))
  
  cat("Environment initialized successfully\n")
}

#' Main workflow function for DESeq2 LRT analysis
#' @export
run_deseq2_lrt_workflow <- function() {
  # Print start time for better logging
  start_time <- Sys.time()
  cat("DESeq2 LRT Step 1 analysis started at: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  
  # Parse command line arguments
  args <- parse_cli_args()
  
  # Print parameters for debugging
  cat("Analysis parameters:\n")
  cat(paste("- Input files:", paste(basename(args$input), collapse=", "), "\n"))
  cat(paste("- Sample names:", paste(args$name, collapse=", "), "\n"))
  cat(paste("- Metadata file:", basename(args$meta), "\n"))
  cat(paste("- Design formula:", args$design, "\n"))
  cat(paste("- Reduced formula:", args$reduced, "\n"))
  cat(paste("- Output prefix:", args$output, "\n"))
  cat(paste("- RPKM cutoff:", args$rpkm_cutoff, "\n"))
  cat(paste("- Batch correction:", args$batchcorrection, "\n"))
  
  # Validate input files
  validate_inputs(args)
  
  # Load and process input data
  cat("Loading and processing input data...\n")
  input_data <- load_expression_data(args)
  metadata <- load_metadata(args)
  
  # Perform quality control and filtering
  cat("Performing quality control and filtering...\n")
  filtered_data <- perform_qc_filtering(input_data, args)
  
  # Run DESeq2 LRT analysis
  cat("Running DESeq2 LRT analysis...\n")
  deseq_results <- run_deseq2_lrt(filtered_data, metadata, args)
  
  # Generate contrasts and perform post-processing
  cat("Generating contrasts and post-processing results...\n")
  contrast_results <- generate_contrasts(deseq_results, args)
  
  # Export results
  cat("Exporting results...\n")
  export_results(contrast_results, args)
  
  # Verify outputs
  cat("Verifying outputs...\n")
  verify_outputs(args)
  
  # Print completion time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  cat("DESeq2 LRT Step 1 analysis completed at: ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total execution time: ", round(as.numeric(duration), 2), " minutes\n")
}

#' Validate that input files and parameters are correct
#' @param args Command line arguments
validate_inputs <- function(args) {
  cat("Validating inputs...\n")
  
  # Check required arguments
  required_args <- c("input", "name", "meta", "design", "reduced")
  missing_args <- required_args[sapply(required_args, function(x) is.null(args[[x]]))]
  
  if (length(missing_args) > 0) {
    stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Check that input files exist
  missing_files <- args$input[!file.exists(args$input)]
  if (length(missing_files) > 0) {
    stop(paste("Input files do not exist:", paste(missing_files, collapse=", ")))
  }
  
  # Check metadata file exists
  if (!file.exists(args$meta)) {
    stop(paste("Metadata file does not exist:", args$meta))
  }
  
  # Check if the number of input files matches the number of sample names
  if (length(args$input) != length(args$name)) {
    stop(paste("Number of input files (", length(args$input), 
               ") does not match number of sample names (", length(args$name), ")", sep=""))
  }
  
  cat("Input validation successful\n")
}

#' Verify that all required outputs were created
#' @param args Command line arguments
verify_outputs <- function(args) {
  output_prefix <- args$output

  # List of expected output files
  expected_files <- c(
    paste0(output_prefix, "_contrasts_table.tsv"),
    paste0(output_prefix, "_contrasts.rds"),
    paste0(output_prefix, "_gene_exp_table.tsv"),
    paste0(output_prefix, "_mds_plot.html"),
    paste0(output_prefix, "_counts_all.gct"),
    paste0(output_prefix, "_counts_filtered.gct"),
    paste0(output_prefix, "_lrt_result.md")
  )
  
  # Check each expected file
  missing_files <- c()
  for (file in expected_files) {
    if (file.exists(file)) {
      cat(paste("  ✓ Output verified:", basename(file), "\n"))
    } else {
      cat(paste("  ✗ Missing expected output:", basename(file), "\n"))
      missing_files <- c(missing_files, file)
    }
  }
  
  # Report missing files
  if (length(missing_files) > 0) {
    cat(paste("WARNING:", length(missing_files), "expected output files are missing\n"))
  } else {
    cat("All expected output files were created successfully\n")
  }
} 