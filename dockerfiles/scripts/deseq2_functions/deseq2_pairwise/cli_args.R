#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for parsing and validating command-line arguments
# for the main DESeq analysis workflow.
#
# Version: 0.1.0

#' Assert and validate command line arguments
#'
#' @param args The parsed arguments from ArgumentParser
#' @return Modified args with validated and processed values
assert_args <- function(args) {
  log_message("Checking input parameters")
  
  # Process aliases if not provided
  if (is.null(args$untreated_sample_names) | is.null(args$treated_sample_names)) {
    log_message("--untreated_sample_names or --treated_sample_names were not set, using default values based on expression file names")
    
    args$untreated_sample_names <- character(0)
    for (i in 1:length(args$untreated_files)) {
      args$untreated_sample_names <- append(args$untreated_sample_names, head(unlist(
        strsplit(basename(args$untreated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
    
    args$treated_sample_names <- character(0)
    for (i in 1:length(args$treated_files)) {
      args$treated_sample_names <- append(args$treated_sample_names, head(unlist(
        strsplit(basename(args$treated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
  } else {
    # Verify correct number of aliases
    if ((length(args$untreated_sample_names) != length(args$untreated_files)) |
        (length(args$treated_sample_names) != length(args$treated_files))) {
      log_error("Not correct number of inputs provided for files and sample names")
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }

  # Check for minimum file requirements
  if (length(args$treated_files) == 1 || length(args$untreated_files) == 1) {
    log_warning("Only one file in a group. DESeq2 requires at least two replicates for accurate analysis.")
    args$batch_file <- NULL # reset batch_file to NULL. We don't need it for DESeq even if it was provided
  }

  # Process batch file if provided
  if (!is.null(args$batch_file)) {
    batch_metadata <- with_error_handling({
      read.table(
        args$batch_file,
        sep = get_file_type(args$batch_file),
        row.names = 1,
        col.names = c("name", "batch"),
        header = FALSE,
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(batch_metadata)) {
      log_error("Failed to read batch metadata file")
      args$batch_file <- NULL
      return(args)
    }
    
    log_message("Loaded batch metadata")
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    
    if (all(is.element(c(args$untreated_sample_names, args$treated_sample_names), rownames(batch_metadata)))) {
      args$batch_file <- batch_metadata # dataframe
    } else {
      log_warning("Missing values in batch metadata file. Skipping multi-factor analysis")
      log_debug(paste("Expected:", paste(c(args$untreated_sample_names, args$treated_sample_names), collapse=", ")))
      log_debug(paste("Found:", paste(rownames(batch_metadata), collapse=", ")))
      args$batch_file <- NULL
    }
  }

  # Convert boolean string values if they came as strings
  for (arg_name in c("use_lfc_thresh")) {
    if (!is.null(args[[arg_name]])) {
      args[[arg_name]] <- convert_to_boolean(args[[arg_name]], FALSE)
    }
  }

  return(args)
}

#' Parse command line arguments for DESeq analysis
#'
#' @return Parsed and validated argument list
get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2)")
  
  # Input file parameters
  parser$add_argument(
    "-u", "--untreated_files",
    help = "Untreated (condition 1) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-t", "--treated_files",
    help = "Treated (condition 2) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-ua", "--untreated_sample_names",
    help = "Unique aliases for untreated (condition 1) expression files. Default: basenames of -u without extensions",
    type = "character",
    nargs = "*"
  )
  parser$add_argument(
    "-ta", "--treated_sample_names",
    help = "Unique aliases for treated (condition 2) expression files. Default: basenames of -t without extensions",
    type = "character",
    nargs = "*"
  )
  
  # Condition naming parameters
  parser$add_argument(
    "-un", "--untreated_name",
    help = "Name for untreated (condition 1), use only letters and numbers",
    type = "character",
    default = "untreated"
  )
  parser$add_argument(
    "-tn", "--treated_name",
    help = "Name for treated (condition 2), use only letters and numbers",
    type = "character",
    default = "treated"
  )
  
  # Batch correction parameters
  parser$add_argument(
    "-bf", "--batch_file",
    help = paste(
      "Metadata file for multi-factor analysis. Headerless TSV/CSV file.",
      "First column - names from --untreated_sample_names and --treated_sample_names, second column - batch group name.",
      "Default: None"
    ),
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.",
      "- 'model' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  
  # Statistical and filtering parameters
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--rpkm_cutoff",
    help = paste(
      "RPKM cutoff for filtering genes. Genes with RPKM values below this threshold will be excluded from the analysis.",
      "Default: NULL (no filtering)"
    ),
    type = "integer",
    default = NULL
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold for greaterAbs and |β| < lfcThreshold for lessAbs, with p-values being two-tailed or maximum of the upper and lower tests, respectively); ",
      "'up' for upregulated genes (β > lfcThreshold in condition2 compared to condition1); ",
      "'down' for downregulated genes (β < -lfcThreshold in condition2 compared to condition1). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  
  # Clustering parameters
  parser$add_argument(
    "--cluster_method",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Specifies the type of scaling to be applied to the expression data.",
      "- 'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].",
      "- 'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.",
      "- Default: zscore"
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--row_distance",
    help = paste(
      "Distance metric for HOPACH row clustering. Ignored if --cluster_method is not",
      "provided. Default: cosangle"
    ),
    type = "character",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--column_distance",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster_method is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--k_hopach",
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax_hopach",
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5.",
    type = "integer",
    default = 5
  )
  
  # Output parameters
  parser$add_argument(
    "-o", "--output_prefix",
    help = "Output prefix. Default: deseq",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-d", "--digits",
    help = "Precision, number of digits to print. Default: 3",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "-p", "--threads",
    help = "Number of threads to use for parallel processing. Default: 1",
    type = "integer",
    default = 1
  )
  
  # Parse arguments with better error handling
  tryCatch({
    args <- parser$parse_args()
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    # Get all command line arguments
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Initialize an empty list for our parsed arguments
    args <- list()
    
    # Initialize arrays for multi-value arguments
    array_args <- c("untreated_files", "treated_files", "untreated_sample_names", "treated_sample_names")
    for (arg in array_args) {
      args[[arg]] <- character(0)
    }
    
    # Make sure to explicitly extract required arguments first
    required_args <- c("untreated_files", "treated_files")
    for (req_arg in required_args) {
      req_flag <- paste0("--", req_arg)
      short_flag <- if(req_arg == "untreated_files") "-u" else if(req_arg == "treated_files") "-t" else NULL
      
      # Try long form first
      arg_indices <- which(all_args == req_flag)
      if (length(arg_indices) > 0) {
        for (idx in arg_indices) {
          if (idx < length(all_args) && !grepl("^-", all_args[idx + 1])) {
            args[[req_arg]] <- c(args[[req_arg]], all_args[idx + 1])
            # Check if more values follow
            i <- idx + 2
            while(i <= length(all_args) && !grepl("^-", all_args[i])) {
              args[[req_arg]] <- c(args[[req_arg]], all_args[i])
              i <- i + 1
            }
          }
        }
      }
      
      # Try short form if available
      if (!is.null(short_flag)) {
        arg_indices <- which(all_args == short_flag)
        if (length(arg_indices) > 0) {
          for (idx in arg_indices) {
            if (idx < length(all_args) && !grepl("^-", all_args[idx + 1])) {
              args[[req_arg]] <- c(args[[req_arg]], all_args[idx + 1])
              # Check if more values follow
              i <- idx + 2
              while(i <= length(all_args) && !grepl("^-", all_args[i])) {
                args[[req_arg]] <- c(args[[req_arg]], all_args[i])
                i <- i + 1
              }
            }
          }
        }
      }
      
      if (length(args[[req_arg]]) > 0) {
        message(paste("Directly extracted required argument:", req_arg, "=", paste(args[[req_arg]], collapse=", ")))
      }
    }
    
    # Process sample names similarly
    optional_array_args <- c(
      "untreated_sample_names" = "--untreated_sample_names", 
      "treated_sample_names" = "--treated_sample_names"
    )
    optional_short_args <- c(
      "untreated_sample_names" = "-ua", 
      "treated_sample_names" = "-ta"
    )
    
    for (arg_name in names(optional_array_args)) {
      long_flag <- optional_array_args[arg_name]
      short_flag <- optional_short_args[arg_name]
      
      # Try long form
      arg_indices <- which(all_args == long_flag)
      if (length(arg_indices) > 0) {
        for (idx in arg_indices) {
          if (idx < length(all_args) && !grepl("^-", all_args[idx + 1])) {
            args[[arg_name]] <- c(args[[arg_name]], all_args[idx + 1])
            # Check if more values follow
            i <- idx + 2
            while(i <= length(all_args) && !grepl("^-", all_args[i])) {
              args[[arg_name]] <- c(args[[arg_name]], all_args[i])
              i <- i + 1
            }
          }
        }
      }
      
      # Try short form
      arg_indices <- which(all_args == short_flag)
      if (length(arg_indices) > 0) {
        for (idx in arg_indices) {
          if (idx < length(all_args) && !grepl("^-", all_args[idx + 1])) {
            args[[arg_name]] <- c(args[[arg_name]], all_args[idx + 1])
            # Check if more values follow
            i <- idx + 2
            while(i <= length(all_args) && !grepl("^-", all_args[i])) {
              args[[arg_name]] <- c(args[[arg_name]], all_args[i])
              i <- i + 1
            }
          }
        }
      }
      
      if (length(args[[arg_name]]) > 0) {
        message(paste("Extracted optional array argument:", arg_name, "=", paste(args[[arg_name]], collapse=", ")))
      }
    }
    
    # Process remaining scalar arguments
    i <- 1
    while (i <= length(all_args)) {
      current_arg <- all_args[i]
      
      # Skip arguments we've already processed
      if (current_arg %in% c("-u", "--untreated_files", "-t", "--treated_files", "-ua", "--untreated_sample_names", "-ta", "--treated_sample_names")) {
        # Skip this flag and its values
        i <- i + 1
        while(i <= length(all_args) && !grepl("^-", all_args[i])) {
          i <- i + 1
        }
        next
      }
      
      # Check if this is a flag argument (starts with -)
      if (grepl("^--", current_arg) || grepl("^-[a-zA-Z]", current_arg)) {
        arg_name <- sub("^--", "", current_arg)
        if (grepl("^-[a-zA-Z]", current_arg)) {
          short_name <- sub("^-", "", current_arg)
          # Map short names to long names
          if (short_name == "tn") arg_name <- "treated_name"
          else if (short_name == "un") arg_name <- "untreated_name"
          else if (short_name == "bf") arg_name <- "batch_file"
          else arg_name <- short_name
        }
        
        # Check if the next item exists and is not a flag
        if (i < length(all_args) && !grepl("^-", all_args[i + 1])) {
          arg_value <- all_args[i + 1]
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
      if (length(args[[arg_name]]) > 1) {
        message(paste0("  ", arg_name, ": [", paste(head(args[[arg_name]], 3), collapse=", "), 
                       if(length(args[[arg_name]]) > 3) "..." else "", "] (", length(args[[arg_name]]), " items)"))
      } else {
        message(paste0("  ", arg_name, ": ", args[[arg_name]]))
      }
    }
    
    return(args)
  })
  
  # Validate arguments and set defaults
  args <- assert_args(args)
  
  # Convert numeric values
  for (arg_name in c("fdr", "lfcthreshold", "threads")) {
    if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
      if (grepl("^[0-9.]+$", args[[arg_name]])) {
        args[[arg_name]] <- as.numeric(args[[arg_name]])
      }
    }
  }
  
  return(args)
} 