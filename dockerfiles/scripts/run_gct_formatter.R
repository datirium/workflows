#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

suppressMessages(library(argparse))
suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape2))



##########################################################################################
#
# v0.0.1
# - use function to covert counts, row metadata, and col metadata into GCT format
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
anno_file <- args[1]
outdir <- args[2]
datadir <- args[3]
genome <- args[4]
threads <- args[5]








export_gct <- function(row_metadata=NULL, col_metadata=NULL, counts_mat, location){
    base::tryCatch(
        expr = {
            if (!is.null(row_metadata)){
                row_metadata <- row_metadata %>%
                                tibble::rownames_to_column("id") %>%
                                dplyr::mutate_at("id", base::as.vector)
                counts_mat <- counts_mat[row_metadata$rid, ]                      # to guarantee the order and number of rows
            }
            if (!is.null(col_metadata)){
                col_metadata <- col_metadata %>%
                                tibble::rownames_to_column("id") %>%
                                dplyr::mutate_at("id", base::as.vector)
                counts_mat <- counts_mat[, col_metadata$cid]                      # to guarantee the order and number of columns
            }
            gct_data <- methods::new(
                "GCT",
                mat=counts_mat,
                rdesc=row_metadata,                                              # can be NULL
                cdesc=col_metadata                                               # can be NULL
            )
            cmapR::write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            base::print(base::paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export GCT data to ", location, " with error - ", e, sep=""))
        }
    )
}


print("Exporting filtered normalized cell read counts to GCT format")
export_gct(                                    # will be used by Morpheus
    row_metadata = row_metadata,    # includes gene list, gene name, and location as row names
    col_metadata = column_metadata, # includes experiment type, sample name, tss window, and data type as row names
    counts_mat = counts_mat,
    location = paste(args$output, "_heatmap.gct", sep = "")
)
