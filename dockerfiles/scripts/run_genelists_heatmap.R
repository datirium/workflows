#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

suppressMessages(library(argparse))
suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(reshape2))
suppressMessages(library(morpheus))
suppressMessages(library(htmlwidgets))


##########################################################################################
#
# v0.0.1
# - use function to covert counts, row metadata, and col metadata into GCT format
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
output_row_metadata <- args[1]  # ex. 'output_row_metadata.tsv'
output_col_metadata <- args[2]  # ex. 'output_col_metadata.tsv'
output_counts <- args[3]        # ex. 'output_counts.tsv'
output_location <- args[4]      # ex. './'

# read in files to data frames (counts need to be in matrix)
row_metadata <- read.table(output_row_metadata,header=T,sep='\t')
column_metadata <- read.table(output_col_metadata,header=T,sep='\t')
df <- read.table(output_counts,header=T,sep='\t')
counts_mat <- acast(df, rid~cid, value.var="value")


# export to GCT format function
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
print("Exporting to GCT format")
export_gct(                                    # will be used by Morpheus
    row_metadata = row_metadata,    # includes gene list, gene name, and location as row names
    col_metadata = column_metadata, # includes experiment type, sample name, tss window, and data type as row names
    counts_mat = counts_mat,
    location = paste(output_location, "/heatmap.gct", sep = "")
)


print("Generating morpheus heatmap")
gct_data <- read.gct(paste(output_location, "/heatmap.gct", sep = ""))
is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}
morpheus_html <- morpheus(
    x=na.omit(gct_data$data),
    rowAnnotations=if(nrow(gct_data$rowAnnotations) == 0) NULL else gct_data$rowAnnotations,
    columnAnnotations=if(nrow(gct_data$columnAnnotations) == 0) NULL else gct_data$columnAnnotations %>% dplyr::mutate_if(is_all_numeric, as.numeric),
    colorScheme=list(scalingMode="fixed", stepped=FALSE, values=list(0, 99, 100,199), colors=list("white", "black", "white", "red")),
    rowSortBy=list(list(field="genelist_name", order=0)),
    columnSortBy=list(
        list(field="data_type", order=0),
        list(field="sample_name", order=0),
        list(field="tss_window", order=0)),
    tools=list(list(
        name="Sort/Group", 
        params=list(group_cols_by="sample_name")
    )),
    rowGroupBy=list(list(field="genelist_name")),
    columnGroupBy=list(list(field="sample_name")),
    rowSize=3,
    columnSize=3
)
html_location <- paste(output_location, "/heatmap.html", sep = "")
print(paste("Saving heatmap to", html_location))
htmlwidgets::saveWidget(
    morpheus_html,
    file=html_location
)
