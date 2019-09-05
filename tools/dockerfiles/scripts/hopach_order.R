#!/usr/bin/env Rscript
options(warn=-1)
options("width"=300)


suppressMessages(library(argparse))
suppressMessages(library(hopach))
suppressMessages(library(gplots))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))


##########################################################################################
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified by file extension
# *.csv - CSV
# *.tsv - TSV
# CSV is used by default
#
# if --method is "both" heatmap is built using expression_data_row. The column order is
# taken from the expression_data_col.
#
##########################################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


load_data_set <- function(filenames, suffixes, target_colname, intersect_by) {
    selected_data <- NULL
    updated_colnames <- c()
    for (i in 1:length(filenames)) {
        raw_data <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, stringsAsFactors=FALSE)
        print(paste("Load ", nrow(raw_data), " rows from ", filenames[i], sep=""))
        new_colname = paste(target_colname, suffixes[i], sep="_")
        updated_colnames = c(updated_colnames, new_colname)
        colnames(raw_data)[colnames(raw_data) == target_colname] <- new_colname
        if (is.null(selected_data)){
            selected_data <- raw_data
        } else {
            if (is.null(intersect_by)){
                selected_data <- cbind(selected_data, raw_data)
            } else {
                print(paste("Combine loaded data by ", paste(intersect_by, collapse=", "), sep=""))
                selected_data <- merge(selected_data, raw_data, by=intersect_by, sort = FALSE)
            }
        }
    }
    return (selected_data[,c(intersect_by, updated_colnames)])
}


cluster <- function(expression_data, logtransform, center, normalize, dist, output, suffix, transpose) {

    if (transpose){
        print("   Transpose expression data")
        expression_data = t(expression_data)
    }

    if (logtransform) { 
        print("   Log2 transform expression data")
        if (any(expression_data==0)){
            print(paste("Zero values are replaced by ", 1, sep=""))
            expression_data[expression_data==0] <- 1
        }
        expression_data = log2(expression_data)
    }

    if (!is.null(center)) {
        print(paste("   Center expression data by ", center, sep=""))
        if (center == "mean"){
            expression_data = expression_data - rowMeans(expression_data)    
        } else {
            expression_data = expression_data - rowMedians(data.matrix(expression_data))    
        }
    }

    if (normalize) {
        print("   Normalize expression data")
        std = sqrt(rowSums(expression_data^2))
        expression_data = expression_data/std
        expression_data = replace(expression_data, is.na(expression_data), 0)
    }

    print("   Create distance matrix")
    distance_matrix <- distancematrix(expression_data, dist)

    print("   Run hopach clustering")
    hopach_results <- hopach(expression_data, dmat=distance_matrix)

    distance_matrix_name = paste(output, "_", suffix, "_dist_matrix.png", sep="")
    png(distance_matrix_name,
        width = 5*300,
        height = 5*300,
        res = 300,
        pointsize = 8)
    dplot(distance_matrix, hopach_results, ord="final", main=paste("Distance matrix (", dist, ", ", suffix, ")", sep=""), showclusters=FALSE, col=colorRampPalette(args$palette)(n = 299))
    print(paste("   Save distance matrix plot to ", distance_matrix_name, sep=""))

    print("   Reorder expression data")
    expression_data <- expression_data[hopach_results$final$order,]
    
    if (transpose){
        print("   Transpose expression data")
        expression_data = t(expression_data)
    }

    labels_df = as.data.frame(hopach_results$final$labels)
    colnames(labels_df) = "label"
    labels_df = cbind(labels_df, "L" = outer(labels_df$label, 10^c((nchar(trunc(labels_df$label))[1]-1):0), function(a, b) a %/% b %% 10))
    labels_df = labels_df[, c(-1), drop = FALSE]

    return( list(expression=expression_data, labels=labels_df) )

}

# Parser
parser <- ArgumentParser(description='Hopach Clustering')
parser$add_argument("--input",           help='Input CSV/TSV files',                                    type="character", required="True", nargs='+')
parser$add_argument("--name",            help='Input aliases, the order corresponds to --input order. Default: basename of --input files', type="character", nargs='+')
parser$add_argument("--target",          help='Target column to be used by hopach clustering. Default: Rpkm',      type="character", default="Rpkm")
parser$add_argument("--combine",         help='Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand', type="character", nargs='+', default=c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand"))
parser$add_argument("--method",          help='Cluster method. Default: both',                          type="character", choices=c("row","column","both"), default="both")
parser$add_argument("--palette",         help='Palette color names. Default: red, black, green',        type="character", nargs='+', default=c("red", "black", "green"))
parser$add_argument("--output",          help='Output prefix. Default: hopach',                         type="character", default="./hopach")


parser$add_argument("--rowmin",          help='Exclude rows from clustering by the min value of a target column. Default: 0',          type="double",    default=0)
parser$add_argument("--rowkeep",         help='Append excluded rows to the output table after clustering is finished. Default: false', action='store_true')


parser$add_argument("--rowdist",         help='Distance metric for row clustering. Default: cosangle',  type="character", choices=c("cosangle","abscosangle","euclid","abseuclid","cor","abscor"), default="cosangle")
parser$add_argument("--coldist",         help='Distance metric for column clustering. Default: euclid', type="character", choices=c("cosangle","abscosangle","euclid","abseuclid","cor","abscor"), default="euclid")


parser$add_argument("--rowlogtransform", help='Log2 transform input data prior to running row clustering. Default: false',    action='store_true')
parser$add_argument("--collogtransform", help='Log2 transform input data prior to running column clustering. Default: false', action='store_true')


parser$add_argument("--rowcenter",       help='Center rows prior to running row clustering. Default: not centered',    type="character", choices=c("mean", "median"))
parser$add_argument("--colcenter",       help='Center columns prior to running column clustering. Default: not centered', type="character", choices=c("mean", "median"))


parser$add_argument("--rownorm",         help='Normalize rows prior to running row clustering. Default: not normalized',    action='store_true' )
parser$add_argument("--colnorm",         help='Normalize columns prior to running column clustering. Default: not normalized', action='store_true' )


args <- parser$parse_args(commandArgs(trailingOnly = TRUE))


# Set default value for --name if it wasn't provided
if(is.null(args$name)){
    for (i in 1:length(args$input)) {
        args$name = append(args$name, head(unlist(strsplit(basename(args$input[i]), ".", fixed = TRUE)), 1))
    }
}


# Load and combine input data
input_data <- load_data_set(args$input, args$name, args$target, args$combine)
expression_columns = c((length(args$combine)+1):ncol(input_data))

# Get rows to be filtered out from input data
print(paste("Apply filter to input data ", args$target, " >= ", args$rowmin, sep=""))
filtered_rows <- rownames(input_data[!rowSums(input_data[,expression_columns] < args$rowmin),])
print(paste("Number of rows after filtering ", length(filtered_rows), sep=""))

# Extract expression data
print("Extract expression data")
expression_data = input_data[filtered_rows, expression_columns]


if (args$method == "both" || args$method == "row"){
    print("Clustering by row")
    cluster_row = cluster(expression_data, args$rowlogtransform, args$rowcenter, args$rownorm, args$rowdist, args$output, "row", FALSE)
    expression_data_row = cluster_row$expression
    labels_row = cluster_row$labels
}

if (args$method == "both" || args$method == "column"){
    print("Clustering by columns")
    cluster_col = cluster(expression_data, args$collogtransform, args$colcenter, args$colnorm, args$coldist, args$output, "column", TRUE)
    expression_data_col = cluster_col$expression
    labels_col = cluster_col$labels
}


# Reorder rows and columns of input data
if (args$method == "both"){
    output_data <- input_data[rownames(expression_data_row), append(colnames(input_data)[1:length(args$combine)], colnames(expression_data_col))]
} else if (args$method == "row"){
    output_data <- input_data[rownames(expression_data_row),]
} else if (args$method == "column"){
    output_data <- input_data[, append(colnames(input_data)[1:length(args$combine)], colnames(expression_data_col))]
}


# Append row cluster labels
if (args$method == "both" || args$method == "row"){
    output_data <- cbind(output_data, labels_row)
}



# Add discarded rows at the bottom of output data
if (args$rowkeep){
    print("Append discarded rows at the bottom of the output table")
    if (args$method == "both" || args$method == "column"){
        discarded_data = input_data[!rownames(input_data) %in% filtered_rows, append(colnames(input_data)[1:length(args$combine)], colnames(expression_data_col))]
    } else {
        discarded_data = input_data[!rownames(input_data) %in% filtered_rows,]   
    }
    if (args$method == "both" || args$method == "row"){
        discarded_data[, colnames(labels_row)] <- 0
    }
    output_data <- rbind(output_data, discarded_data)
}


# Save main clustering results
output_data_name = paste(args$output, "_clustering.tsv", sep="")
write.table(output_data,
            file=output_data_name,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)
print(paste("Export hopach clustering results to ", output_data_name, sep=""))


# Save column cluster labels
if (args$method == "both" || args$method == "column"){
    labels_col = cbind(experiment = colnames(expression_data_col), labels_col)
    labels_col_name = paste(args$output, "_column_clustering_labels.tsv", sep="")
    write.table(labels_col,
                file=labels_col_name,
                sep="\t",
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE)
    print(paste("Export hopach column clustering results to ", labels_col_name, sep=""))
}


# Create and export heatmap to png file
heatmap_name = paste(args$output, "_heatmap.png", sep="")
if (args$method == "both"){
    heatmap_data <- data.matrix(expression_data_row[,colnames(expression_data_col)])
} else if (args$method == "row"){
    heatmap_data <- data.matrix(expression_data_row)
} else if (args$method == "column"){
    heatmap_data <- data.matrix(expression_data_col)
}
pheatmap(heatmap_data,
        cluster_row=FALSE,
        cluster_cols=FALSE,
        treeheight_col = 0,
        main = paste("Heatmap (clustering: ", args$method, ")", sep=""),
        color=colorRampPalette(args$palette)(n = 299),
        scale="none",
        border_color=FALSE,
        show_rownames=FALSE,
        filename=heatmap_name)
print(paste("Export heatmap to ", heatmap_name, sep=""))


dev.off()