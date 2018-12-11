#!/usr/bin/env Rscript
options(warn=-1)
options("width"=300)

suppressMessages(library(argparse))
suppressMessages(library(hopach))

get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}

# Parser
parser <- ArgumentParser(description='Hopach Ordering')
parser$add_argument("-i", "--input",     help='Input CSV/TSV genelist file',     type="character", required="True")
parser$add_argument("-d", "--dist",      help='Distance metric',                 type="character", choices=c("cosangle","abscosangle","euclid","abseuclid","cor","abscor"), default="euclid")
parser$add_argument("-t", "--transform", help='Log2 transform input data prior running hopach',    action='store_true')
parser$add_argument("-k", "--keep",      help='Keep discarded values at the end of the file',      action='store_true')
parser$add_argument("-m", "--min",       help='Min RPKM value',                  type="double",    default=0)
parser$add_argument("-o", "--output",    help='Output sorted genelist TSV file', type="character", default="./sorted_genelist.tsv")
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Load and order data
print(paste("Load data from ", args$input, sep=""))
original_data <- read.table(args$input, sep=get_file_type(args$input), header=TRUE, stringsAsFactors=FALSE)
print(paste("Number of rows: ", nrow(original_data), sep=""))

selected_data <- original_data[!rowSums(original_data[,c(7:ncol(original_data))] < args$min) == (ncol(original_data)-6),]
print(paste("Discarded from calculation due to RPKM filtering with min=", args$min, ": ", nrow(original_data)-nrow(selected_data), sep=""))

expression_data = selected_data[, c(7:ncol(selected_data))]

if (args$transform)
{
    print("Log2 transform of expression data")
    if (any(expression_data==0)){
        print(paste("Zero values are replaced by ", 1, sep=""))
        expression_data[expression_data==0] <- 1
    }
    expression_data = log2(expression_data)
}

print("Run hopach")
orderonly <- hopach(expression_data, clusters="none", d=args$dist)

expression_data <- expression_data[orderonly$final$order,]
result_data <- original_data[rownames(expression_data),]

if (args$keep)
{
    print("Append discarded data")
    discarded_data = original_data[!rownames(original_data) %in% rownames(expression_data),]
    result_data <- rbind(result_data, discarded_data)
}

# Export results
write.table(result_data,
            file=args$output,
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

print(paste("Results are saved to ", args$output, sep=""))