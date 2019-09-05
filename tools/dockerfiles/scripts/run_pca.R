#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
suppressMessages(library(argparse))
suppressMessages(library(scatterplot3d))


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
##########################################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = ","
    if (ext == "tsv"){
        separator = "\t"
    }
    return (separator)
}


load_data_set <- function(filenames, prefixes, target_colname, intersect_by) {
    selected_data <- NULL
    for (i in 1:length(filenames)) {
        raw_data <- read.table(filenames[i], sep=get_file_type(filenames[i]), header=TRUE, stringsAsFactors=FALSE)
        print(paste("Load ", nrow(raw_data), " rows from ", filenames[i], sep=""))
        colnames(raw_data)[colnames(raw_data) == target_colname] <- prefixes[i]
        if (is.null(selected_data)){
            selected_data <- raw_data
        } else {
            if (is.null(intersect_by)){
                selected_data <- cbind(selected_data, raw_data)
            } else {
                selected_data <- merge(selected_data, raw_data, by=intersect_by, sort = FALSE)
            }
        }
    }
    return (selected_data[,prefixes])
}


# Parser
parser <- ArgumentParser(description='Run BioWardrobe PCA')
parser$add_argument("-i", "--input",     help='Input CSV/TSV files',                   type="character", required="True", nargs='+')
parser$add_argument("-n", "--name",      help='Input aliases, the order corresponds to --input order. Default: basename of --input files', type="character", nargs='+')
parser$add_argument("-t", "--target",    help='Target column name to be used by PCA',  type="character", default="Rpkm")
parser$add_argument("-c", "--combine",   help='Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand', type="character", nargs='+', default=c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand"))
parser$add_argument("-o", "--output",    help='Output prefix. Default: pca_',          type="character", default="./pca_")
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Set default value for --name if it wasn't provided
if(is.null(args$name)){
    for (i in 1:length(args$input)) {
        args$name = append(args$name, head(unlist(strsplit(basename(args$input[i]), ".", fixed = TRUE)), 1))
    }
}

png(filename=paste(args$output, "%03d.png", sep=""), width=800, height=800)
target_data <- load_data_set(args$input, args$name, args$target, args$combine)
filtered_target_data <- target_data[rowSums(target_data) != 0,]

icolor <- colorRampPalette(c("red", "black", "green", "yellow", "blue", "pink", "brown"))(length(colnames(filtered_target_data)))

pca <- prcomp(t(filtered_target_data), cor=TRUE, scale.=T)
result <- pca$x

result_df <- as.data.frame(result)
result_df <- cbind(exp=rownames(result_df), result_df)

write.table(result_df,
            file=paste(args$output, "pca.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

plot(result[,1], result[,2], col=icolor, xlab="PCA1", ylab="PCA2", main="")
legend("bottomright", text.col=icolor, bg="white", legend = args$name, yjust=0, horiz=F, bty='n', cex=0.8)

plot(result[,2], result[,3], col=icolor, xlab="PCA2", ylab="PCA3", main="")
legend("bottomright", text.col=icolor, bg="white", legend = args$name, yjust=0, horiz=F, bty='n', cex=0.8)

plot(pca, type="lines")

s3d <- scatterplot3d(result[,1], result[,2], result[,3], xlab="PC1", ylab="PC2", zlab="PC3", main="",
	                 color=icolor, col.axis="blue", sub="", box=T, lwd=5, type="p")

legend("bottomright", inset=c(0.03,0.03), text.col=icolor, bg="white", legend=args$name, yjust=0, horiz=F, bty='n', cex=0.8)

graphics.off()
