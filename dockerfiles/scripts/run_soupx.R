#!/usr/bin/env Rscript
options("warn"=-1)
options("width"=300)


suppressMessages(library(argparse))
suppressMessages(library(SoupX))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(viridis))
suppressMessages(library(gridExtra))
suppressMessages(library(DropletUtils))


load_target_genes <- function(filename) {
    target_genes_vector <- scan(filename, character())
    return (target_genes_vector)
}


export_gene_expression_plots <- function(barcode_matrix, metadata, rootname, target_genes, threshold, width=1200, height=1200, resolution=72){
    tryCatch(
        expr = {
            plot_data = metadata[colnames(barcode_matrix),]
            mids = aggregate(cbind(tSNE1, tSNE2) ~ clusters, data=plot_data, FUN=mean)
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            plots = list()
            for(i in 1:length(target_genes)) {
                target_gene = target_genes[i]
                plot_data[,target_gene] = barcode_matrix[target_gene,] > threshold
                plots[[i]] = ggplot(plot_data, aes(tSNE1, tSNE2, color=as.factor(clusters))) +
                    ggtitle(target_gene) +
                    geom_point(size=3, shape=19, aes_string(alpha=target_gene)) +
                    xlab("t-SNE-1") +
                    ylab("t-SNE-2") +
                    geom_label_repel(data=mids, aes(label=clusters), point.padding=0.5, box.padding=0.5, check_overlap = TRUE, size = 10, show.legend = FALSE) +
                    scale_color_viridis(discrete=TRUE, name="Clusters") +
                    guides(alpha=guide_legend(title=paste("expression >", threshold, sep=" ")))
            }
            do.call(grid.arrange, plots)
            dev.off()

            cat(paste("\nExport gene expression plots to ", rootname, ".pdf", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export gene expression plots to ", rootname, ".pdf", "\n",  sep=""))
        }
    )
}


export_ratio_plots <- function(soup_data, rootname, target_genes, fdr, width=1200, height=1200, resolution=72){
    tryCatch(
        expr = {
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            plots = list()
            for(i in 1:length(target_genes)) {
                target_gene = target_genes[i]
                plots[[i]] = plotMarkerMap(soup_data, target_gene, FDR=fdr) + 
                ggtitle(target_gene) +
                guides(colour=guide_legend(title=paste("FDR <", fdr, sep=" ")))
            }
            do.call(grid.arrange, plots)
            dev.off()
            cat(paste("\nExport ratio plots to ", rootname, ".pdf", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export ratio plots to ", rootname, ".pdf", "\n",  sep=""))
        }
    )
}


export_soup_fraction_plots <- function(soup_data, corrected_counts, rootname, target_genes, width=1200, height=1200, resolution=72){
    tryCatch(
        expr = {
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            plots = list()
            for(i in 1:length(target_genes)) {
                target_gene = target_genes[i]
                plots[[i]] = plotChangeMap(soup_data, corrected_counts, target_gene) + 
                ggtitle(target_gene)
            }
            do.call(grid.arrange, plots)
            dev.off()
            cat(paste("\nExport soup fraction plots to ", rootname, ".pdf", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export soup fraction plots to ", rootname, ".pdf", "\n",  sep=""))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="Runs SoupX - an R package for the estimation and removal of cell free mRNA contamination")
    parser$add_argument("--counts",    help="Path to the output folder produced by 'cellranger count' command", type="character", required="True")
    parser$add_argument("--genes",     help="Path to the file with target genes (headerless, one gene per line)", type="character")
    parser$add_argument("--threshold", help="Expression threshold for displaying target genes on a plot (expression > threshold). Default: 0.0", type="double", default=0.0)
    parser$add_argument("--fdr",       help="FDR cutoff for expression ratio plots. Default: 0.05", type="double", default=0.05)
    parser$add_argument("--round",     help="Round adjusted counts to integers", action="store_true", default = FALSE)
    parser$add_argument("--format",    help='Output matrix format version. Default: 3 (corresponds to the latest Cell Ranger matrix format)', type="character", choices=c("2","3"), default="3")
    parser$add_argument("--output",    help='Output prefix. Default: soupx', type="character", default="./soupx")
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()


soup_data = load10X(args$counts)
pdf(file=paste(args$output, "contamination_estimation_plot.pdf", sep="_"))
soup_data = autoEstCont(soup_data)
dev.off()
adjusted_soup_counts = adjustCounts(soup_data, roundToInt=args$round)
DropletUtils:::write10xCounts(
    path=paste(args$output, "adjusted_counts", sep="_"),
    x=adjusted_soup_counts,
    type="sparse",
    version=args$format,
    overwrite=TRUE
)
DropletUtils:::write10xCounts(
    path=paste(args$output, "adjusted_counts.h5", sep="_"),
    x=adjusted_soup_counts,
    type="HDF5",
    version=args$format,
    overwrite=TRUE
)


if (!is.null(args$genes)){
    target_genes = load_target_genes(args$genes)
    export_gene_expression_plots(
        soup_data$toc,
        soup_data$metaData,
        paste(args$output, "raw_gene_expression_plots", sep="_"),
        target_genes,
        args$threshold
    )
    export_gene_expression_plots(
        adjusted_soup_counts,
        soup_data$metaData,
        paste(args$output, "adjusted_gene_expression_plots", sep="_"),
        target_genes,
        args$threshold
    )
    export_ratio_plots(
        soup_data,
        paste(args$output, "raw_gene_expression_to_pure_soup_ratio_plots", sep="_"),
        target_genes,
        args$fdr
    )
    export_soup_fraction_plots(
        soup_data,
        adjusted_soup_counts,
        paste(args$output, "raw_to_adjusted_gene_expression_ratio_plots", sep="_"),
        target_genes
    )
}