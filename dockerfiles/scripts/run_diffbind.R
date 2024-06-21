#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

suppressMessages(library(argparse))
suppressMessages(library(DiffBind))


##########################################################################################
#
# v0.0.14
# - add --color parameter
#
# v0.0.13
# - add --blockfile to set multiple groups for multi-factor analysis
#
# v0.0.12
# - export all graphics to pdf too
# - allow to filter out intervals with low raw read counts
#
# v0.0.11
# - add occupancy based consensus peak selection
#
# v0.0.10
# - suppress scientific notation when exporting to TSV
#
# v0.0.9
# - export not filtered TSV results
#
# v0.0.8
# - supports blocking analyses for DESeq2 and EdgeR
#
# v0.0.7
# - add tryCatch to all optional outputs
#
# v0.0.6
# - filtering by P-value or FDR
#
# v0.0.5
# - add P-value cutoff for reported results
#
# v0.0.4
# - increased default padding for generated heatmaps
#
# v0.0.3
# - allows to control threads number
#
# v0.0.2
# - exports
#   * peak overlap correlation heatmap
#   * counts correlation heatmap
#   * correlation heatmap based on all normalized data
#   * correlation heatmap based on DB sites only
#   * PCA plot using affinity data for only differentially bound sites
#   * MA plot
#   * volcano plot
#   * box plots of read distributions for significantly differentially bound (DB) sites
# - allows to choose from deseq2 or edger
#
# v0.0.1
# - use DiffBind with default parameters
# - use only condition option in comparison
# - export results as TSV
#
##########################################################################################


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_data_set <- function(diff_dba, peaks, reads, name, condition, peakformat, replicate, block) {
    if (is.vector(block)){
        cat(paste("\nLoading data for condition '", condition, "' as '", name, "', replicate = ", replicate, "\n - ", reads, "\n - ", peaks, "\n", sep=""))
        diff_dba <- dba.peakset(DBA=diff_dba, sampID=name, peaks=peaks, bamReads=reads, replicate=replicate, condition=condition, peak.caller=peakformat)
    } else {  # load group from dataframe
        group <- block[name, "group"]
        cat(paste("\nLoading data for condition '", condition, "' as '", name, "', replicate = ", replicate, ", group = '", group, "'\n - ", reads, "\n - ", peaks, "\n", sep=""))
        diff_dba <- dba.peakset(DBA=diff_dba, sampID=name, peaks=peaks, bamReads=reads, replicate=replicate, condition=condition, tissue=group, peak.caller=peakformat)
    }
    return (diff_dba)
}


assert_args <- function(args){
    if (is.null(args$name1) & is.null(args$name2)){
        if ( (length(args$read1) != length(args$peak1)) | (length(args$read2) != length(args$peak2)) ){
            cat("\nNot correct number of inputs provided as -r1, -r2, -p1, -p1")
            quit(save = "no", status = 1, runLast = FALSE)
        } else {
            for (i in 1:length(args$read1)) {
                args$name1 = append(args$name1, head(unlist(strsplit(basename(args$read1[i]), ".", fixed = TRUE)), 1))
            }
            for (i in 1:length(args$read2)) {
                args$name2 = append(args$name2, head(unlist(strsplit(basename(args$read2[i]), ".", fixed = TRUE)), 1))
            }
        }
    } else {
        if ( (length(args$read1) != length(args$peak1)) |
             (length(args$name1) != length(args$peak1)) | 
             (length(args$read2) != length(args$peak2)) |
             (length(args$name2) != length(args$peak2)) ){
            cat("\nNot correct number of inputs provided as -r1, -r2, -p1, -p1, -n1, -n2")
            quit(save = "no", status = 1, runLast = FALSE)
        }
    }
    if (args$method == "deseq2"){
        args$method <- DBA_DESEQ2
    } else if (args$method == "edger") {
        args$method <- DBA_EDGER
    } else if (args$method == "all") {
        args$method <- DBA_ALL_METHODS_BLOCK
    }

    if (args$cparam == "fdr"){
        args$cparam <- FALSE
    } else if (args$cparam == "pvalue") {
        args$cparam <- TRUE
    }

    if (args$usecommon){
        args$minoverlap = 1
    }

    if (!is.null(args$blockfile)){
        block_metadata <- read.table(
            args$blockfile,
            sep=get_file_type(args$blockfile),
            row.names=1,
            col.names=c("name", "group"),
            header=FALSE,
            stringsAsFactors=FALSE
        )
        if (all(is.element(c(args$name1, args$name2), rownames(block_metadata)))){
            args$block <- block_metadata  # dataframe
        } else {
            cat("\nMissing values in block metadata file. Using non-blocked analysis\n")
            args$block = rep(FALSE,length(args$read1)+length(args$read2))  # boolean vector
        }
    } else if (!is.null(args$block)){
        blocked_attributes = as.logical(args$block)
        if (any(is.na(blocked_attributes)) | length(blocked_attributes) != length(args$read1)+length(args$read2) ){
            args$block = is.element(c(args$name1, args$name2), args$block)  # boolean vector
        } else {
            args$block = blocked_attributes  # boolean vector
        }
    } else {
        args$block = rep(FALSE,length(args$read1)+length(args$read2))  # boolean vector
    }

    return (args)
}


export_raw_counts_correlation_heatmap <- function(dba_data, rootname, padding, color_scheme, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotHeatmap(dba_data, margin=padding, score=DBA_SCORE_READS, colScheme=color_scheme)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotHeatmap(dba_data, margin=padding, score=DBA_SCORE_READS, colScheme=color_scheme)
            dev.off()

            cat(paste("\nExport raw counts correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export raw counts correlation heatmap to ", rootname, ".(png/pdf)",  sep=""))
        }
    )
}


export_peak_overlap_correlation_heatmap <- function(dba_data, rootname, padding, color_scheme, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotHeatmap(dba_data, margin=padding, colScheme=color_scheme)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotHeatmap(dba_data, margin=padding, colScheme=color_scheme)
            dev.off()

            cat(paste("\nExport peak overlap correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export peak overlap correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_peak_overlap_rate_plot <- function(peak_overlap_rate, rootname, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            plot(peak_overlap_rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets', pch=19, cex=2)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            plot(peak_overlap_rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets', pch=19, cex=2)
            dev.off()

            cat(paste("\nExport peak overlap rate plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export peak overlap rate plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_normalized_counts_correlation_heatmap <- function(dba_data, rootname, method, padding, color_scheme, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotHeatmap(dba_data, contrast=1, colScheme=color_scheme, th=th, bUsePval=use_pval, method=method, margin=padding)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotHeatmap(dba_data, contrast=1, colScheme=color_scheme, th=th, bUsePval=use_pval, method=method, margin=padding)
            dev.off()

            cat(paste("\nExport normalized counts correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export normalized counts correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_binding_heatmap <- function(dba_data, rootname, method, padding, color_scheme, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {
            
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotHeatmap(dba_data, contrast=1, colScheme=color_scheme, correlations=FALSE, th=th, bUsePval=use_pval, method=method, margin=padding, scale="row")
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotHeatmap(dba_data, contrast=1, colScheme=color_scheme, correlations=FALSE, th=th, bUsePval=use_pval, method=method, margin=padding, scale="row")
            dev.off()

            cat(paste("\nExport binding heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export binding heatmap to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_pca_plot <- function(dba_data, rootname, method, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotPCA(dba_data, attributes=DBA_CONDITION, contrast=1, th=th, bUsePval=use_pval, label=DBA_ID, method=method)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotPCA(dba_data, attributes=DBA_CONDITION, contrast=1, th=th, bUsePval=use_pval, label=DBA_ID, method=method)
            dev.off()

            cat(paste("\nExport PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_ma_plot <- function(dba_data, rootname, method, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {
            
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotMA(dba_data, contrast=1, method=method, th=th, bUsePval=use_pval)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotMA(dba_data, contrast=1, method=method, th=th, bUsePval=use_pval)
            dev.off()

            cat(paste("\nExport MA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export MA plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_volcano_plot <- function(dba_data, rootname, method, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotVolcano(dba_data, contrast=1, method=method, th=th, bUsePval=use_pval)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotVolcano(dba_data, contrast=1, method=method, th=th, bUsePval=use_pval)
            dev.off()

            cat(paste("\nExport volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_box_plot <- function(dba_data, rootname, method, th=1, use_pval=FALSE, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotBox(dba_data, method=method, th=th, bUsePval=use_pval)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotBox(dba_data, method=method, th=th, bUsePval=use_pval)
            dev.off()

            cat(paste("\nExport box plot to", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export box plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_consensus_peak_venn_diagram <- function(dba_data, rootname, width=800, height=800, res=72){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=res)
            dba.plotVenn(dba_data, mask=dba_data$masks$Consensus)
            dev.off()

            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/res), height=round(height/res))
            dba.plotVenn(dba_data, mask=dba_data$masks$Consensus)
            dev.off()

            cat(paste("\nExport consensus peak venn diagram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export consensus peak venn diagram to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_results <- function(dba_data, filename, method, th=1, use_pval=FALSE){
    tryCatch(
        expr = {
            diff_dba.DB <- dba.report(dba_data, contrast=1, DataType=DBA_DATA_FRAME, method=method, bCalled=TRUE, bCounts=TRUE, th=th, bUsePval=use_pval)
            if (!is.null(diff_dba.DB)){
                write.table(diff_dba.DB, file=filename, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
                cat(paste("\nExport differential binding analysis results as TSV to", filename, sep=" "))
            } else {
                cat(paste("\nSkip exporting differential binding analysis results as TSV to", filename, "[no data]", sep=" "))
            }
        },
        error = function(e){
            cat(paste("\nFailed to export differential binding analysis results as TSV to", filename, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Differential binding analysis of ChIP-Seq experiments using affinity (read count) data')
    parser$add_argument("-r1", "--read1",        help='Read files for condition 1. Minimim 2 files in BAM format', type="character", required="True",  nargs='+')
    parser$add_argument("-r2", "--read2",        help='Read files for condition 2. Minimim 2 files in BAM format', type="character", required="True",  nargs='+')
    parser$add_argument("-p1", "--peak1",        help='Peak files for condition 1. Minimim 2 files in format set with -pf', type="character", required="True",  nargs='+')
    parser$add_argument("-p2", "--peak2",        help='Peak files for condition 2. Minimim 2 files in format set with -pf', type="character", required="True",  nargs='+')

    parser$add_argument("-n1", "--name1",        help='Sample names for condition 1. Default: basenames of -r1 without extensions', type="character", nargs='*')
    parser$add_argument("-n2", "--name2",        help='Sample names for condition 2. Default: basenames of -r2 without extensions', type="character", nargs='*')

    parser$add_argument("-bl", "--block",        help='Blocking attribute for multi-factor analysis. Minimum 2. Either names from --name1 or/and --name2 or array of bool based on [read1]+[read2]. Default: not applied', type="character", nargs='*')
    parser$add_argument("-bf", "--blockfile",    help='Blocking attribute metadata file for multi-factor analysis. Headerless TSV/CSV file. First column - names from --name1 and --name2, second column - group name. --block is ignored', type="character")

    parser$add_argument("-pf", "--peakformat",   help='Peak files format. One of [raw, bed, narrow, macs, bayes, tpic, sicer, fp4, swembl, csv, report]. Default: macs', type="character", choices=c("raw","bed","narrow","macs","bayes","tpic","sicer","fp4","swembl","csv","report"), default="macs")

    parser$add_argument("-c1", "--condition1",   help='Condition 1 name, single word with letters and numbers only. Default: condition_1', type="character", default="condition_1")
    parser$add_argument("-c2", "--condition2",   help='Condition 2 name, single word with letters and numbers only. Default: condition_2', type="character", default="condition_2")
    parser$add_argument("-me", "--method",       help='Method by which to analyze differential binding affinity. Default: all', type="character", choices=c("edger","deseq2","all"), default="all")
    parser$add_argument("-mo", "--minoverlap",   help='Min peakset overlap. Only include peaks in at least this many peaksets when generating consensus peakset. Default: 2', type="integer", default=2)
    parser$add_argument("-uc", "--usecommon",    help='Derive consensus peaks only from the common peaks within each condition. Min peakset overlap and min read counts are ignored. Default: false', action='store_true')
    parser$add_argument(
        "--summits",
        help=paste(
            "Width in bp to extend peaks around their summits in both directions",
            "and replace the original ones. Set it to 100 bp for ATAC-Seq and 200",
            "bp for ChIP-Seq datasets. To skip peaks extension and replacement, set",
            "it to negative value. Default: 200 bp (results in 401 bp wide peaks)"
        ),
        type="integer", default=200
    )
    parser$add_argument("-cu", "--cutoff",       help='Cutoff for reported results. Applied to the parameter set with -cp. Default: 0.05', type="double",    default=0.05)
    parser$add_argument("-cp", "--cparam",       help='Parameter to which cutoff should be applied (fdr or pvalue). Default: fdr',         type="character", choices=c("pvalue","fdr"), default="fdr")

    parser$add_argument("-co", "--color",        help='Color scheme. Default: Greens', type="character", choices=c("Reds", "Greens", "Blues", "Greys", "YlOrRd", "Oranges"), default="Greens")
    parser$add_argument("-th", "--threads",      help='Threads to use',                                     type="integer",   default=1)
    parser$add_argument("-pa", "--padding",      help='Padding for generated heatmaps. Default: 20',        type="integer",   default=20)
    parser$add_argument("-o",  "--output",       help='Output prefix. Default: diffbind',                   type="character", default="./diffbind")
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    return (args)
}


args <- get_args()


diff_dba <- NULL
for (i in 1:length(args$read1)){
    diff_dba <- load_data_set(diff_dba, args$peak1[i], args$read1[i], args$name1[i], args$condition1, args$peakformat, i, args$block)
}
for (i in 1:length(args$read2)){
    diff_dba <- load_data_set(diff_dba, args$peak2[i], args$read2[i], args$name2[i], args$condition2, args$peakformat, i, args$block)
}


diff_dba$config$cores <- args$threads


cat("\nLoaded data\n - chrM removed\n - intersected peaks merged\n\n")
print(diff_dba)


mask_cond_1 <- dba.mask(diff_dba, attribute=DBA_CONDITION, value=args$condition1)
mask_cond_2 <- dba.mask(diff_dba, attribute=DBA_CONDITION, value=args$condition2)
mask_cond_both <- as.logical(mask_cond_1 + mask_cond_2)

if (args$usecommon){
    cat("\nDeriving consensus peaks only from the common peaks within each condition. Min peakset overlap is ignored.\n")
    diff_dba_consensus <- dba.peakset(diff_dba, consensus=DBA_CONDITION, minOverlap=0.99)  # use 0.99 to include all intersected peaks within condition
    diff_dba_consensus <- dba(diff_dba_consensus, mask=diff_dba_consensus$masks$Consensus, minOverlap=args$minoverlap)
    cat("\nDerived consensus data\n - intersected peaks merged\n\n")
    print(diff_dba_consensus)
    consensus_peaks <- dba.peakset(diff_dba_consensus, bRetrieve=TRUE, minOverlap=args$minoverlap)
    cat("\nConsensus peaks\n\n")
    print(consensus_peaks)
    export_consensus_peak_venn_diagram(diff_dba_consensus, paste(args$output, "_consensus_peak_venn_diagram", sep=""))
}


# Get peak overlap rates, export plots
cat("\n")
cat(paste("\nPeak overlap rate for", args$condition1, "\n", sep=" "))
peak_overlap_rate_cond_1 = dba.overlap(diff_dba, mask=mask_cond_1, mode=DBA_OLAP_RATE)
print(peak_overlap_rate_cond_1)

cat(paste("\nPeak overlap rate for", args$condition2, "\n", sep=" "))
peak_overlap_rate_cond_2 = dba.overlap(diff_dba, mask=mask_cond_2, mode=DBA_OLAP_RATE)
print(peak_overlap_rate_cond_2)

export_peak_overlap_rate_plot(peak_overlap_rate_cond_1, paste(args$output, "_condition_1_peak_overlap_rate", sep=""))
export_peak_overlap_rate_plot(peak_overlap_rate_cond_2, paste(args$output, "_condition_2_peak_overlap_rate", sep=""))

if (!args$usecommon){
    peak_overlap_rate_all = dba.overlap(diff_dba, mode=DBA_OLAP_RATE)
    cat("\n\nPeak overlap rate for all peaksets\n")
    print(peak_overlap_rate_all)
    export_peak_overlap_rate_plot(peak_overlap_rate_all, paste(args$output, "_all_peak_overlap_rate", sep=""))
}


# Export peak overlap correlation heatmap
export_peak_overlap_correlation_heatmap(diff_dba, paste(args$output, "_peak_overlap_correlation_heatmap", sep=""), args$padding, args$color)


# Count reads in binding site intervals
cat(paste("\n\nCount reads using", diff_dba$config$cores, "threads. Min peakset overlap is set to", args$minoverlap, "\n", sep=" "))

if (args$usecommon){
    diff_dba <- dba.count(
        diff_dba,
        peaks=consensus_peaks,
        summits=ifelse(args$summits >= 0, args$summits, FALSE),
        minOverlap=args$minoverlap
    )
} else {
    diff_dba <- dba.count(
        diff_dba,
        summits=ifelse(args$summits >= 0, args$summits, FALSE),
        minOverlap=args$minoverlap
    )
}

cat("\nCounted data\n\n")
print(diff_dba)


# Export raw counts correlation heatmap
export_raw_counts_correlation_heatmap(diff_dba, paste(args$output, "_raw_counts_correlation_heatmap", sep=""), args$padding, args$color)


# Export consensus peak venn diagram
if (!args$usecommon){
    export_consensus_peak_venn_diagram(diff_dba, paste(args$output, "_consensus_peak_venn_diagram", sep=""))
}


# Establish contrast
metadata = dba.show(diff_dba)
diff_dba$contrasts <- NULL
cat(paste("\n\nEstablish contrast [", paste(metadata[metadata["Condition"]==args$condition1, "ID"], collapse=", "), "] vs [", paste(metadata[metadata["Condition"]==args$condition2, "ID"], collapse=", "), "]", "\n\n", sep=""))
if (is.vector(args$block)){
    cat(paste("Blocked attributes [", paste(metadata[args$block, "ID"], collapse=", "), "]", "\n\n", sep=""))
    diff_dba <- dba.contrast(diff_dba, group1=mask_cond_1, group2=mask_cond_2, name1=args$condition1, name2=args$condition2, block=args$block)
} else {
    cat("Blocked attributes\n")
    print(args$block)
    cat("\n")
    diff_dba <- dba.contrast(diff_dba, group1=mask_cond_1, group2=mask_cond_2, name1=args$condition1, name2=args$condition2, block=DBA_TISSUE)
}


diff_dba <- dba.analyze(diff_dba, method=args$method)


cat("\nAnalyzed data\n\n")
print(diff_dba)


# Export all normalized counts correlation heatmaps
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_all_normalized_counts_correlation_heatmap_deseq", sep=""),
                                             DBA_DESEQ2,
                                             args$padding,
                                             args$color)
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_all_normalized_counts_correlation_heatmap_deseq_block", sep=""),
                                             DBA_DESEQ2_BLOCK,
                                             args$padding,
                                             args$color)
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_all_normalized_counts_correlation_heatmap_edger", sep=""),
                                             DBA_EDGER,
                                             args$padding,
                                             args$color)
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_all_normalized_counts_correlation_heatmap_edger_block", sep=""),
                                             DBA_EDGER_BLOCK,
                                             args$padding,
                                             args$color)


# Export filtered normalized counts correlation heatmaps
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_filtered_normalized_counts_correlation_heatmap_deseq", sep=""),
                                             DBA_DESEQ2,
                                             args$padding,
                                             args$color,
                                             args$cutoff,
                                             args$cparam)
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_filtered_normalized_counts_correlation_heatmap_deseq_block", sep=""),
                                             DBA_DESEQ2_BLOCK,
                                             args$padding,
                                             args$color,
                                             args$cutoff,
                                             args$cparam)                                             
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_filtered_normalized_counts_correlation_heatmap_edger", sep=""),
                                             DBA_EDGER,
                                             args$padding,
                                             args$color,
                                             args$cutoff,
                                             args$cparam)
export_normalized_counts_correlation_heatmap(diff_dba,
                                             paste(args$output, "_filtered_normalized_counts_correlation_heatmap_edger_block", sep=""),
                                             DBA_EDGER_BLOCK,
                                             args$padding,
                                             args$color,
                                             args$cutoff,
                                             args$cparam)


# Export filtered binding heatmaps
export_binding_heatmap(diff_dba,
                       paste(args$output, "_filtered_binding_heatmap_deseq", sep=""),
                       DBA_DESEQ2,
                       args$padding,
                       args$color,
                       args$cutoff,
                       args$cparam)
export_binding_heatmap(diff_dba,
                       paste(args$output, "_filtered_binding_heatmap_deseq_block", sep=""),
                       DBA_DESEQ2_BLOCK,
                       args$padding,
                       args$color,
                       args$cutoff,
                       args$cparam)
export_binding_heatmap(diff_dba,
                       paste(args$output, "_filtered_binding_heatmap_edger", sep=""),
                       DBA_EDGER,
                       args$padding,
                       args$color,
                       args$cutoff,
                       args$cparam)
export_binding_heatmap(diff_dba,
                       paste(args$output, "_filtered_binding_heatmap_edger_block", sep=""),
                       DBA_EDGER_BLOCK,
                       args$padding,
                       args$color,
                       args$cutoff,
                       args$cparam)


# Export not filtered PCA plots
export_pca_plot(diff_dba,
                paste(args$output, "_all_pca_plot_deseq", sep=""),
                DBA_DESEQ2)
export_pca_plot(diff_dba,
                paste(args$output, "_all_pca_plot_deseq_block", sep=""),
                DBA_DESEQ2_BLOCK)
export_pca_plot(diff_dba,
                paste(args$output, "_all_pca_plot_edger", sep=""),
                DBA_EDGER)
export_pca_plot(diff_dba,
                paste(args$output, "_all_pca_plot_edger_block", sep=""),
                DBA_EDGER_BLOCK)


# Export filtered PCA plots
export_pca_plot(diff_dba,
                paste(args$output, "_filtered_pca_plot_deseq", sep=""),
                DBA_DESEQ2,
                args$cutoff,
                args$cparam)
export_pca_plot(diff_dba,
                paste(args$output, "_filtered_pca_plot_deseq_block", sep=""),
                DBA_DESEQ2_BLOCK,
                args$cutoff,
                args$cparam)
export_pca_plot(diff_dba,
                paste(args$output, "_filtered_pca_plot_edger", sep=""),
                DBA_EDGER,
                args$cutoff,
                args$cparam)
export_pca_plot(diff_dba,
                paste(args$output, "_filtered_pca_plot_edger_block", sep=""),
                DBA_EDGER_BLOCK,
                args$cutoff,
                args$cparam)


# Export filtered MA plot
export_ma_plot(diff_dba,
               paste(args$output, "_filtered_ma_plot_deseq", sep=""),
               DBA_DESEQ2,
               args$cutoff,
               args$cparam)
export_ma_plot(diff_dba,
               paste(args$output, "_filtered_ma_plot_deseq_block", sep=""),
               DBA_DESEQ2_BLOCK,
               args$cutoff,
               args$cparam)
export_ma_plot(diff_dba,
               paste(args$output, "_filtered_ma_plot_edger", sep=""),
               DBA_EDGER,
               args$cutoff,
               args$cparam)
export_ma_plot(diff_dba,
               paste(args$output, "_filtered_ma_plot_edger_block", sep=""),
               DBA_EDGER_BLOCK,
               args$cutoff,
               args$cparam)


# Export filtered volcano plots
export_volcano_plot(diff_dba,
                    paste(args$output, "_filtered_volcano_plot_deseq", sep=""),
                    DBA_DESEQ2,
                    args$cutoff,
                    args$cparam)
export_volcano_plot(diff_dba,
                    paste(args$output, "_filtered_volcano_plot_deseq_block", sep=""),
                    DBA_DESEQ2_BLOCK,
                    args$cutoff,
                    args$cparam)
export_volcano_plot(diff_dba,
                    paste(args$output, "_filtered_volcano_plot_edger", sep=""),
                    DBA_EDGER,
                    args$cutoff,
                    args$cparam)
export_volcano_plot(diff_dba,
                    paste(args$output, "_filtered_volcano_plot_edger_block", sep=""),
                    DBA_EDGER_BLOCK,
                    args$cutoff,
                    args$cparam)


# Export filtered box plots
export_box_plot(diff_dba,
                paste(args$output, "_filtered_box_plot_deseq", sep=""),
                DBA_DESEQ2,
                args$cutoff,
                args$cparam)
export_box_plot(diff_dba,
                paste(args$output, "_filtered_box_plot_deseq_block", sep=""),
                DBA_DESEQ2_BLOCK,
                args$cutoff,
                args$cparam)
export_box_plot(diff_dba,
                paste(args$output, "_filtered_box_plot_edger", sep=""),
                DBA_EDGER,
                args$cutoff,
                args$cparam)
export_box_plot(diff_dba,
                paste(args$output, "_filtered_box_plot_edger_block", sep=""),
                DBA_EDGER_BLOCK,
                args$cutoff,
                args$cparam)


# Export filtered results
export_results(diff_dba,
               paste(args$output, "_filtered_report_deseq.tsv", sep=""),
               DBA_DESEQ2,
               args$cutoff,
               args$cparam)
export_results(diff_dba,
               paste(args$output, "_filtered_report_deseq_block.tsv", sep=""),
               DBA_DESEQ2_BLOCK,
               args$cutoff,
               args$cparam)
export_results(diff_dba,
               paste(args$output, "_filtered_report_edger.tsv", sep=""),
               DBA_EDGER,
               args$cutoff,
               args$cparam)
export_results(diff_dba,
               paste(args$output, "_filtered_report_edger_block.tsv", sep=""),
               DBA_EDGER_BLOCK,
               args$cutoff,
               args$cparam)


# Export not filtered results
export_results(diff_dba,
               paste(args$output, "_all_report_deseq.tsv", sep=""),
               DBA_DESEQ2)
export_results(diff_dba,
               paste(args$output, "_all_report_deseq_block.tsv", sep=""),
               DBA_DESEQ2_BLOCK)
export_results(diff_dba,
               paste(args$output, "_all_report_edger.tsv", sep=""),
               DBA_EDGER)
export_results(diff_dba,
               paste(args$output, "_all_report_edger_block.tsv", sep=""),
               DBA_EDGER_BLOCK)