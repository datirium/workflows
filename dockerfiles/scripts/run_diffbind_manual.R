#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(Glimma))
suppressMessages(library(hopach))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(forcats))
suppressMessages(library(argparse))
suppressMessages(library(DiffBind))
suppressMessages(library(tidyverse))
suppressMessages(library(htmlwidgets))
suppressMessages(library(BiocParallel))
suppressMessages(library(GenomicRanges))
suppressMessages(library(EnhancedVolcano))


theme_set(theme_classic())                       # set classic theme for all ggplot generated graphics


D40_COLORS <- c("#FF6E6A", "#71E869", "#6574FF", "#F3A6B5", "#FF5AD6", "#6DDCFE", "#FFBB70", "#43A14E", "#D71C7C", "#E1E333", "#8139A8", "#00D8B6", "#B55C00", "#7FA4B6", "#FFA4E3", "#B300FF", "#9BC4FD", "#FF7E6A", "#9DE98D", "#BFA178", "#E7C2FD", "#8B437D", "#ADCDC0", "#FE9FA4", "#FF53D1", "#D993F9", "#FF47A1", "#FFC171", "#625C51", "#4288C9", "#9767D4", "#F2D61D", "#8EE6FD", "#B940B1", "#B2D5F8", "#9AB317", "#C70000", "#AC8BAC", "#D7D1E4", "#9D8D87")
ALL_CRITERIA <- c("Tissue", "Factor", "Condition", "Treatment", "Replicate")


get_file_type <- function (filename){
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


get_metadata <- function(args){
    dba_metadata <- read.table(
        args$metadata,
        sep=get_file_type(args$metadata),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )  %>% dplyr::rename("SampleID"=1) %>% mutate_at(colnames(.)[2:length(colnames(.))], factor)

    if (!is.null(args$base)){
        print(
            paste(
                "Attempting to relevel metadata table based on",
                paste(args$base, collapse=", "), "values"
            )
        )
        for (i in 1:length(args$base)) {
            current_base_level <- args$base[i]
            tryCatch(
                expr = {
                    current_column <- colnames(dba_metadata)[i+1]
                    dba_metadata[[current_column]] <- relevel(dba_metadata[[current_column]], current_base_level)
                    print(paste("Setting", current_base_level, "as a base level for", current_column))
                },
                error = function(e){
                    print(
                        paste(
                            "Failed to set", current_base_level, "as a base level",
                            "due to", e
                        )
                    )
                }
            )
        }
    }

    raw_metadata <- dba_metadata %>%                     # need it for an easier access
                    remove_rownames() %>%
                    column_to_rownames("SampleID")
    print("Raw metadata")
    print(raw_metadata)

    all_metadata_vars <- unique(colnames(raw_metadata))
    if(length(intersect(all_metadata_vars, ALL_CRITERIA)) != length(all_metadata_vars)){
        print("Metadata includes not supported column names")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    locations_df <- data.frame(
        SampleID=args$aliases,
        bamReads=args$alignments,
        Peaks=args$peaks,
        check.names=FALSE
    )
    
    if (!all(sort(dba_metadata$SampleID) == sort(unique(locations_df$SampleID)))){
        print("Metadata file is malformed. Exiting.")
        quit(save="no", status=1, runLast=FALSE)
    }
    dba_metadata <- dba_metadata %>% dplyr::left_join(locations_df, by="SampleID")
    print("DBA metadata")
    print(dba_metadata)

    return (list(raw=raw_metadata, dba=dba_metadata))
}


get_contrast <- function(design_formula, metadata, args){
    model <- model.matrix(design_formula, metadata)
    print("Model matrix")
    print(model)
    for (i in 1:length(colnames(metadata))){
        current_column <- colnames(metadata)[i]
        unique_keys <- unique(metadata[[current_column]])
        for (j in 1:length(unique_keys)){
            key <- as.character(unique_keys[j])
            print(paste("Subsetting model matrix for", current_column, "column", "with value", key))
            subset <- colMeans(model[metadata[[current_column]] == key, ])
            print(subset)
            assign(key, subset)
        }
    }
    print(paste("Evaluating contrast", args$contrast))
    contrast <- eval(parse(text=args$contrast))
    print(contrast)
    return (contrast)
}


get_clustered_data <- function(data, center, dist, transpose) {

    if (transpose){
        print("Transposing expression data")
        data = t(data)
    }
    if (!is.null(center)) {
        print(paste("Centering expression data by ", center, sep=""))
        if (center == "mean"){
            data = data - rowMeans(data)    
        } else {
            data = data - rowMedians(data.matrix(data))    
        }
    }
    print("Creating distance matrix")
    distance_matrix <- distancematrix(data, dist)
    print("Running HOPACH")
    hopach_results <- hopach(data, dmat=distance_matrix)

    if (transpose){
        print("Transposing expression data")
        data = t(data)
    }

    print("Parsing cluster labels")
    options(scipen=999)                                           # need to temporary disable scientific notation, because nchar gives wrong answer
    clusters = as.data.frame(hopach_results$clustering$labels)
    colnames(clusters) = "label"
    clusters = cbind(
        clusters,
        "HCL"=outer(
            clusters$label,
            10^c((nchar(trunc(clusters$label))[1]-1):0),
            function(a, b) {
                paste0("c", a %/% b)
            }
        )
    )
    clusters = clusters[, c(-1), drop=FALSE]
    options(scipen=0)                                            # setting back to the default value

    return (
        list(
            order=as.vector(hopach_results$clustering$order),
            counts=data,
            clusters=clusters
        )
    )
}


export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    tryCatch(
        expr = {
            row_metadata <- row_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            col_metadata <- col_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            print(paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export GCT data to", location, sep=" "))
        }
    )
}


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE, digits=NULL){
    tryCatch(
        expr = {
            if (!is.null(digits)){
                data <- format(data, digits=digits)
            }
            data <- data %>%                                             # to remove leading and trailing spaces when exporting to TSV
                    mutate(across(everything(), as.character)) %>%
                    mutate(across(everything(), str_trim))
            write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


export_overlap_plot <- function(data, rootname, plot_title, plot_subtitle, highlight=NULL, pdf=FALSE, width=600, height=600, resolution=100){
    tryCatch(
        expr = {
            colors <- rep("darkseagreen", length(data))
            if(!is.null(highlight) && highlight > 0 && highlight <= length(data)){
                colors[highlight] <- "darkseagreen4"
            }
            peaksets <- factor(
                as.character(1:length(data)),                   # need to have it as character, otherwise some labels will be skipped,
                levels=as.character(1:length(data))             # but the levels shouldn't be alphabetical, so we force the proper order
            )
            data_df <- data.frame(
                overlaps=data,
                peaksets=peaksets,
                color=colors,
                check.names=FALSE
            )
            vjust <- 1.6
            hjust <- 0.5
            angle <- 0
            if(length(data) > 5){
                vjust <- 0.5
                hjust <- -0.2
                angle <- 90
            }
            plot <- ggplot(data_df, aes(x=peaksets, y=overlaps, fill=peaksets)) +
                    geom_bar(stat="identity", show.legend=FALSE) +
                    scale_fill_manual(values=data_df$color) +
                    geom_text(
                        aes(label=scales::comma(overlaps)),
                        vjust=vjust,
                        hjust=hjust,
                        angle=angle,
                        y=0,
                        color="darkblue", size=3.5
                    ) +
                    xlab("Min peakset overlap") +
                    ylab("Number of consensus peaks") +
                    ggtitle(plot_title, subtitle=plot_subtitle)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Exporting overlap plot to ", rootname, ".(png/pdf)", sep=""))

        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(
                paste0(
                    "Failed to export overlap plot to ", rootname,
                    ".(png/pdf) with error - ", e
                )
            )
        }
    )
}


export_corr_heatmap <- function(data, rootname, plot_title, score=NULL, padj=1, colors="Greens", pdf=FALSE, margin=20, width=800, height=800, resolution=100){
    tryCatch(
        expr = {

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            if (!is.null(score)){
                dba.plotHeatmap(dba_data, correlations=TRUE, th=padj, score=score, main=plot_title, colScheme=colors, margin=margin, density.info="none")
            } else {
                print("Correlation values")
                print(dba.plotHeatmap(dba_data, correlations=TRUE, breaks=seq(0, 1, length.out=257), th=padj, main=plot_title, colScheme=colors, margin=margin, density.info="none"))
            }
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                if (!is.null(score)){
                    dba.plotHeatmap(dba_data, correlations=TRUE, th=padj, score=score, main=plot_title, colScheme=colors, margin=margin, density.info="none")
                } else {
                    dba.plotHeatmap(dba_data, correlations=TRUE, breaks=seq(0, 1, length.out=257), th=padj, main=plot_title, colScheme=colors, margin=margin, density.info="none")
                }
                dev.off()
            }

            print(paste("Exporting correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(
                paste0(
                    "Failed to export correlation heatmap to ", rootname,
                    ".(png/pdf) with error - ", e
                )
            )
        }
    )
}


export_profile_heatmap <- function(data, rootname, sites, samples=NULL, scores="rank", abs_scores=FALSE, decreasing=TRUE, split_by=NULL, max_sites=1000, colors=c("black", "yellow"), pdf=FALSE, width=800, height=800, resolution=100){
    tryCatch(
        expr = {

            profiles_data <- dba.plotProfile(
                data,
                sites=sites,
                samples=samples,                      # used only when computing profile, not building the plot
                scores=scores,
                absScores=abs_scores,
                maxSites=max_sites,
                matrices_color=colors,
                decreasing=decreasing,
                merge=NULL
            )

            if(!is.null(split_by)){                   # used only to fix a bug with not proper Binding Sites groups order and names
                rowData(profiles_data)$"Binding Sites" <- factor(
                    as.vector(as.character(rowData(profiles_data)[[split_by]])),
                    levels=names(sites)               # to make sure the order of clusters is correct, because sites now is GRangesList split by cluster
                )
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            dba.plotProfile(
                profiles_data,
                scores=scores,
                absScores=abs_scores,
                maxSites=max_sites,
                matrices_color=colors,
                decreasing=decreasing,
                merge=NULL
            )
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                dba.plotProfile(
                    profiles_data,
                    scores=scores,
                    absScores=abs_scores,
                    maxSites=max_sites,
                    matrices_color=colors,
                    decreasing=decreasing,
                    merge=NULL
                )
                dev.off()
            }

            print(paste("Exporting profile heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(
                paste0(
                    "Failed to export profile heatmap to ", rootname,
                    ".(png/pdf) with error - ", e
                )
            )
        }
    )
}


export_ma_plot <- function(data, rootname, fdr_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, pdf=FALSE, width=800, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ggmaplot(data,
                        fdr=fdr_cutoff,
                        fc=y_cutoff,
                        top=0,
                        size=2,
                        alpha=0.6,
                        main=plot_title,
                        subtitle=plot_subtitle,
                        caption=caption,
                        xlab=x_label,
                        ylab=y_label,
                        palette=c("red", "red", "grey30")
                    ) +
                    theme_classic() +
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Exporting MA-plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(
                paste0(
                    "Failed to export MA-plot to ", rootname,
                    ".(png/pdf) with error - ", e
                )
            )
        }
    )
}


export_volcano_plot <- function(data, rootname, x_axis, y_axis, x_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, features=NA, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- EnhancedVolcano(
                        data,
                        x=x_axis,
                        y=y_axis,
                        lab=rownames(data),
                        selectLab=features,
                        FCcutoff=x_cutoff,
                        pCutoff=y_cutoff,
                        xlab=x_label,
                        ylab=y_label,
                        title=plot_title,
                        subtitle=plot_subtitle,
                        caption=caption,
                        labSize=4,
                        labFace="bold",
                        labCol="red4",
                        colAlpha=0.6,
                        col=c("grey30", "grey30", "grey30", "red"),
                        drawConnectors=TRUE,
                        widthConnectors=0.75
                    ) +
                    scale_y_log10() + annotation_logticks(sides="l", alpha=0.3) +
                    theme_classic() +
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(paste("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


export_mds_html_plot <- function(data, groups, location){
    tryCatch(
        expr = {
            saveWidget(
                glimmaMDS(
                    x=data,
                    groups=groups,
                    labels=rownames(groups)
                ),
                file=location
            )
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}


get_pca_data <- function(counts_data){
    tryCatch(
        expr = {
            print("Computing PCA for counts data")
            counts_data <- counts_data %>%
                           filter_all(any_vars(. != 0))   # remove rows with only zeros, otherwise prcomp may fail
            pca_raw <- prcomp(
                t(counts_data),
                center=TRUE,
                scale.=TRUE
            )
            pca_scores <- as.data.frame(pca_raw$x) %>%
                          rownames_to_column(var="group")
            pca_variance <- round(pca_raw$sdev / sum(pca_raw$sdev) * 100, 2)
            return (list(scores=pca_scores, variance=pca_variance))
        },
        error = function(e){
            print(paste("Failed to compute PCA for counts data due to", e))
        }
    )
}


export_pca_plot <- function(data, pcs, rootname, plot_title, legend_title, color_by="label", label_size=5, pt_size=8, pt_shape=19, alpha=0.75, palette_colors=D40_COLORS, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            x_score_column <- paste0("PC", pcs[1])
            y_score_column <- paste0("PC", pcs[2])
            x_variance <- data$variance[pcs[1]]
            y_variance <- data$variance[pcs[2]]
            plot <- ggplot(
                        data$scores,
                        aes_string(x=x_score_column, y=y_score_column, color=color_by)
                    ) +
                    geom_point(size=pt_size, shape=pt_shape, alpha=alpha) +
                    xlab(paste0(x_score_column, ": ", x_variance, "% variance")) +
                    ylab(paste0(y_score_column, ": ", y_variance, "% variance")) + 
                    geom_label_repel(
                        aes_string(label=color_by),
                        size=label_size,
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title)) +
                    scale_color_manual(values=palette_colors) +

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Exporting PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(paste("Failed to export PCA plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}


export_plots <- function(data, metadata, args){
    export_volcano_plot(
        data=data,
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=0,
        y_cutoff=args$padj,
        x_label=bquote(~Log[2]~"fold change"),
        y_label=bquote("-"~Log[10]~"padj"),
        plot_title="Volcano plot for differentially bound sites",
        plot_subtitle=paste0("Differentially bound sites with padj", "<=", args$padj),
        caption=paste(nrow(data), "features"),
        rootname=paste(args$output, "diff_vlcn", sep="_"),
        pdf=args$pdf
    )

    export_ma_plot(
        data=data,
        fdr_cutoff=args$padj,
        y_cutoff=0,
        x_label=bquote(~Log[2]~"base mean"),
        y_label=bquote(~Log[2]~"fold change"),
        plot_title="MA-plot for differentially bound sites",
        plot_subtitle=paste0("Differentially bound sites with padj", "<=", args$padj),
        caption=paste(nrow(data), "features"),
        rootname=paste(args$output, "diff_ma", sep="_"),
        pdf=args$pdf
    )

    pca_data <- get_pca_data(data[, rownames(metadata)])               # adds 'group' column to identify the datasets

    export_pca_plot(
        data=pca_data,
        pcs=c(1, 2),
        plot_title="PCA (1,2) of not filtered normalized counts",
        legend_title="Dataset",
        color_by="group",
        width=ifelse(length(args$aliases) > 13, 1600, 1200),       # need to make it wider if we have a lot of samples
        rootname=paste(args$output, "nr_rds_pca_1_2", sep="_"),
        pdf=args$pdf
    )

    export_pca_plot(
        data=pca_data,
        pcs=c(2, 3),
        plot_title="PCA (2,3) of not filtered normalized counts",
        legend_title="Dataset",
        color_by="group",
        width=ifelse(length(args$aliases) > 13, 1600, 1200),       # need to make it wider if we have a lot of samples
        rootname=paste(args$output, "nr_rds_pca_2_3", sep="_"),
        pdf=args$pdf
    )

    export_mds_html_plot(
        data=as.matrix(data[, rownames(metadata)]),
        groups=metadata,
        location=paste(args$output, "nr_rds_mds.html", sep="_")
    )
}


assert_args <- function(args){
    if (length(args$alignments) != length(args$peaks) || length(args$alignments) != length(args$aliases)){
        print("--alignments, --peaks and/or --aliases parameters have different number of values")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    if(args$method == "edger" && is.null(args$contrast)){
        print("EdgeR can't be used without --contrast. Select DESeq2 or provide contrast.")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    args$method <- switch(
        args$method,
        "deseq2" = DBA_DESEQ2,
        "edger"  = DBA_EDGER
    )

    args$scoreby <- switch(
        args$scoreby,
        "pvalue" = 7,           # -log10(pvalue)
        "qvalue" = 9            # -log10(qvalue)
    )

    args$norm <- switch(
        args$norm,
        "auto" = DBA_NORM_NATIVE,
        "rle"  = DBA_NORM_RLE,
        "tmm"  = DBA_NORM_TMM,
        "lib"  = DBA_NORM_LIB
    )

    print("Converting the --score parameter to -log10 form")
    args$score <- -log10(args$score)             # need to convert it to -log10, as both pvalue and qvalue are -log10

    all_design_vars <- unique(all.vars(as.formula(args$design)))
    if(length(intersect(all_design_vars, ALL_CRITERIA)) != length(all_design_vars)){
        print("--design parameter includes not supported values")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    if(!is.null(args$groupby)){
        all_groupby_vars <- unique(args$groupby)
        if(length(intersect(all_groupby_vars, ALL_CRITERIA)) != length(all_groupby_vars)){
            print("--groupby parameter includes not supported values")
            quit(save = "no", status = 1, runLast = FALSE)
        }
        dba_groupby <- c()
        for (i in 1:length(all_groupby_vars)){
            dba_groupby <- append(
                dba_groupby,
                switch(
                    all_groupby_vars[i],
                    "Tissue"     = DBA_TISSUE,
                    "Factor"     = DBA_FACTOR,
                    "Condition"  = DBA_CONDITION,
                    "Treatment"  = DBA_TREATMENT,
                    "Replicate"  = DBA_REPLICATE
                )
            )
        }
        args$groupby_names <- args$groupby
        args$groupby_codes <- dba_groupby
    }

    return (args)
}


get_args <- function(){
    parser <- ArgumentParser(description="DiffBind Multi-factor Analysis")
    parser$add_argument(
        "--alignments",
        help="Sorted and indexed alignment files in bam format.",
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--peaks",
        help=paste(
            "Peak files in the MACS2 xls format. Number and order of the",
            "files should correspond to the files provided in --alignments",
            "parameter."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--aliases",
        help=paste(
            "Unique names for datasets provided in --alignments and --peaks",
            "parameters, no special characters or spaces are allowed. Number",
            "and order of the names should correspond to the values provided",
            "in --alignments and --peaks parameters."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "TSV/CSV metadata file to describe datasets provided in --alignments",
            "and --peaks parameters. First column should have the name 'sample',",
            "all other columns names should be selected from the following list:",
            "Tissue, Factor, Condition, Treatment, Replicate. The values from the",
            "'sample' column should correspond to the values provided in --aliases",
            "parameter. For a proper --contrast intepretation, values defined in",
            "each metadata column should not be used in any of the other columns.",
            "All metadata columns are treated as factors (no covariates are supported)."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--scoreby",
        help=paste(
            "Score metrics to exclude low quality peaks based on the",
            "threshold provided in the --score parameter.",
            "Default: pvalue"                                                              # https://support.bioconductor.org/p/66628/
        ),
        type="character",
        choices=c("pvalue", "qvalue"),
        default="pvalue"
    )
    parser$add_argument(
        "--score",
        help=paste(
            "Filtering threshold to keep only those peaks where the metric selected",
            "in --scoreby parameter is less than or equal to the provided value.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--minrpkm",
        help=paste(
            "Filtering threshold to keep only those peaks where the max RPKM for",         # in fact, this filter will be applied to DBA_SCORE_RPKM_MINUS,
            "all datasets is bigger than or equal to the provided value. Default: 1"       # but we don't have any controls so it's the same as DBA_SCORE_RPKM 
        ),
        type="double", default=1
    )
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
    parser$add_argument(
        "--minoverlap",
        help=paste(
            "Filtering threshold to keep only those peaks that are present in at",
            "least this many datasets when generating consensus set of peaks used",
            "in differential analysis. If this threshold has a value between zero",
            "and one, only those peaks will be included that are present in at least",
            "this proportion of datasets. When combined with --groupby parameter,",
            "--minoverlap threshold is applied per group, then union of the resulted",
            "peaks are used in the differential analysis. Default: 2"
        ),
        type="double", default=2
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column(s) from the metadata table to define datasets grouping. --minoverlap",
            "filtering threshold will be applied within each group independently. Union",
            "of the resulted peaks from each of the groups will be used in the differential",
            "analysis. Default: apply --minoverlap filtering threshold for all datasets",
            "jointly"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--design",
        help=paste(
            "Design formula comprised of the metadata columns names.",
            "It should start with ~."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--contrast",
        help=paste(
            "Contrast applied to the analysis results when calculating log2 fold changes.",
            "It should be formatted as a mathematical formula of values present in the",
            "metadata table. It is a required parameter if --method is set to edger. If not",
            "provided and --method is set to deseq2, the last term from the design formula",
            "will be used."
        ),
        type="character"
    )
    parser$add_argument(
        "--base",
        help=paste(
            "Base levels for each of the metadata columns. Number and order of the provided",
            "values should correspond to the metadata columns. Default: define base levels",
            "alphabetically."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--method",
        help=paste(
            "Method used in the differential binding analysis. Should be equal to either edger",
            "or deseq2. Default: deseq2"
        ),
        type="character", choices=c("edger", "deseq2"),
        default="deseq2"
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Normalization technique applied to the read counts before running differential",
            "binding analysis. When set to auto selects rle for deseq2 and tmm for edger.",
            "Default: auto"
        ),
        type="character", choices=c("auto", "rle", "tmm", "lib"),
        default="auto"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "Filtering threshold to report only differentially bound sites with adjusted",
            "P-value less than or equal to the provided value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on normalized read counts. Default:",
            "do not run clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--rowdist",
        help=paste(
            "Distance metric for HOPACH row clustering. Ignored if --cluster is not",
            "provided. Default: cosangle"
        ),
        type="character", default="cosangle",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--columndist",
        help=paste(
            "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
            "provided. Default: euclid"
        ),
        type="character", default="euclid",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for normalized read counts prior to running",
            "clustering by row. Ignored when --cluster is not row or both.",
            "Default: do not centered"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help=paste(
            "Output prefix for generated files"
        ),
        type="character", default="./diffbind"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    return (args)
}


# Parse arguments
args <- get_args()

print("Input parameters")
print(args)

print(paste("Loading metadata from", args$metadata))
metadata <- get_metadata(args)

# peaks will be filtered by --scoreby column to include only >= args$score values
# We already transformed args$score to -log10, that's why we use >=
# Scores column will be normalized from 0 to 1
dba_default_params <- list(
    minOverlap=1,                               # we should set it to 1, otherwise minOverlap is getting applied even when we load data
    peakCaller="macs",                          # we don't set peakFormat and scoreCol as they will be derived from the peakCaller value
    peakFormat="macs",
    scoreCol=args$scoreby,                      # either 7 (-log10 pvalue) or 9 (-log10 qvalue)
    filter=args$score,
    bLowerScoreBetter=FALSE,                    # both score columns are -log10, so lower values are not better
    config=list(                                # default config with certain parameters overwritten
        AnalysisMethod=args$method,
        th=args$padj,
        DataType=DBA_DATA_FRAME,                # to return data as dataframe
        RunParallel=TRUE, 
        minQCth=15,
        fragmentSize=125,                       # this value will be ignored when running dba.count with bUseSummarizeOverlaps=TRUE
        bCorPlot=FALSE,
        reportInit="DBA", 
        bUsePval=FALSE,
        design=TRUE,
        doBlacklist=FALSE,
        doGreylist=FALSE,
        cores=args$cpus
    )
)

print("Constructing DBA object. Forcing --minoverlap to 1")
dba_data <- do.call(dba, append(dba_default_params, list(sampleSheet=metadata$dba)))
print(
    paste0(
        "This DBA object includes merged with minOverlap=1 peaks ",
        "preliminary filtered by ", args$scoreby, " column to ",
        "contain only values >= ", args$score
    )
)
print(dba_data)
print(str(dba_data))

sample_groups <- c(paste(metadata$dba$SampleID, collapse="@"))
if(!is.null(args$groupby_names)){
    sample_groups <- append(
        sample_groups,
        metadata$dba %>% dplyr::group_by(across(all_of(args$groupby_names))) %>%
                         mutate(sample=paste(SampleID, collapse="@")) %>%
                         ungroup() %>%
                         select(sample) %>%
                         distinct() %>%
                         pull(sample)
    )
}
print("Sample groups for peakset overlap rate plots")
print(sample_groups)
for (i in 1:length(sample_groups)){
    current_samples <- unlist(strsplit(sample_groups[i], "@"))
    print(
        paste(
            "Processing", paste(current_samples, collapse=", "),
            "group of samples"
        )
    )
    current_mask = c()
    for (j in 1:length(current_samples)) {
        current_mask <- append(
            current_mask,
            which(as.vector(as.character(metadata$dba$SampleID)) == current_samples[j])
        )
    }
    export_overlap_plot(
        data=dba.overlap(dba_data, mask=current_mask, mode=DBA_OLAP_RATE),
        plot_title="Peakset overlap rate",
        plot_subtitle=paste(current_samples, collapse=", "),
        rootname=paste(args$output, "_pk_vrlp_s", i, sep=""),
        pdf=args$pdf
    )
}

# correlation heatmap is build based on the -log10(pvalue) or
# -log10(qvalue) columns from the loaded peak files. All peaks
# are added to the same list where the overlapping or adjacent
# peaks are merged. If several peaks from the specific peakset
# overlap the same merged peak, the best score is taken. If a
# certain peakset doesn't include a peak that overlaps with the
# merged peak, negative score is assigned.

export_corr_heatmap(                              # how it is calculated https://support.bioconductor.org/p/63034/
    data=dba_data,
    plot_title="Datasets correlation (all peaks)",
    rootname=paste(
        args$output, "all_pk_scr_corr", sep="_"
    ),
    pdf=args$pdf
)

print("Counting reads in peaks")
dba_count_default_params <- list(
    DBA=dba_data,
    bRemoveDuplicates=FALSE,                                  # set to FALSE because with bUseSummarizeOverlaps duplicates should be removed in advance
    bUseSummarizeOverlaps=TRUE,                               # this is a default value but we set it here explicitely, fragmentSize will be ignored
    summits=ifelse(args$summits >= 0, args$summits, FALSE),   # will make all consensus peaks have the same size 2*summits+1
    filter=args$minrpkm,                                      # keep only those peaks where max RPKM for all datasets is bigger or equal to this value
    bParallel=TRUE
)

if(!is.null(args$groupby_codes)){
    print(
        paste(
            "Searching for the consensus peaks with minimum peakset overlap",
            args$minoverlap, "applied to each dataset group separately"
        )
    )
    dba_peakset_default_params <- list(
        peak.caller="macs",
        peak.format="macs",
        scoreCol=args$scoreby,
        bLowerScoreBetter=FALSE,
        filter=args$score,
        DataType=DBA_DATA_FRAME
    )
    dba_data_temp <- do.call(
        dba.peakset,
        append(
            dba_peakset_default_params,
            list(
                DBA=dba_data,
                consensus=args$groupby_codes,
                minOverlap=args$minoverlap                               # applying --minoverlap per group
            )
        )
    )
    print("Temporary DBA object to obtain consensus peaks")
    print(dba_data_temp)
    print(str(dba_data_temp))

    if (!("Consensus" %in% names(dba_data_temp$masks))){
        print(
            paste(
                "--minoverlap parameter it too high for the provided --groupby",
                "or each group ended up having exactly one dataset. Try to use",
                "smaller --minoverlap or set it in a fraction form, alternatively,",
                "update --groupby parameter."
            )
        )
        quit(save = "no", status = 1, runLast = FALSE)
    }

    dba_data_temp <- do.call(
        dba,
        append(
            dba_default_params,
            list(
                DBA=dba_data_temp,
                mask=dba_data_temp$masks$Consensus
            )
        )
    )
    print("Temporary DBA object made of only consensus peaks")
    print(dba_data_temp)
    print(str(dba_data_temp))

    consensus_peaks <- do.call(
        dba.peakset,
        append(
            dba_peakset_default_params,
            list(
                DBA=dba_data_temp,
                minOverlap=1,                                             # force it to 1 to get the union of peaks from each of the groups
                bRetrieve=TRUE
            )
        )
    )
    print(paste("Found", nrow(consensus_peaks), "consensus peaks"))
    if(nrow(consensus_peaks) == 0){
        print("Found 0 consensus peaks. Exiting.")
        quit(save = "no", status = 1, runLast = FALSE)
    }
    print(head(consensus_peaks))

    dba_data <- do.call(
        dba.count,
        append(
            dba_count_default_params,
            list(
                minOverlap=1,                                             # force it to 1 to get the union of peaks from each of the groups
                peaks=consensus_peaks
            )
        )
    )
    export_corr_heatmap(                                                  # how it is calculated https://support.bioconductor.org/p/63034/
        data=dba_data,
        plot_title=paste0(
            "Datasets correlation (",
            ifelse(args$summits >= 0, paste0(args$summits, "bp rec. "), ""),
            "cons. peaks)"
        ),
        rootname=paste(
            args$output, "cns_pk_scr_corr", sep="_"
        ),
        pdf=args$pdf
    )
} else {
    print(
        paste(
            "Searching for the consensus peaks with minimum peakset overlap",
            args$minoverlap, "applied to all datasets jointly"
        )
    )
    dba_data <- do.call(
        dba.count,
        append(
            dba_count_default_params,
            list(
                minOverlap=args$minoverlap                                # to force using user-provided value
            )
        )
    )
    export_corr_heatmap(                                                  # how it is calculated https://support.bioconductor.org/p/63034/
        data=dba_data,
        plot_title=paste0(
            "Datasets correlation (",
            ifelse(args$summits >= 0, paste0(args$summits, "bp rec. "), ""),
            "cons. peaks)"
        ),
        rootname=paste(
            args$output, "cns_pk_scr_corr", sep="_"
        ),
        pdf=args$pdf
    )
}
print(
    paste0(
        "This DBA object includes peaks with ",
        "minOverlap=", args$minoverlap,
        " applied either to all or grouped datasets"
    )
)
print(dba_data)
print(str(dba_data))

export_corr_heatmap(
    data=dba_data,
    score=DBA_SCORE_READS,
    plot_title="Datasets correlation (raw reads)",
    rootname=paste(
        args$output, "rw_rds_corr", sep="_"
    ),
    pdf=args$pdf
)

print("Defining correct base levels to reorder DBA metadata")
correct_base_levels <- list()
for (i in 1:length(colnames(metadata$raw))){
    current_column <- colnames(metadata$raw)[i]
    current_levels <- levels(metadata$raw[[current_column]])
    correct_base_levels[[current_column]] <- current_levels[1]
}
print(correct_base_levels)

print("Reordering DBA metadata levels")
dba_data <- dba.contrast(                                # we need it only to make DBA object include correct metadata
    dba_data,                                            # otherwise it will use some default order which we can't reproduce
    design=args$design,
    reorderMeta=correct_base_levels
)
print("This is DBA object with reordered metadata levels")
print(dba_data)
print(str(dba_data))

print("Defining the contrast")
if(is.null(args$contrast)){
    print(
        paste(
            "Contrast is not provided. The last term",
            "from the design formula will be used."
        )
    )
    available_contrasts <- dba.contrast(
        dba_data,
        design=args$design,
        bGetCoefficients=TRUE
    )
    print("Available contrasts")
    print(available_contrasts)
    selected_contrast <- available_contrasts[length(available_contrasts)]
    print("Selected contrast")
    print(selected_contrast)
} else {
    print("Using user-provided contrast")
    selected_contrast <- get_contrast(
        design_formula=as.formula(args$design),
        metadata=metadata$raw,
        args=args
    )
}

dba_data <- dba.contrast(
    dba_data,
    design=args$design,                           # just in case set it here as well
    contrast=selected_contrast
)
print("This is DBA object after we added contrast")
print(dba_data)
print(str(dba_data))

print("Normalizing counts")
dba_data <- dba.normalize(                        # this doesn't know anything about design formula yet, so it's not take into account when normalizing
    dba_data,                                     # normalization somehow updates Score column from peaks?
    method=args$method,
    normalize=args$norm
)
print("This is DBA object after we normalized counts")
print(dba_data)
print(str(dba_data))

export_corr_heatmap(
    data=dba_data,
    score=DBA_SCORE_NORMALIZED,
    plot_title="Datasets correlation (norm. reads)",
    rootname=paste(
        args$output, "nr_rds_corr", sep="_"
    ),
    pdf=args$pdf
)

print("Performing the differential analysis")
dba_data <- dba.analyze(
    dba_data,
    method=args$method,
    bParallel=TRUE,
    bReduceObjects=FALSE                             # we want to have a complete DESeq2 or edgeR object
)
print("This is DBA object after we performed differential analysis")
print(dba_data)
print(str(dba_data))

norm_counts_data <- dba.report(
    dba_data,
    method=args$method,
    th=1,                                            # to include all differentially bound sites
    DataType=DBA_DATA_FRAME,
    bCalled=FALSE,
    bCounts=TRUE,
    bNormalized=TRUE
) %>% drop_na() %>%
      filter_all(all_vars(!is.infinite(.))) %>%
      mutate(feature=paste(Chr, paste(Start, End, sep="-"), sep=":")) %>%
      dplyr::rename("log2FoldChange"="Fold") %>%
      dplyr::rename("padj"="FDR") %>%
      dplyr::rename("baseMean"="Conc") %>%
      dplyr::rename("pvalue"="p-value") %>%
      remove_rownames() %>% column_to_rownames("feature")
print(head(norm_counts_data))

export_plots(norm_counts_data, metadata$raw, args)

print("Exporting differential binding sites")
export_data(
    norm_counts_data,                                                          # this is not filtered differentially expressed features
    location=paste(args$output, "diff_sts.tsv", sep="_"),
    digits=5
)

row_metadata <- norm_counts_data %>%
                select(log2FoldChange, padj) %>%
                filter(.$padj <= args$padj) %>%
                arrange(desc(log2FoldChange))
col_metadata <- metadata$raw %>%
                mutate_at(colnames(.), as.vector)                              # need to convert to vector, because in our metadata everything was a factor

print(
    paste(
        "Filtering normalized read counts matrix to include",
        "only differential binding sites with padj <= ", args$padj
    )
)
norm_counts_mat <- as.matrix(norm_counts_data[as.vector(rownames(row_metadata)), rownames(col_metadata)])
print("Size of the normalized read counts matrix after filtering")
print(dim(norm_counts_mat))
print(head(norm_counts_mat))

if (!is.null(args$cluster)){
    if (args$cluster == "column" || args$cluster == "both") {
        print("Clustering filtered read counts by columns")
        clustered_data = get_clustered_data(
            data=norm_counts_mat,
            center=NULL,                                              # centering doesn't influence on the samples order
            dist=args$columndist,
            transpose=TRUE
        )
        col_metadata <- cbind(col_metadata, clustered_data$clusters)  # adding cluster labels
        col_metadata <- col_metadata[clustered_data$order, ]          # reordering samples order based on the HOPACH clustering resutls
        print("Reordered samples")
        print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
        print("Clustering filtered normalized read counts by rows")
        clustered_data = get_clustered_data(
            data=norm_counts_mat,
            center=if(args$center) "mean" else NULL,                  # about centering normalized data https://www.biostars.org/p/387863/
            dist=args$rowdist,
            transpose=FALSE
        )
        norm_counts_mat <- clustered_data$counts                      # can be different because of centering by rows mean
        row_metadata <- cbind(row_metadata, clustered_data$clusters)  # adding cluster labels
        row_metadata <- row_metadata[clustered_data$order, ]          # reordering features order based on the HOPACH clustering results
        print("Reordered features")
        print(head(row_metadata))
    }
}

selected_sites <- makeGRangesFromDataFrame(
    row_metadata %>%
        mutate(rank=row_number()) %>%                                 # we use rank to maintain the  proper order of peaks
        mutate(ranges=rownames(.)) %>%                                # need to save rownames, if we remove them profileplyr and dba.plotProfile messes up peaks order
        separate(
            col="ranges",
            into=c("chr", "start", "end"),
            sep=paste0(":", "|", "-")
    ),
    keep.extra.columns=TRUE                                           # to keep all metadata columns
)

split_by <- NULL
if ("HCL" %in% colnames(row_metadata)){
    split_by <- "HCL"
} else if ("HCL.1" %in% colnames(row_metadata)){
    split_by <- "HCL.1"
}

if(!is.null(split_by)){
    print(paste("Splitting selected sites for profile heatmap by", split_by))
    selected_sites <- split(
        selected_sites,
        fct_inorder(as.vector(as.character(mcols(selected_sites)[[split_by]])))            # reorder levels by first appearance in a vector - gives the proper cluster order
    )
} else {
    print("Splitting selected sites for profile heatmap by log2FoldChange")
    selected_sites <- GRangesList(
        gain=selected_sites[selected_sites$log2FoldChange >= 0, ],
        loss=selected_sites[selected_sites$log2FoldChange < 0, ]
    )
}

print(head(selected_sites))

selected_samples <- list()                                                                 # correct order of sample based on col_metadata
for (i in 1:length(rownames(col_metadata))) {
    current_sample <- rownames(col_metadata)[i]
    selected_samples[[current_sample]] <- which(
        as.vector(as.character(metadata$dba$SampleID)) == current_sample
    )
}
print(head(selected_samples))

export_profile_heatmap(
    data=dba_data,
    sites=selected_sites,
    decreasing=FALSE,                                                                      # rank is a row number so we want to sort in increasing order
    samples=selected_samples,
    max_sites=Inf,                                                                         # to show all binding sites
    split_by=split_by,                                                                     # we need it only to reorder Binging sites if selected_sites split by cluster
    rootname=paste(args$output, "pk_prfl", sep="_"),
    width=ifelse(length(args$aliases) > 5, length(args$aliases) * 150, 800),
    pdf=args$pdf
)

# we do not reorder norm_counts_mat based on the clustering order
# because when exporting it to GCT we use row_metadata and col_metadata
# to force the proper order of rows and columns
print("Exporting normalized read counts to GCT format")
export_gct(
    counts_mat=norm_counts_mat,
    row_metadata=row_metadata,                                        # includes features as row names
    col_metadata=col_metadata,                                        # includes samples as row names
    location=paste(args$output, "nr_rds.gct", sep="_")
)