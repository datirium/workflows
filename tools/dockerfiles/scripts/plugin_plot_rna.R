#!/usr/bin/env Rscript
options(warn=-1)

suppressMessages(library(optparse))
suppressMessages(library(sqldf))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(BiocParallel))


get_coverage <- function(ranges, isoforms, bam_file, is_pair, is_dutp, threads) {
	plt_len <- 200
	dummy <- rep(0, plt_len + 1)
	is_ranges<-sqldf("SELECT ranges.chrom, ranges.txStart, ranges.txEnd, ranges.name, 0, ranges.strand, ranges.cdsStart,
							 ranges.cdsEnd, 0, ranges.exonCount, ranges.exonStarts, ranges.exonEnds
					  FROM   ranges, isoforms
					  WHERE  ranges.chrom=isoforms.Chrom AND
							 ranges.strand=isoforms.Strand AND
							 ranges.txStart=isoforms.TxStart AND
							 ranges.txEnd=isoforms.TxEnd AND
							 ranges.name=isoforms.RefseqId AND
							 ranges.txEnd-ranges.txStart > 1000 AND
							 isoforms.Rpkm > 10")
	is_count <- length(is_ranges$exonCount)
	if(is_count < 5) return(dummy)
	namesl <- list()
	irl <- list()
	for(i in seq_len(is_count)){
		st <- as.numeric(strsplit(is_ranges$exonStarts[i], ',')[[1]])
		en <- as.numeric(strsplit(is_ranges$exonEnds[i], ',')[[1]])
		names <- paste0(rep(is_ranges$chrom[i], is_ranges$exonCount[i]), ":", st, "-", en)
		ir <- IRanges(start=st, end=en, names=names)
		metadata(ir) <- list(strand=is_ranges$strand[i])
		namesl <- append(namesl, is_ranges$chrom[i])
		irl <- append(irl, list(ir))
	}
	what <- c("pos","strand")
	which <- RangesList(irl)
	names(which) <- namesl
	flags <- scanBamFlag(isProperPair = is_pair,
                         isUnmappedQuery = FALSE,
                         hasUnmappedMate = FALSE,
                         isFirstMateRead = is_pair)
	param <- ScanBamParam(which = which, what = what, flag = flags)
	bam <- scanBam(bam_file, index = bam_file, param = param)
	cov <- bplapply(irl,
                    function(ir){
                        genewidth <- sum(width(ir))
                        bin <- genewidth / plt_len
                        cov <- c(rep(0, plt_len + 1))
                        for(i in seq_len(length(ir))){
                            abs_pos <- 0
                            idx <- 1
                            cwidth <- 0
                            # Strand specificity
                            if(is_dutp) {
                                strand <- "+"
                                if(metadata(ir)$strand == "+")
                                    strand <- "-"
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos[bam[[names(ir)[i]]]$strand == strand])
                            } else {
                                tmp1 <- rle(bam[[names(ir)[i]]]$pos)
                            }
                            for(j in seq_len(length(tmp1$lengths))){
                                pos <- tmp1$values[j]
                                if(pos > end(ir)[idx]){
                                    cwidth <- sum(width(ir[end(ir) < pos]))
                                    idx <- sum(start(ir) < pos)
                                if(idx > length(ir))
                                    break
                            }
                                if(pos - start(ir)[idx] + 1 > width(ir)[idx])
                                    next
                                abs_pos <- (pos - start(ir)[idx]) + cwidth
                                il <- floor(abs_pos / bin)
                                cov[il+1] = cov[il+1] + tmp1$lengths[j]
                            }
                        }
                        cov <- cov / bin
                        if(metadata(ir)$strand == "+") { cov } else { rev(cov) }
                    },
                    BPPARAM = MulticoreParam(workers = threads))
    return(Reduce("+", cov) / is_count)
}

option_list <- list(make_option(c("-a", "--annotation"), type="character", help="Path to annotation file"),
                    make_option(c("-b", "--bam"),        type="character", help="Path to BAM file (+BAI)"),
                    make_option(c("-i", "--isoforms"),   type="character", help="Path to isoforms file"),
                    make_option(c("-s", "--stat"),       type="character", help="Path to statistics file"),
                    make_option(c("-o", "--output"),     type="character", help="Output file prefix", default="./"),
                    make_option(c("-p", "--pair"),       type="logical",   help="Is paired end", action="store_true", default = FALSE),
                    make_option(c("-d", "--dutp"),       type="logical",   help="Is dUTP", action="store_true", default = FALSE),
                    make_option(c("-t", "--threads"),    type="integer",   help="Threads number", default = 1));

opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);

png(filename=paste(args$output, "%03d.png", sep=""))

colvar <- 5
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f"))(colvar)

ranges <- read.table(args$annotation, sep="\t", header=TRUE, stringsAsFactors=FALSE)
isoforms <- read.table(args$isoforms, sep=",", header=TRUE, stringsAsFactors=FALSE)
stat <- read.table(args$stat, sep=" ", header=FALSE, stringsAsFactors=FALSE)

tags_mapped = as.numeric(stat[2])
cov_norm <- get_coverage(ranges, isoforms, args$bam, args$pair, args$dutp, args$threads)/(tags_mapped/1000000)
write.table(cov_norm,
            file = paste(args$output, "cov.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
plot(cov_norm,
     type="l",
     xaxt = "n",
     main="Gene body average tag density",
     ylab="Average Tag Density (per percentile)",
     xlab="Gene body percentile (5'->3')",
     lwd=3, col=icolor[floor(runif(1)*colvar)+1])
axis(1, at=seq(0,200,40), labels=seq(0,100,20), las=1)

rpkm_hist <- hist(isoforms$Rpkm[isoforms$Rpkm>2 & isoforms$Rpkm<500],
                  main="RPKM distribution",
                  breaks=1000,
                  xlab="rpkm>2 & rpkm<500",
                  col=icolor[floor(runif(1)*colvar)+1])

write.table(data.frame(mids=rpkm_hist$mids, counts=rpkm_hist$counts, density=rpkm_hist$density),
            file = paste(args$output, "rpkm.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)


graphics.off()