#!/usr/bin/env Rscript
options(warn=-1)

suppressMessages(library(optparse))
suppressMessages(library(Rsamtools))


option_list <- list(make_option(c("-i", "--islands"), type="character", help="Path to the islands file (Iaintersect or MACS2 xls output)"),
                    make_option(c("-b", "--bam"),     type="character", help="Path to BAM file (+BAI)"),
                    make_option(c("-o", "--output"),  type="character", help="Output file prefix", default="./"));


opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
png(filename=paste(args$output, "%03d.png", sep=""))


colvar<-5
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f"))(colvar)

islands <- read.table(args$islands, sep="\t", header=TRUE, stringsAsFactors=FALSE)

sorted_pileup <- sort(islands$pileup, decreasing=T)
write.table(sorted_pileup,
            file = paste(args$output, "pileup.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
plot(sorted_pileup,
     main=paste("Rank plot"),
     xlab="Ranked peak number",
     ylab="Pileup",
     type="l",
     lwd=3,
     col=icolor[floor(runif(1)*colvar)+1])

length_hist <- hist(islands$length,
                   breaks=100,
                   main=paste("Islands distribution"),
                   xlab="Island length",
                   col=icolor[floor(runif(1)*colvar)+1])

write.table(data.frame(mids=length_hist$mids, counts=length_hist$counts, density=length_hist$density),
            file = paste(args$output, "length.tsv", sep=""),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)

what <- c("rname", "strand", "pos", "qwidth","seq")
param <- ScanBamParam(what=what)
bam <- scanBam(args$bam, param=param)
bam.dataframe = lapply(bam, function(x)do.call("DataFrame", x))
bam.dataframe2 = bam.dataframe[[1]][!is.na(bam.dataframe[[1]][["rname"]]), ]
chrs = unique(bam.dataframe2[["rname"]])
reads.mat = t(data.frame(read = rep(0, length(chrs))))
colnames(reads.mat) = chrs
for(chr in chrs){
    reads.mat[1, chr] = 100 * (nrow(bam.dataframe2[bam.dataframe2[["rname"]] == chr, ]) / nrow(bam.dataframe2))
}

write.table(t(reads.mat),
            file = paste(args$output, "reads.tsv", sep=""),
            sep="\t",
            row.names=TRUE,
            col.names=FALSE,
            quote=FALSE)
bp <- barplot(reads.mat,
              names.arg=colnames(reads.mat),
              col="yellowgreen",
              ylab="%",
              main="Chromosomal Distribution of Reads",
              las=2,
              mar=c(8,4,6,2))
text(c(bp),
     c(reads.mat),
     labels=round(c(reads.mat),digits=2),
     pos=3,
     offset=1,
     xpd=NA,
     srt=90)