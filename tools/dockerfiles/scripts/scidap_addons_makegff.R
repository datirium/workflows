#!/usr/bin/env Rscript
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)

BASE_FILENAME = as.character(args[1])
GFF_FILENAME = as.character(args[2])
CONTROL_FILENAME = as.character(args[3])

base_islands <- read.table(BASE_FILENAME, sep="\t", header=TRUE)

CHROM_COL_NAME = "chrom"
if ( !(CHROM_COL_NAME %in% colnames(base_islands)) )
{
    CHROM_COL_NAME = "chr"
}

base_islands = base_islands[base_islands[ ,CHROM_COL_NAME] != "chrM", ]

if(!is.na(CONTROL_FILENAME)){
    control_islands <- read.table(CONTROL_FILENAME, sep="\t", header=TRUE)
    control_islands = control_islands[control_islands[ ,CHROM_COL_NAME] != "chrM", ]

    extended = 100000
    control_islands.granges = GRanges(seqnames=control_islands[, CHROM_COL_NAME], ranges = IRanges(start = control_islands[, "start"] - extended, end = control_islands[, "end"] + extended))
    control_islands.granges = reduce(control_islands.granges)
    base_islands.granges = GRanges(seqnames=base_islands[, CHROM_COL_NAME], ranges = IRanges(start = base_islands[, "start"], end= base_islands[, "end"]))
    mark = as.matrix(findOverlaps(base_islands.granges, control_islands.granges, ignore.strand = T))
    base_islands = base_islands[-mark[ ,1], ]
}

gff_islands = cbind.data.frame(base_islands[,CHROM_COL_NAME],rownames(base_islands),"", base_islands[ ,"start"],base_islands[ ,"end"],"",".","",rownames(base_islands))
write.table(gff_islands, GFF_FILENAME, quote=FALSE, col.names=FALSE, row.names=FALSE, append=FALSE, sep="\t")