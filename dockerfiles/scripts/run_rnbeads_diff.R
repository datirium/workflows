#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

#suppressMessages(library(argparse))
suppressMessages(library(RnBeads))


##########################################################################################
#
# v0.0.1
# - use RnBeads differential methylation analysis
# - site specific and region analysis enabled
# - computes enrichments for Gene Ontology (GO) terms
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
anno_file <- args[1]
outdir <- args[2]
datadir <- args[3]
genome <- args[4]
threads <- args[5]

# Directory where your data is located
#sample.annotation <- read.table(anno_file, header=T, sep=",")  # manual specifies this should be a df, but it is not accepted as input as part of data.source 
sample.annotation <- anno_file
class(sample.annotation)
sample.annotation

# Directory where the output should be written to
analysis.dir <- outdir
# Directory where the report files should be written to
report.dir <- file.path(analysis.dir, "reports")
rnb.initialize.reports(report.dir)
# Set some analysis options
rnb.options(filtering.sex.chromosomes.removal = TRUE, identifiers.column="Sample_ID", differential.comparison.columns="condition", import.bed.style="bismarkCov", assembly=genome)

# characterize differentially methylated regions by computing enrichments
# for Gene Ontology (GO) terms and for custom genomic annotations using the LOLA tool.
rnb.options("differential.enrichment.go"=TRUE)
#rnb.options("differential.enrichment.lola"=TRUE)   # LOLA download fails, leaving out for now

## Restrict logging to the console only
logger.start(fname = NA)
# setup parallel processing
parallel.isEnabled()
num.cores <- as.numeric(threads)
parallel.setup(num.cores)
parallel.isEnabled()
if (parallel.isEnabled()) parallel.getNumWorkers()

## Data import
#        ‘"bs.bed.dir"’            ‘list’ or ‘character’                ‘1..3’
#   (1) Directory with BED files each giving a DNA methylation profile of a sample;
#   (2) a sample annotation table as a ‘data.frame’ or the name of the corresponding file; 
#       In case only the sample sheet is provided as the second element of the data.source list
#       (the first element can be set to NULL), the provided sample sheet should contain absolute
#       paths to the bed files.
#data.dir <- NULL   # apparently this doesn't work! Contacted devs and they confirmed this option is not currently supported by RnBeads (working on a fix)
#Error in logger.error(txt) : 
#  invalid data.source parameter, bed.dir is not found, or is not directory
#Calls: rnb.run.import ... rnb.step.import -> rnb.execute.import -> rnb.error -> logger.error
#Execution halted

data.dir <- file.path(datadir)
data.dir
data.source <- c(data.dir, sample.annotation)
result <- rnb.run.import(data.source=data.source, data.type="bs.bed.dir", dir.reports=report.dir)
rnb.set <- result$rnb.set

## Quality Control
rnb.run.qc(rnb.set, report.dir)

## Preprocessing
rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports=report.dir)$rnb.set

## Data export
rnb.options(export.to.csv = TRUE)
rnb.run.tnt(rnb.set, report.dir)

## Differential methylation
#   min.group.size=2, if either group has <2 samples per condition, diffmeth will not be run
rnb.run.differential(rnb.set, report.dir)

## disable parallel processing
parallel.teardown()