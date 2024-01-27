#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(scipen=999)

##########################################################################################
#
# v0.0.1
# - normalizes RNA-Seq data using spike-in of ERCC ExFold mix 1, dilution factor, and uL per million cells
#
##########################################################################################

# character vector of positional arguments
args <- commandArgs(trailingOnly = TRUE)
dilution_factor <- as.numeric(args[1])                  # dilution factor used for sample
uL_per_M_cells <- as.numeric(args[2])                   # volume of spike-in per 1M cells in sample
d <- read.table(args[3], sep='\t', header=T)            # counts of ercc aligned reads
expected_d <- read.table(args[4], sep='\t', header=T)   # expected ercc counts
counts <- read.table(args[5], sep=',', header=T)        # target organism transcript count data

### EXAMPLE INPUTS  ###
#dilution_factor <- 0.10
#uL_per_M_cells <- 1
#d <- read.table('ercc_counts.tsv', sep='\t', header=T)
#expected_d <- read.table('/dockerdata/ercc_exfold_mix1_expected_counts.tsv', sep='\t', header=T)
#counts <- read.table('/data/scidap_workflow_testing/workflow_rnaseq_ercc_spikein/testing_docker_image/read_1.isoforms.csv', sep=',', header=T)

#   transform 'molecules per uL' to 'molecules per cell' for mix1 expected molecules
#       multiply by dilution factor
#       multiply by uL/10^6 cells
#       rename column
my_function <- function(x) (x*dilution_factor*uL_per_M_cells)/1000000
expected_d[c('molecules_per_uL_mix1')] <- lapply(expected_d[c('molecules_per_uL_mix1')], my_function)
colnames(expected_d)[2] <- 'molecules_per_cell'
#       merge expected_d (x) and d (y) counts on ercc id, retain rows where molecules_per_cell > 1
ab <- merge(expected_d, d, by="ERCC_ID")
ab <- ab[ab$count >= 1 ,]
#ab[is.na(ab)] <- 0  # fill missing id "count" with NA, then replace with 0 (NOT DOING THIS)
#       make linear regression model where lm(Y ~ X, data=data.frame)
model <- lm(log10(count) ~ log10(molecules_per_cell), data=ab)
print('SUMMARY OF LOG10 REGRESSION MODEL:')
summary(model)
#   just checking r^2 for raw data
model_raw <- lm(count ~ molecules_per_cell, data=ab)
print('SUMMARY OF RAW REGRESSION MODEL:')
summary(model_raw)

#   extract model coefficients and r^2 (FOR PLOT)
m <- signif(model$coefficients[2] ,digits=3)
b <- signif(model$coefficients[1] ,digits=3)
summary <- summary(model)
r2 <- signif(summary$r.squared ,digits=3)


#   scatter plot
#       https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/plot
pdf(file = "ercc_expected_v_actual_count_plot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6.4)
plot(log10(ab$molecules_per_cell), log10(ab$count), main="ERCC molecules per cell counts (log10)\nExpected vs Observed", xlab="log10(expected molecules per cell)", ylab="log10(observed counts)", pch=19)
#       Add fit lines
abline(lm(log10(count) ~ log10(molecules_per_cell), data=ab), col="red")
legend('bottomright', bty="n", legend=c(paste('linear regression\nr^2=',r2,'; ','y=',m,'x+',b,sep="")),
        col=c('red'), pch=19)
dev.off()



#   extract model coefficients and r^2 (FOR NORMALIZATION)
m <- model$coefficients[2]
b <- model$coefficients[1]


#   normalize transcript RPKM counts to spike-in
rpkm_norm <- counts
#       transform rpkm to log10
my_function <- function(x) log10(x)
counts[c('Rpkm')] <- lapply(rpkm_norm[c('Rpkm')], my_function)
#       apply each count to linear function to predict normalized count
my_function <- function(x) (m*x)+b
counts[c('Rpkm')] <- lapply(counts[c('Rpkm')], my_function)
#       convert out of log10
my_function <- function(x) as.integer(10^x)
counts[c('Rpkm')] <- lapply(counts[c('Rpkm')], my_function)
#       rename column header
#colnames(rpkm_norm)[8] <- 'ercc_norm_rpkm'     # ACTUALL DON'T, this will cause group_isoforms step to fail

#   normalize transcript Total counts to spike-in
total_norm <- counts
#       transform rpkm to log10
my_function <- function(x) log10(x)
counts[c('TotalReads')] <- lapply(total_norm[c('TotalReads')], my_function)
#       apply each count to linear function to predict normalized count
my_function <- function(x) (m*x)+b
counts[c('TotalReads')] <- lapply(counts[c('TotalReads')], my_function)
#       convert out of log10
my_function <- function(x) as.integer(10^x)
counts[c('TotalReads')] <- lapply(counts[c('TotalReads')], my_function)

#       save to new csv file
write.csv(counts, "isoforms.ercc_norm_rpkm.csv-hasquotes", row.names=FALSE)
