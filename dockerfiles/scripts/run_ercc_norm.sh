#!/bin/bash
printf "$(date)\nLog file for run_ercc_norm.sh\n\n"
# export paths (just in case docker failed to do so)
export PATH="$PATH:/dockerdata/apps/src/bowtie2-2.5.2"
export PATH="$PATH:/dockerdata/apps/src/samtools-1.12"

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_ercc_norm.sh\`:
    Wrapper for building linear regression function from ERCC ExFold mix 1 RPKM (molecule per cell vs RPKM), and applying this for normalization of RNA-Seq RPKM count data.

    Primary Output files:
     - unaligned_pairs-to-ERCC.sam
     - ercc_counts.tsv
     - isoforms.ercc_norm_rpkm.csv
     - ercc_expected_v_actual_count_plot.pdf


PARAMS:
 -h  help	show this message
 -t  INT	number of threads
 -u  FILE   array of unaligned "R1,R2" reads post-primary alignment
 -d  FLOAT  dilution factor used for ERCC ExFold mix 1 before spike-in
 -m  FLOAT  volume of ERCC ExFold mix 1 spike-in to sample per million cells
 -c  FILE   csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)

____________________________________________________________________________________________________
References:
    Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.
    Twelve years of SAMtools and BCFtools. Danecek et al. GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:u:v:d:m:c:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		u) R1_UNALIGNED=$OPTARG ;;
		v) R2_UNALIGNED=$OPTARG ;;
        d) DILTUOIN_FACTOR=$OPTARG ;;
		m) UL_PER_M_CELLS=$OPTARG ;;
		c) RNASEQ_COUNTS=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$R1_UNALIGNED" == "" ]]; then echo "error: required param missing (-u), array of unaligned R1 reads post-primary alignment"; exit; fi
if [[ "$R2_UNALIGNED" == "" ]]; then echo "error: required param missing (-v), array of unaligned R2 reads post-primary alignment"; exit; fi
if [[ "$DILTUOIN_FACTOR" == "" ]]; then echo "warning: required param missing (-d), dilution factor used for ERCC ExFold mix 1 before spike-in"; exit; fi
if [[ "$UL_PER_M_CELLS" == "" ]]; then echo "warning: required param missing (-m), volume of ERCC ExFold mix 1 spike-in to sample per million cells"; exit; fi
if [[ "$RNASEQ_COUNTS" == "" ]]; then echo "warning: required param missing (-c), csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)"; exit; fi
# print inputs to stdout
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-u, \$R1_UNALIGNED, $R1_UNALIGNED\n"
printf "\t-u, \$R2_UNALIGNED, $R2_UNALIGNED\n"
printf "\t-d, \$DILTUOIN_FACTOR, $DILTUOIN_FACTOR\n"
printf "\t-m, \$UL_PER_M_CELLS, $UL_PER_M_CELLS\n"
printf "\t-c, \$RNASEQ_COUNTS, $RNASEQ_COUNTS\n"
printf "\n\n"



#	MAIN
#===============================================================================
# map unaligned reads to ERCC sequences
bowtie2 -p $THREADS -x /dockerdata/refs/ERCC92.fa -1 $R1_UNALIGNED -2 $R2_UNALIGNED -S unaligned_pairs-to-ERCC.sam
# get counts of each ERCC sequence
printf "ERCC_ID\tcount\n" > ercc_counts.tsv
samtools view -F 260 unaligned_pairs-to-ERCC.sam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' >> ercc_counts.tsv
# normalize RPKM values in the input rna-seq count file
run_ercc_regression.R $DILTUOIN_FACTOR $UL_PER_M_CELLS ercc_counts.tsv /dockerdata/ercc_exfold_mix1_expected_counts.tsv $RNASEQ_COUNTS


#	OUTPUTS
#===============================================================================

# Rscript outputs:
#   ercc_expected_v_actual_count_plot.pdf (grab as cwl workflow output)
#   isoforms.ercc_norm_rpkm.csv-hasquotes

#   format for diffexp integration (probably edgeR)
sed 's/"//g' isoforms.ercc_norm_rpkm.csv-hasquotes > isoforms.ercc_norm_rpkm.csv



# clean up
rm -rf "$workdir"