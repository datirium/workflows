#!/bin/bash

printf "$(date)\nLog file for run_ercc_norm.sh\n\n"

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_ercc_norm.sh\`:
    Wrapper for building linear regression function from ERCC ExFold mix 1 RPKM (molecule per cell vs RPKM), and applying this for normalization of RNA-Seq RPKM count data.

    Primary Output files:
     - rpkm_erccnorm_counts.tsv


PARAMS:
 -h  help	show this message
 -t  INT	number of threads
 -u  FILE   array of unaligned "R1,R2" reads post-primary alignment
 -d  FLOAT  dilution factor used for ERCC ExFold mix 1 before spike-in
 -m  FLOAT  volume of ERCC ExFold mix 1 spike-in to sample per million cells
 -c  FILE   csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)


BismarkCov formatted bed:
    https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
    The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

____________________________________________________________________________________________________
References:

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:u:d:m:c:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		u) R1R2_UNALIGNED=$OPTARG ;;
        d) DILTUOIN_FACTOR=$OPTARG ;;
		m) UL_PER_M_CELLS=$OPTARG ;;
		c) RNASEQ_COUNTS=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$R1R2_UNALIGNED" == "" ]]; then echo "error: required param missing (-a), array of unaligned R1,R2 reads post-primary alignment"; exit; fi
if [[ "$DILTUOIN_FACTOR" == "" ]]; then echo "warning: required param missing (-d), dilution factor used for ERCC ExFold mix 1 before spike-in"; exit; fi
if [[ "$UL_PER_M_CELLS" == "" ]]; then echo "warning: required param missing (-m), volume of ERCC ExFold mix 1 spike-in to sample per million cells"; exit; fi
if [[ "$RNASEQ_COUNTS" == "" ]]; then echo "warning: required param missing (-c), csv file containing isoform counts (format: RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm)"; exit; fi
# print inputs to stdout
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-u, \$R1R2_UNALIGNED, $R1R2_UNALIGNED\n"
printf "\t-d, \$DILTUOIN_FACTOR, $DILTUOIN_FACTOR\n"
printf "\t-m, \$UL_PER_M_CELLS, $UL_PER_M_CELLS\n"
printf "\t-c, \$RNASEQ_COUNTS, $RNASEQ_COUNTS\n"
printf "\n\n"



#	MAIN
#===============================================================================
R1_UNALIGNED=$(echo $R1R2_UNALIGNED | sed 's/,/\n/' | head -1)
R2_UNALIGNED=$(echo $R1R2_UNALIGNED | sed 's/,/\n/' | tail -1)
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