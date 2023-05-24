#!/bin/bash

#   Shell wrapper for genelists from atac/chip/crt and rna data
#
##########################################################################################
#
# v1.0.0
# - 2023
#
##########################################################################################
printf "$(date)\nLog file for run_genelists.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_genelists.sh\`:
Shell wrapper for producing data file for gene list(s) heatmap. Uses both ATAC/ChIP/CRT and RNA-Seq data.

Primary Output files:
 - 

PARAMS:
    SECTION 1: general
    -h	help		show this message
    -t  INT			number of threads
	-a	ARRAY		array of genelist sample names (no commas in names)
	-b  FILE ARRAY	array of associated annotation files for each gene list from (-c), with header
    -c  FILE ARRAY	array of filtered gene list TSVs (must be headerless, columns are: chr, txStart, txEnd, geneID, L2FC, Strand)
	-d	ARRAY		array of sample names from NA binding experiments (no commas in names)
	-e	ARARY		array of sample names from RNA-Seq experiments (no commas in names)
	-f	FILE ARRAY	array of BAM files from NA binding experiments
	-g	FILE ARARY	array of expression table files from RNA-Seq experiments	



____________________________________________________________________________________________________
References:
 - 
    
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:a:b:c:d:e:f:g:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		a) GENELIST_NAMES=$OPTARG ;;
		b) GENELIST_ANNOTATION_FILES=$OPTARG ;;
		c) GENELIST_FILTERED_FILES=$OPTARG ;;
		d) NAMES_NABIND=$OPTARG ;;
		e) NAMES_RNASEQ=$OPTARG ;;
		f) FILES_NABIND_BAM=$OPTARG ;;
		g) FILES_RNASEQ_EXP=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$GENELIST_NAMES" == "" ]]; then echo "error: required param missing (-a), array of genelist sample names"; exit; fi
if [[ "$GENELIST_ANNOTATION_FILES" == "" ]]; then echo "error: required param missing (-b), array of genelist annotation (feature_file.tsv) files"; exit; fi
if [[ "$GENELIST_FILTERED_FILES" == "" ]]; then echo "error: required param missing (-c), array of genelist filtered (feature_file) files"; exit; fi
if [[ "$NAMES_NABIND" == "" ]]; then echo "error: required param missing (-d), array of NA binding sample names"; exit; fi
if [[ "$NAMES_RNASEQ" == "" ]]; then echo "error: required param missing (-e), array of RNA-Seq sample names"; exit; fi
if [[ "$FILES_NABIND_BAM" == "" ]]; then echo "error: required param missing (-f), array of NA binding BAM files"; exit; fi
if [[ "$FILES_RNASEQ_EXP" == "" ]]; then echo "error: required param missing (-g), array of RNA-Seq expression files"; exit; fi


# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-a, \$GENELIST_NAMES, $GENELIST_NAMES\n"
printf "\t-b, \$GENELIST_ANNOTATION_FILES, $GENELIST_ANNOTATION_FILES\n"
printf "\t-c, \$GENELIST_FILTERED_FILES, $GENELIST_FILTERED_FILES\n"
printf "\t-d, \$NAMES_NABIND, $NAMES_NABIND\n"
printf "\t-e, \$NAMES_RNASEQ, $NAMES_RNASEQ\n"
printf "\t-f, \$FILES_NABIND_BAM, $FILES_NABIND_BAM\n"
printf "\t-g, \$FILES_RNASEQ_EXP, $FILES_RNASEQ_EXP\n"
printf "\n\n"


#	MAIN
#===============================================================================
printf "generating master samplesheet from all inputs\n"
#	start counter for genelist arrays
list_counter=0
#	turn genelist annotion files and genelist names into arrays
annotations_array=($(echo "$GENELIST_ANNOTATION_FILES" | sed 's/,/\n/g'))
names_array=($(echo "$GENELIST_NAMES" | sed 's/,/\n/g'))
for f in $(echo "$GENELIST_FILTERED_FILES" | sed 's/,/\n/g'); do
	# loop through NA binding names
	#	start counter for name array
	name_counter=0
	#	turn bam files into an array
	bam_array=($(echo "$FILES_NABIND_BAM" | sed 's/,/\n/g'))
	for n in $(echo "$NAMES_NABIND" | sed 's/,/\n/g'); do
		# print formatted samplesheet row
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "na-binding" $n ${bam_array[name_counter]}
		((name_counter++))
	done

	# loop through RNA names
	#	start counter for name array
	name_counter=0
	#	turn bam files into an array
	exp_array=($(echo "$FILES_RNASEQ_EXP" | sed 's/,/\n/g'))
	for n in $(echo "$NAMES_RNASEQ" | sed 's/,/\n/g'); do
		# print formatted samplesheet row
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "rna-seq" $n ${exp_array[name_counter]}
		((name_counter++))
	done

	((list_counter++))
done > master_samplesheet.tsv


# initialize output files with headers
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tavg_depth\n" > output_na-binding.tmp
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tTotalReads\tRpkm\n" > output_rna-seq.tmp
# function for building tag density window data and heatmap expression data for each experiment type
get_data()
{
	printf "\tget_data(), processing $1\n"
	genelist_number=$(printf "$1" | cut -f1)
	genelist_name=$(printf "$1" | cut -f2)
	genelist_file=$(printf "$1" | cut -f3)
	genelist_annotation_file=$(printf "$1" | cut -f4)
	experiment_type=$(printf "$1" | cut -f5)	# determines data extraction method of file at $6 ("no-binding" or "rna-seq")
	sample_name=$(printf "$1" | cut -f6)
	sample_data=$(printf "$1" | cut -f7)			# either bam or expression table

	if [[ $experiment_type == "na-binding" ]]; then
		# user input variables
		bam=$sample_data
		bn=$(basename $bam | sed 's/\..*//')
		window="5000"   # this is +/- $window of the TSS (for +strand TSS=$from, for -strand TSS=$to)
		# for each gene in filtered genelist file, get tag denisty (average read depth) per window step
		cat $genelist_file | while read filtered_gene; do
			timestamp=$(date +%s)
			printf "\t\tna-binding, processing gene:\t$filtered_gene\n"
			chr=$(printf "$filtered_gene" | cut -f1)  # genelist col1
			txStart=$(printf "$filtered_gene" | cut -f2)  # genelist col2
			txEnd=$(printf "$filtered_gene" | cut -f3)  # genelist col3
			geneid=$(printf "$filtered_gene" | cut -f4)  # genelist col4
			strand=$(printf "$filtered_gene" | cut -f6)  # genelist col6
			#   generate 3 column tsv of chr, position, depth that is +/-5Kbp of TSS
			if [[ "$strand" == "+" ]]; then
				from=$(echo "$txStart" | awk -v w=$window '{print($0-w)}')
				to=$(echo "$txStart" | awk -v w=$window '{print($0+w)}')
				samtools depth -a -r $chr:$from-$to $bam > $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth
			else
				from=$(echo "$txEnd" | awk -v w=$window '{print($0-w)}')
				to=$(echo "$txEnd" | awk -v w=$window '{print($0+w)}')
				samtools depth -a -r $chr:$from-$to $bam > $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth
			fi
			#   reduce to 50 windows
			#    - 100 bp per window frame in steps of 50 bp to cover the 10,000 bp +/-5Kbp of TSS
			#    - i.e. there will be a 50 bp overlap between 100 bp windows to smoothed depth average values, and total of 99 windows
			seq 100 50 5000 | while read step; do
				# calculate average depth for window
				avg_depth=$(head -$step $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth | tail -100 | awk -F'\t' '{x+=$3}END{printf("%.2f",x/NR)}')
				# print formatted output for na-binding experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\t$step\t$avg_depth\n"
			done >> output_na-binding.tmp
			rm $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth
		done
	elif [[ $experiment_type == "rna-seq" ]]; then
		# user input variables
		exp=$sample_data
		bn=$(basename $exp | sed 's/\..*//')
		# for each gene in filtered genelist file, get TotalReads and Rpkm values
		cat $genelist_file | while read filtered_gene; do
			printf "\t\trna-seq, processing gene:\t $filtered_gene\n"
			chr=$(printf "$filtered_gene" | cut -f1)  # genelist col1
			txStart=$(printf "$filtered_gene" | cut -f2)  # genelist col2
			txEnd=$(printf "$filtered_gene" | cut -f3)  # genelist col3
			geneid=$(printf "$filtered_gene" | cut -f4)  # genelist col4
			strand=$(printf "$filtered_gene" | cut -f6)  # genelist col6
			# search for each gene
			g=$(grep -P -m 1 "\t$geneid\t$chr\t$txStart\t$txEnd\t$strand" $exp)
			if [[ $g == "" ]]; then
				printf "\t\t\tWARNING, gene does not exist in expression data file: $geneid,$chr,$txStart,$txEnd,$strand"
			else
				# still include the steps (tss_window) to pad for heatmap
				seq 100 50 5000 | while read step; do
					TotalReads=$(printf "$g" | cut -f7)
					Rpkm=$(printf "$g" | cut -f8)
					# print formatted output for rna-seq experiments
					printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\t$step\t$TotalReads\t$Rpkm\n"
				done >> output_rna-seq.tmp
			fi
		done
	fi
}
export -f get_data


printf "\n\n"
printf "running master_samplesheet.tsv through parallel\n"
parallel --arg-file master_samplesheet.tsv --jobs="$THREADS" get_data
# ensure rows are unique
head -1 output_na-binding.tmp > output_na-binding.tsv
head -1 output_rna-seq.tmp > output_rna-seq.tsv
tail -n+2 output_na-binding.tmp | sort | uniq >> output_na-binding.tsv
tail -n+2 output_rna-seq.tmp | sort | uniq >> output_rna-seq.tsv
rm output_na-binding.tmp
rm output_rna-seq.tmp



# merge na-binding and rna-seq outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window
#		row9 - data_type (avg_depth, TotalReads, or Rpkm)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (avg_depth, TotalReads, Rpkm)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
#	merge rna-seq and na-binding data (13 total columns)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"TotalReads",$11)}}' output_rna-seq.tsv >> output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"Rpkm",$12)}}' output_rna-seq.tsv >> output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"avg_depth",$11)}}' output_na-binding.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_rna-seq.tsv > output_row_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_rna-seq.tsv >> output_row_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_na-binding.tsv >> output_row_metadata.tmp
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\ttss_window\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"TotalReads",$3,$4,$10,"TotalReads")}}' output_rna-seq.tsv > output_col_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"Rpkm",$3,$4,$10,"Rpkm")}}' output_rna-seq.tsv >> output_col_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"avg_depth",$3,$4,$10,"avg_depth")}}' output_na-binding.tsv >> output_col_metadata.tmp
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv











exit

#	SUMMARY/OUTPUTS
#===============================================================================
printf "\n\nGenerating metrics and formatting overview.md file\n"
printf "\tRead metrics...\n"
reads_processed=$(grep "reads processed" bowtie.log | grep -o "[0-9]\+")
reads_gt0_alignment=$(grep "reads with at least one alignment" bowtie.log | sed -e 's/.*: //')

printf "\tmiRNA metrics...\n"
mirs_total_novel=$(tail -n+2 mirs_novel.tsv | wc -l)
mirs_total_known=$(tail -n+2 mirs_known.tsv | wc -l)
mirs_total_exosome=$(tail -n+2 mirs_known_exocarta_deepmirs.tsv | wc -l)



printf "\tformatting...\n"

printf "## Results Interpretation\n" > overview.md
printf "\n" >> overview.md
printf "#### For the identification of miRDeep2 novel miRNA candidates, the following may be used as a filtering guideline:\n\n" >> overview.md
printf "1. miRDeep score >4\n" >> overview.md
printf "\n" >> overview.md

printf "## ALIGNMENT & miRNA METRICS\n" >> overview.md
printf "-" >> overview.md
printf " Total reads processed: $reads_processed\n" >> overview.md
printf "-" >> overview.md
printf " %s\n" "Reads with at least one alignment: $reads_gt0_alignment" >> overview.md
printf "-" >> overview.md
printf " Novel miRNA detected: $mirs_total_novel\n" >> overview.md
printf "-" >> overview.md
printf " Known miRNA detected: $mirs_total_known\n" >> overview.md
printf "-" >> overview.md
printf " Exosome miRNA detected: $mirs_total_exosome\n" >> overview.md

printf "\n\nWorkflow script run_mirdeep2.sh complete!\n"