#!/bin/bash

#   Shell wrapper for producing morpheus heatmap from genelists from atac/chip/crt (peak) and rna (expression) data
#
##########################################################################################
#
# v1.0.0
# - 20230612, initial release
#
##########################################################################################
printf "$(date)\nLog file for run_genelists.sh\n\n"

# increase stack size from default (8MB) to 32MB (issue with saving from R)
printf "Adjusting ulimit stack size\n"
printf "\tdefault value\n"
R --slave -e 'Cstack_info()["size"]'
printf "\tsetting to 32MB value\n"
ulimit -s 32000
printf "\tcheck new value\n"
R --slave -e 'Cstack_info()["size"]'


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_genelists.sh\`:
Shell wrapper for producing a GCT data file for the morpheus heatmap.
Uses both ATAC/ChIP/CRT (NA [nucleic acid] binding) and RNA-Seq data to derive visualization data.
NA binding data in the form of BAM files per sample is processed to output an average read depth per window +/-5Kbp of each gene's TSS (transcription start site).
RNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.
These data are then integrated into a single count matrix, a row, and a column metadata file as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap viewer.


Primary Output files:
 - heatmap.gct, GCT formatted peak and expression data for morpheus viewer

Secondary Output files:
 - master_samplesheet.tsv, contains formatted information of the input data and files
 - output_row_metadata.tsv, row metadata for GCT formatter
 - output_col_metadata.tsv, column metadata for GCT formatter
 - output_counts.tsv, peak average read depth per TSS window and gene expression counts matrix

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

NOTE:
	All arrays need to be comma separated.

____________________________________________________________________________________________________
References:
 - Tange, O. (2023, May 22). GNU Parallel 20230522 ('Charles'). Zenodo. https://doi.org/10.5281/zenodo.7958356
    
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
		if [[ ${bam_array[name_counter]} != "" ]]; then
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "na-binding" $n ${bam_array[name_counter]}
		fi
		((name_counter++))
	done

	# loop through RNA names
	#	start counter for name array
	name_counter=0
	#	turn bam files into an array
	exp_array=($(echo "$FILES_RNASEQ_EXP" | sed 's/,/\n/g'))
	for n in $(echo "$NAMES_RNASEQ" | sed 's/,/\n/g'); do
		# print formatted samplesheet row
		if [[ ${exp_array[name_counter]} != "" ]]; then
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "rna-seq" $n ${exp_array[name_counter]}
		fi
		((name_counter++))
	done

	((list_counter++))
done > master_samplesheet.tsv


# initialize output files with headers
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tavg_depth\n" > output_na-binding.tmp
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tTotalReads\tRpkm\n" > output_rna-seq.tmp
# function for building tag density window data and heatmap expression data for each experiment type



printf "\n\n"
printf "running master_samplesheet.tsv through simple loop...\n"
while read master; do
	printf "\tget_data(), processing $master\n"
	genelist_number=$(printf "$master" | cut -f1)
	genelist_name=$(printf "$master" | cut -f2)
	genelist_file_tmp=$(printf "$master" | cut -f3)
	# ensure only 1 gene name per row (macs2 will sometimes report 2+ comma separated genes on a single row)
	awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' $genelist_file_tmp > genelist.tsv
	genelist_file="genelist.tsv"
	genelist_annotation_file=$(printf "$master" | cut -f4)
	experiment_type=$(printf "$master" | cut -f5)	# determines data extraction method of file at $6 ("no-binding" or "rna-seq")
	sample_name=$(printf "$master" | cut -f6)
	sample_data=$(printf "$master" | cut -f7)			# either bam or expression table
	window=5000   # this is +/- $window of the TSS (for +strand TSS=$txStart, for -strand TSS=$txEnd), total seq length of plotted area is 2*$window
	total_window_size=$(printf $window | awk '{print($0*2)}')
	window_size=500
	step_size=250
	if [[ $experiment_type == "na-binding" ]]; then
		# user input variables
		bam=$sample_data
		bn=$(basename $bam | sed 's/\..*//')
		samtools sort -@ $THREADS $bam > $bn.sorted.bam
		samtools index -@ $THREADS $bn.sorted.bam
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
			else
				from=$(echo "$txEnd" | awk -v w=$window '{print($0-w)}')
				to=$(echo "$txEnd" | awk -v w=$window '{print($0+w)}')
			fi
			# in case tss is less than the window size away from the start of the sequence, set from to zero
			if [[ "$from" -lt 0 ]]; then from=0; fi
			samtools depth -@ $THREADS -a -r $chr:$from-$to $bn.sorted.bam > $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth
			#   reduce to x windows
			seq $window_size $step_size $total_window_size | while read step; do
				# calculate average depth for window
				avg_depth=$(head -$step $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth | tail -$window_size | awk -F'\t' '{x+=$3}END{printf("%.2f",x/NR)}')
				# print formatted output for na-binding experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\t$step\t$avg_depth\n"
			done >> output_na-binding.tmp
			rm $bn-$timestamp.$geneid.$chr-$txStart-$txEnd.depth
		done
		rm $bn.sorted.bam*
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
				TotalReads=$(printf "$g" | cut -f7)
				Rpkm=$(printf "$g" | cut -f8)
				# print formatted output for rna-seq experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\tNA\t$TotalReads\t$Rpkm\n" >> output_rna-seq.tmp
			fi
		done
	fi
done < master_samplesheet.tsv




# ensure rows are unique
head -1 output_na-binding.tmp > output_na-binding_raw.tsv
head -1 output_rna-seq.tmp > output_rna-seq.tsv
tail -n+2 output_na-binding.tmp | sort | uniq >> output_na-binding_raw.tsv
tail -n+2 output_rna-seq.tmp | sort | uniq >> output_rna-seq.tsv
#rm output_na-binding.tmp
#rm output_rna-seq.tmp







#	HEATMAP 95th percentile
# normalize peak data within each sample and scale from 0-99 per sample (for better visualization)
#	for each sample, find the 95th percentile average depth (pad), and apply normalization by changing:
#		values >= pad to pad value
#		values < 0-pad remain unchanged
printf "\n\n"
printf "running each na-binding data sample through normalization to 95th percentile...\n"
# loop for each sample, grep sample rows, calc percentile, and apply normalization
head -1 output_na-binding_raw.tsv > output_na-binding.tsv
tail -n+2 output_na-binding_raw.tsv | cut -f4 | sort | uniq | while read sample_name; do
	grep "$sample_name" output_na-binding_raw.tsv > peak_norm.tmp
	#	awk explanation:
	#		Sort the file numerically
	#		drop the top 5%
	#		pick the next value
	pad=$(cut -f11 peak_norm.tmp | sort -n | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p=(c/100*5); p=p%1?int(p)+1:p; print a[c-p-1]}')
	#	apply pad normalization and scale from 0-99
	awk -F'\t' -v pad=$pad '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10); if($11 >= pad){printf("%s\n",99)}else{printf("%s\n",99*($11/pad))}}' peak_norm.tmp
done >> output_na-binding.tsv

# scale peak data between 0-99 for better visualization (already done in previous step)
cp output_na-binding.tsv output_na-binding-forheatmap.tsv

# scale RPKM data between 100-199 for better visualization
max=$(tail -n+2 output_rna-seq.tsv | cut -f12 | awk 'BEGIN{max=0};{if ($1 > max) max=$1}END{print max}')
head -1 output_rna-seq.tsv | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12)}' > output_rna-seq-forheatmap.tsv
tail -n+2 output_rna-seq.tsv  | awk -F'\t' -v max=$max '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,(99*($12/max))+100)}' >> output_rna-seq-forheatmap.tsv

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
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"Rpkm",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"avg_depth",$11)}}' output_na-binding-forheatmap.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_na-binding-forheatmap.tsv >> output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq > output_row_metadata.tmp.uniq
# add row number for acast aggregation error for assumed uniqueness
awk -F'\t' '{printf("%s:%s\n",NR,$0)}' output_row_metadata.tmp.uniq >> output_row_metadata.tsv
# update rid in output_counts.tsv file
cp output_counts.tsv output_counts.tmp
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR==FNR){rid[$2]=$1}else{printf("%s:%s\n",rid[$1],$0)}}' <(tail -n+2 output_row_metadata.tsv | cut -f1 | sed 's/:/\t/') <(tail -n+2 output_counts.tmp) >> output_counts.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\ttss_window\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"Rpkm",$3,$4,$10,"Rpkm")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"avg_depth",$3,$4,$10,"avg_depth")}}' output_na-binding-forheatmap.tsv >> output_col_metadata.tmp
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
run_genelists_heatmap.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
# make a copy of R inputs
cp output_row_metadata.tsv output_row_metadata-95p.tsv
cp output_col_metadata.tsv output_col_metadata-95p.tsv
cp output_counts.tsv output_counts-95p.tsv

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF

mv heatmap.html heatmap_peaknorm95.html









#	HEATMAP 99th percentile
# normalize peak data within each sample and scale from 0-99 per sample (for better visualization)
#	for each sample, find the 99th percentile average depth (pad), and apply normalization by changing:
#		values >= pad to pad value
#		values < 0-pad remain unchanged
printf "\n\n"
printf "running each na-binding data sample through normalization to 99th percentile...\n"
# loop for each sample, grep sample rows, calc percentile, and apply normalization
head -1 output_na-binding_raw.tsv > output_na-binding.tsv
tail -n+2 output_na-binding_raw.tsv | cut -f4 | sort | uniq | while read sample_name; do
	grep "$sample_name" output_na-binding_raw.tsv > peak_norm.tmp
	#	awk explanation:
	#		Sort the file numerically
	#		drop the top 1%
	#		pick the next value
	pad=$(cut -f11 peak_norm.tmp | sort -n | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p=(c/100*1); p=p%1?int(p)+1:p; print a[c-p-1]}')
	#	apply pad normalization and scale from 0-99
	awk -F'\t' -v pad=$pad '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10); if($11 >= pad){printf("%s\n",99)}else{printf("%s\n",99*($11/pad))}}' peak_norm.tmp
done >> output_na-binding.tsv

# scale peak data between 0-99 for better visualization (already done in previous step)
cp output_na-binding.tsv output_na-binding-forheatmap.tsv

# scale RPKM data between 100-199 for better visualization
max=$(tail -n+2 output_rna-seq.tsv | cut -f12 | awk 'BEGIN{max=0};{if ($1 > max) max=$1}END{print max}')
head -1 output_rna-seq.tsv | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12)}' > output_rna-seq-forheatmap.tsv
tail -n+2 output_rna-seq.tsv  | awk -F'\t' -v max=$max '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,(99*($12/max))+100)}' >> output_rna-seq-forheatmap.tsv

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
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"Rpkm",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"avg_depth",$11)}}' output_na-binding-forheatmap.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_na-binding-forheatmap.tsv >> output_row_metadata.tmp
sort output_row_metadata.tmp | uniq > output_row_metadata.tmp.uniq
# add row number for acast aggregation error for assumed uniqueness
awk -F'\t' '{printf("%s:%s\n",NR,$0)}' output_row_metadata.tmp.uniq >> output_row_metadata.tsv
# update rid in output_counts.tsv file
cp output_counts.tsv output_counts.tmp
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR==FNR){rid[$2]=$1}else{printf("%s:%s\n",rid[$1],$0)}}' <(tail -n+2 output_row_metadata.tsv | cut -f1 | sed 's/:/\t/') <(tail -n+2 output_counts.tmp) >> output_counts.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\ttss_window\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"Rpkm",$3,$4,$10,"Rpkm")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"avg_depth",$3,$4,$10,"avg_depth")}}' output_na-binding-forheatmap.tsv >> output_col_metadata.tmp
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
run_genelists_heatmap.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
# make a copy of R inputs
cp output_row_metadata.tsv output_row_metadata-99p.tsv
cp output_col_metadata.tsv output_col_metadata-99p.tsv
cp output_counts.tsv output_counts-99p.tsv

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF

mv heatmap.html heatmap_peaknorm99.html








#	HEATMAP no percentile normalization
# normalize peak data within each sample before scaling among ALL samples
#	for each sample, find the 95th percentile average depth (pad), and apply normalization by changing:
#		values >= pad to pad value
#		values < 0-pad remain unchanged
printf "\n\n"
printf "copying na-binding data (no percentile normalization)\n"
# loop for each sample, grep sample rows, calc percentile, and apply normalization
cp output_na-binding_raw.tsv output_na-binding.tsv

# scale peak data between 0-99 for better visualization
#	get max peak value
max=$(tail -n+2 output_na-binding.tsv | cut -f11 | awk 'BEGIN{max=0};{if ($1 > max) max=$1}END{print max}')
head -1 output_na-binding.tsv > output_na-binding-forheatmap.tsv
tail -n+2 output_na-binding.tsv  | awk -F'\t' -v max=$max '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,99*($11/max))}' >> output_na-binding-forheatmap.tsv

# scale RPKM data between 100-199 for better visualization
max=$(tail -n+2 output_rna-seq.tsv | cut -f12 | awk 'BEGIN{max=0};{if ($1 > max) max=$1}END{print max}')
head -1 output_rna-seq.tsv | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12)}' > output_rna-seq-forheatmap.tsv
tail -n+2 output_rna-seq.tsv  | awk -F'\t' -v max=$max '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,(99*($12/max))+100)}' >> output_rna-seq-forheatmap.tsv

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
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"Rpkm",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$3,$4,$10,"avg_depth",$11)}}' output_na-binding-forheatmap.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$1,$2,$5,$6,$7,$8,$9)}}' output_na-binding-forheatmap.tsv >> output_row_metadata.tmp
sort output_row_metadata.tmp | uniq > output_row_metadata.tmp.uniq
# add row number for acast aggregation error for assumed uniqueness
awk -F'\t' '{printf("%s:%s\n",NR,$0)}' output_row_metadata.tmp.uniq >> output_row_metadata.tsv
# update rid in output_counts.tsv file
cp output_counts.tsv output_counts.tmp
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR==FNR){rid[$2]=$1}else{printf("%s:%s\n",rid[$1],$0)}}' <(tail -n+2 output_row_metadata.tsv | cut -f1 | sed 's/:/\t/') <(tail -n+2 output_counts.tmp) >> output_counts.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\ttss_window\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"Rpkm",$3,$4,$10,"Rpkm")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s\t%s\t%s\t%s\t%s\n",$3,$4,$10,"avg_depth",$3,$4,$10,"avg_depth")}}' output_na-binding-forheatmap.tsv >> output_col_metadata.tmp
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
run_genelists_heatmap.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
# make a copy of R inputs
cp output_row_metadata.tsv output_row_metadata-100p.tsv
cp output_col_metadata.tsv output_col_metadata-100p.tsv
cp output_counts.tsv output_counts-100p.tsv

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF

printf "\n\nWorkflow script run_genelists.sh complete!\n"