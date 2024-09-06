#!/bin/bash

#   Shell wrapper for producing morpheus heatmap from genelists from rna (expression) data only. Adapted from the combined version (run_genelists.sh)
#
##########################################################################################
#
# v1.0.0
# - 20240603, initial release
#
##########################################################################################
printf "$(date)\nLog file for run_genelists_rnaseq.sh\n\n"

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
Help message for \`run_genelists_rnaseq.sh\`:
Shell wrapper for producing a GCT data file for the morpheus heatmap.
Uses only RNA-Seq data to derive visualization data.
RNA-Seq data in the form of gene expression count matrices are processed to output TotalReads and Rpkm values per gene.
Additionally, TotalReads are used to produce a separate count matrix of VST values, and VST z-scores.
This data is then integrated into a single count matrix (one for scaled, a separate matrix for VST, and another for VST z-scores), a row, and a column metadata file.
These are used as input to an Rscript that will format the 3 files into GCT format for morpheus heatmap generation.


  Primary Output files:
  - heatmap_TotalReads.html, html of morpheus heatmap with preconfigured settings, TotalReads, no data scaling
  - heatmap_vst.html, html of morpheus heatmap with preconfigured settings, VST values, no data scaling
  - heatmap_vst_zscore.html, html of morpheus heatmap with preconfigured settings, Z-scores of VST values, no data scaling
  - heatmap_Rpkm.html, html of morpheus heatmap with preconfigured settings, RPKM, no data scaling
  - heatmap_scaled100.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, no percentile cutoff
  - heatmap_scaled99.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled 0-99, max value set to 99th percentile RPKM value
  - heatmap_scaled95.html, html of morpheus heatmap with preconfigured settings, RPKM, data scaled  0-99, max value set to 95th percentile RPKM value

  Secondary Output files:
  - master_samplesheet_scaled.tsv, contains formatted information of the input data and files for raw and scaled heatmaps
  - master_samplesheet_vst.tsv, contains formatted information of the input data and files for vst normalized heatmaps
  - master_samplesheet_vst_zscore.tsv, contains formatted information of the input data and files for vst normalized heatmaps

PARAMS:
    SECTION 1: general
	-h	help		show this message
	-t  INT			number of threads
	-a	ARRAY		array of genelist sample names (no commas in names)
	-b  FILE ARRAY	array of associated annotation files for each gene list from (-c), with header
	-c  FILE ARRAY	array of filtered gene list TSVs (must be headerless, columns are: chr, txStart, txEnd, geneID, L2FC, Strand)
	-e	ARARY		array of sample names from RNA-Seq experiments (no commas in names)
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
while getopts "ht:a:b:c:e:g:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		a) GENELIST_NAMES=$OPTARG ;;
		b) GENELIST_ANNOTATION_FILES=$OPTARG ;;
		c) GENELIST_FILTERED_FILES=$OPTARG ;;
		e) NAMES_RNASEQ=$OPTARG ;;
		g) FILES_RNASEQ_EXP=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$GENELIST_NAMES" == "" ]]; then echo "error: required param missing (-a), array of genelist sample names"; exit; fi
if [[ "$GENELIST_ANNOTATION_FILES" == "" ]]; then echo "error: required param missing (-b), array of genelist annotation (feature_file.tsv) files"; exit; fi
if [[ "$GENELIST_FILTERED_FILES" == "" ]]; then echo "error: required param missing (-c), array of genelist filtered (feature_file) files"; exit; fi
if [[ "$NAMES_RNASEQ" == "" ]]; then echo "warning: optional param missing (-e), array of RNA-Seq sample names"; fi
if [[ "$FILES_RNASEQ_EXP" == "" ]]; then echo "warning: optional param missing (-g), array of RNA-Seq expression files"; fi


# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-a, \$GENELIST_NAMES, $GENELIST_NAMES\n"
printf "\t-b, \$GENELIST_ANNOTATION_FILES, $GENELIST_ANNOTATION_FILES\n"
printf "\t-c, \$GENELIST_FILTERED_FILES, $GENELIST_FILTERED_FILES\n"
printf "\t-e, \$NAMES_RNASEQ, $NAMES_RNASEQ\n"
printf "\t-g, \$FILES_RNASEQ_EXP, $FILES_RNASEQ_EXP\n"
printf "\n\n"


#	MAIN
#===============================================================================




#########################################################################################
#					HEATMAPS FROM VST COUNT DATA										#
#					EXPRESSION VALUES not scaled										#
#########################################################################################

printf "generating count matrix from all input sample expression data for VST\n"
#	make list of all unique genes among all input samples expression tables
#oIFS=$IFS
IFS=$','
printf "\tmaking unique geneid list\n"
for n in ${FILES_RNASEQ_EXP[@]}; do tail -n+2 $n | cut -f2; done | sort | uniq > uniq.geneids.tsv
#	make matrix
#		initialize with headers
printf "\tmaking TotalReads count matrix\n"
printf "geneid" > expression_matrix_totalreads.tsv
#for n in ${NAMES_RNASEQ[@]}; do printf "\t$n"; done >> expression_matrix_totalreads.tsv
for n in ${NAMES_RNASEQ[@]}; do x=$(echo $n | sed 's/ /_/g'); printf "\t$x"; done >> expression_matrix_totalreads.tsv
printf "\n" >> expression_matrix_totalreads.tsv
#		loop through each gene, pull TotalReads from each sample expression tsv
while read geneid; do printf "%s" $geneid; for n in ${FILES_RNASEQ_EXP[@]}; do printf "\t%s" $(grep -m 1 "`printf '\t'`$geneid`printf '\t'`" ${n} | cut -f7); done; printf "\n"; done < uniq.geneids.tsv >> expression_matrix_totalreads.tsv
# 		reset IFS
IFS=$' \t\n'
#IFS=$oIFS

# run Rscript to produce variance stabilized transformed (VST) count data matrix (vst_counts_matrix.tsv) and table (vst_counts_table.tsv)
#		also produces z-score of VST matrix (vst_z-score_matrix.tsv) and table (vst_zscore_table.tsv)
printf "\trunning vst normalization\n"
Rscript /usr/local/bin/run_vst_norm.R expression_matrix_totalreads.tsv






###		VST heatmap

printf "generating master samplesheet from all inputs using VST counts\n"
#	start counter for genelist arrays (needs to start at 0 for getting correct values of names and annotations array indices)
list_counter=0
# 	set "internal field separator" to use comma on array string 'itemSeparator: ","'inputs from cwl (default is IFS=$' \t\n')
IFS=$','
#	turn genelist annotation files and genelist names into arrays
annotations_array=($(echo "$GENELIST_ANNOTATION_FILES"))
names_array=($(echo "$GENELIST_NAMES"))
for f in ${GENELIST_FILTERED_FILES[@]}; do
	# loop through RNA names
	#	start counter for name array
	name_counter=0
	#	turn expression tsv files into an array
	exp_array=($(echo "$FILES_RNASEQ_EXP"))
	for n in ${NAMES_RNASEQ[@]}; do
		# need to replace spaces with underscores in sample names before grepping the sample name in table output from 'run_vst_norm.R'
		x=$(echo $n | sed 's/ /_/g')
        # replace TotalReads and RPKM with VST values, make new 'exp_array' file (genes.tsv) for samplesheet
		#		replace or match geneids with other metadata from each individual sample file, while making the master samplesheet
		#		these are then used in place of the 'exp_array=($(echo "$FILES_RNASEQ_EXP"))' in this section's 'master samplesheet')
        printf "RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tVST\n" > "sample${name_counter}-vst_"$(basename ${exp_array[name_counter]})
        awk -F'\t' '{if(NR==FNR){vst[$1]=$3}else{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,vst[$2])}}' <(grep "${x}" vst_counts_table.tsv) <(tail -n+2 ${exp_array[name_counter]}) >> "sample${name_counter}-vst_"$(basename ${exp_array[name_counter]})
		# print formatted samplesheet row
		if [[ ${exp_array[name_counter]} != "" ]]; then
			# 20231102 - ensure sample name uniqueness, possible fix for duplicate rows > acast > aggregate length issue
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "rna-seq" ${name_counter}_${n} "sample${name_counter}-vst_"$(basename ${exp_array[name_counter]})
		fi
		((name_counter++))
	done
	((list_counter++))
done > master_samplesheet_vst.tsv
# reset IFS
IFS=$' \t\n'




# initialize output files with headers
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tVST\n" > output_rna-seq.tmp
printf "\n\n"
printf "\trunning master_samplesheet_vst.tsv through loop, organizing and filtering expression data...\n"
while read master; do
	printf "\t\tprocessing $master\n"
	genelist_number=$(printf "$master" | cut -f1)
	genelist_name=$(printf "$master" | cut -f2)
	genelist_file_tmp=$(printf "$master" | cut -f3)
	# ensure only 1 gene name per row (iaintersect can sometimes report 2+ comma separated genes on a single row)
	awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' $genelist_file_tmp > genelist.tmp
	# also ensure only a single geneid per chr (different chr can share a geneid)
	cut -f1,4 genelist.tmp | sort | uniq | while read chr_geneid; do
		chr=$(printf "$chr_geneid" | cut -f1)
		geneid=$(printf "$chr_geneid" | cut -f2)
		grep "$chr" genelist.tmp | grep -m1 "$geneid"
	done > genelist.tsv
	genelist_file="genelist.tsv"
	genelist_annotation_file=$(printf "$master" | cut -f4)
	experiment_type=$(printf "$master" | cut -f5)	# determines data extraction method of file at $6 ("na-binding" or "rna-seq" - all in this script should be "rna-seq")
	sample_name=$(printf "$master" | cut -f6)
	sample_data=$(printf "$master" | cut -f7)		# should all be *.genes.tsv files (vst expression data)

	if [[ $experiment_type == "rna-seq" ]]; then
		# user input variables
		exp=$sample_data
		bn=$(basename $exp | sed 's/\..*//')
		# for each gene in filtered genelist file, get VST values
		cat $genelist_file | while read filtered_gene; do
			printf "\t\trna-seq, processing gene:\t $filtered_gene\n"
			chr=$(printf "$filtered_gene" | cut -f1)  # genelist col1
			txStart=$(printf "$filtered_gene" | cut -f2)  # genelist col2
			txEnd=$(printf "$filtered_gene" | cut -f3)  # genelist col3
			geneid=$(printf "$filtered_gene" | cut -f4)  # genelist col4
			strand=$(printf "$filtered_gene" | cut -f6)  # genelist col6
			# search for each gene, only use $chr and $geneid due to peak annotations having multiple geneids per line (when split coordinates for each are not accurate)
			#	this also addresses the same geneid on multiple chromosomes without using exact coordinates, since the nearest peaks to multiple genes much be co-localized (on same chr)
			g=$(grep -P -m 1 "\t$geneid\t$chr\t" $exp)
			if [[ $g == "" ]]; then
				printf "\t\t\tWARNING, gene does not exist in expression data file: $geneid"
			else
				VST=$(printf "$g" | cut -f7)
				# print formatted output for rna-seq experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\tNA\t$VST\n" >> output_rna-seq.tmp
			fi
		done
	fi
done < master_samplesheet_vst.tsv

cp output_rna-seq.tmp output_rna-seq-vst.tmp






# ensure rows are unique
head -1 output_rna-seq.tmp > output_rna-seq_raw.tsv
tail -n+2 output_rna-seq.tmp | sort | uniq >> output_rna-seq_raw.tsv


# 	CLUSTERING by HOPACH method for VST data, for each genelist
printf "\n\n"
printf "\trunning hopach clustering per genelist\n"
#		apply same rank orders (of each list and data type) to row metadata file

# expression data clustering
#	run cluster per gene list
printf "\t\tdata type: rna-seq\n"
#	print header for new rna-seq data file for heatmaps
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tVST\thopach_rank\n" > output_rna-seq-cluster_data.tmp
#	make another genelist counter, but start at 1, not 0, as leading zeros will be stripped in R after prepending to the cluster rank
list_counter=1
tail -n+2 output_rna-seq_raw.tsv | cut -f1 | sort | uniq | while read genelist_number; do
	head -1 output_rna-seq_raw.tsv > ${genelist_number}-cluster_data.tmp
	grep "$genelist_number" output_rna-seq_raw.tsv >> ${genelist_number}-cluster_data.tmp
	printf "\t\t\tclustering for $genelist_number\n"
	R CMD /usr/local/bin/run_hopach_clustering.R ${genelist_number}-cluster_data.tmp "VST" 1> hopach_stdout.log 2> hopach_stderr.log
	# save a copy of hopach output
	cp hopach_results.out hopach_results.out-rnaseq-${genelist_number}
	# use col2 "UID" (geneid) to add the rank order from col7 "Final.Level.Order" to both ${genelist_number}-cluster_data.tmp
	#	for each sample, each gene should have the same rank order value
	#	in awk, add zero padding up to 5 digits, very unlikely to have >99999 genes in a list
	awk -F'\t' -v listnumber=$list_counter '{if(NR==FNR){rank_order[$2]=$7}else{if(rank_order[$5]!=""){printf("%s\t%s%05d\n",$0,listnumber,rank_order[$5])}}}' <(tail -n+2 hopach_results.out | sed 's/"//g') ${genelist_number}-cluster_data.tmp >> output_rna-seq-cluster_data.tmp
	((list_counter++))
done





# setup heatmap file variables
data_rnaseq="output_rna-seq-cluster_data.tmp"			# expression rank col12 (for vst)
cp output_rna-seq-cluster_data.tmp output_rna-seq-forheatmap.tsv

#	HEATMAP no percentile normalization
printf "\n\n"
printf "make gct files for vst heatmap generation\n"

# merge outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window		<- (REMOVING FOR EXP-ONLY)
#		row9 - data_type (vst)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (vst)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv


# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-2colors.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF
mv heatmap.html heatmap_vst.html







###		VST Z-score heatmap

printf "generating master samplesheet from all inputs using VST z-scores\n"
#	start counter for genelist arrays (needs to start at 0 for getting correct values of names and annotations array indices)
list_counter=0
# 	set "internal field separator" to use comma on array string 'itemSeparator: ","'inputs from cwl (default is IFS=$' \t\n')
IFS=$','
#	turn genelist annotation files and genelist names into arrays
annotations_array=($(echo "$GENELIST_ANNOTATION_FILES"))
names_array=($(echo "$GENELIST_NAMES"))
for f in ${GENELIST_FILTERED_FILES[@]}; do
	# loop through RNA names
	#	start counter for name array
	name_counter=0
	#	turn expression tsv files into an array
	exp_array=($(echo "$FILES_RNASEQ_EXP"))
	for n in ${NAMES_RNASEQ[@]}; do
		# need to replace spaces with underscores in sample names before grepping the sample name in table output from 'run_vst_norm.R'
		x=$(echo $n | sed 's/ /_/g')
        # replace TotalReads and RPKM with VST z-scores, make new 'exp_array' file (genes.tsv) for samplesheet
		#		replace or match geneids with other metadata from each individual sample file, while making the master samplesheet
		#		these are then used in place of the 'exp_array=($(echo "$FILES_RNASEQ_EXP"))' in this section's 'master samplesheet')
        printf "RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tzscore\n" > "sample${name_counter}-vst_zscore_"$(basename ${exp_array[name_counter]})
        awk -F'\t' '{if(NR==FNR){vzs[$1]=$3}else{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,vzs[$2])}}' <(grep "${x}" vst_zscore_table.tsv) <(tail -n+2 ${exp_array[name_counter]}) >> "sample${name_counter}-vst_zscore_"$(basename ${exp_array[name_counter]})
		# print formatted samplesheet row
		if [[ ${exp_array[name_counter]} != "" ]]; then
			# 20231102 - ensure sample name uniqueness, possible fix for duplicate rows > acast > aggregate length issue
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "rna-seq" ${name_counter}_${n} "sample${name_counter}-vst_zscore_"$(basename ${exp_array[name_counter]})
		fi
		((name_counter++))
	done
	((list_counter++))
done > master_samplesheet_vst_zscore.tsv
# reset IFS
IFS=$' \t\n'




# initialize output files with headers
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tzscore\n" > output_rna-seq.tmp
printf "\n\n"
printf "\trunning master_samplesheet_vst_zscore.tsv through loop, organizing and filtering expression data...\n"
while read master; do
	printf "\t\tprocessing $master\n"
	genelist_number=$(printf "$master" | cut -f1)
	genelist_name=$(printf "$master" | cut -f2)
	genelist_file_tmp=$(printf "$master" | cut -f3)
	# ensure only 1 gene name per row (iaintersect can sometimes report 2+ comma separated genes on a single row)
	awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' $genelist_file_tmp > genelist.tmp
	# also ensure only a single geneid per chr (different chr can share a geneid)
	cut -f1,4 genelist.tmp | sort | uniq | while read chr_geneid; do
		chr=$(printf "$chr_geneid" | cut -f1)
		geneid=$(printf "$chr_geneid" | cut -f2)
		grep "$chr" genelist.tmp | grep -m1 "$geneid"
	done > genelist.tsv
	genelist_file="genelist.tsv"
	genelist_annotation_file=$(printf "$master" | cut -f4)
	experiment_type=$(printf "$master" | cut -f5)	# determines data extraction method of file at $6 ("na-binding" or "rna-seq" - all in this script should be "rna-seq")
	sample_name=$(printf "$master" | cut -f6)
	sample_data=$(printf "$master" | cut -f7)		# should all be *.genes.tsv files (vst z-score data)

	if [[ $experiment_type == "rna-seq" ]]; then
		# user input variables
		exp=$sample_data
		bn=$(basename $exp | sed 's/\..*//')
		# for each gene in filtered genelist file, get VST z-scores
		cat $genelist_file | while read filtered_gene; do
			printf "\t\trna-seq, processing gene:\t $filtered_gene\n"
			chr=$(printf "$filtered_gene" | cut -f1)  # genelist col1
			txStart=$(printf "$filtered_gene" | cut -f2)  # genelist col2
			txEnd=$(printf "$filtered_gene" | cut -f3)  # genelist col3
			geneid=$(printf "$filtered_gene" | cut -f4)  # genelist col4
			strand=$(printf "$filtered_gene" | cut -f6)  # genelist col6
			# search for each gene, only use $chr and $geneid due to peak annotations having multiple geneids per line (when split coordinates for each are not accurate)
			#	this also addresses the same geneid on multiple chromosomes without using exact coordinates, since the nearest peaks to multiple genes much be co-localized (on same chr)
			g=$(grep -P -m 1 "\t$geneid\t$chr\t" $exp)
			if [[ $g == "" ]]; then
				printf "\t\t\tWARNING, gene does not exist in expression data file: $geneid"
			else
				zscore=$(printf "$g" | cut -f7)
				# print formatted output for rna-seq experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\tNA\t$zscore\n" >> output_rna-seq.tmp
			fi
		done
	fi
done < master_samplesheet_vst_zscore.tsv

cp output_rna-seq.tmp output_rna-seq-vst-zscore.tmp






# ensure rows are unique
head -1 output_rna-seq.tmp > output_rna-seq_raw.tsv
tail -n+2 output_rna-seq.tmp | sort | uniq >> output_rna-seq_raw.tsv




# 	CLUSTERING by HOPACH method for VST z-scores, for each genelist
printf "\n\n"
printf "\trunning hopach clustering per genelist\n"
#		apply same rank orders (of each list and data type) to row metadata file

# expression data clustering
#	run cluster per gene list
printf "\t\tdata type: rna-seq\n"
#	print header for new rna-seq data file for heatmaps
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tzscore\thopach_rank\n" > output_rna-seq-cluster_data.tmp
#	make another genelist counter, but start at 1, not 0, as leading zeros will be stripped in R after prepending to the cluster rank
list_counter=1
tail -n+2 output_rna-seq_raw.tsv | cut -f1 | sort | uniq | while read genelist_number; do
	head -1 output_rna-seq_raw.tsv > ${genelist_number}-cluster_data.tmp
	grep "$genelist_number" output_rna-seq_raw.tsv >> ${genelist_number}-cluster_data.tmp
	printf "\t\t\tclustering for $genelist_number\n"
	R CMD /usr/local/bin/run_hopach_clustering.R ${genelist_number}-cluster_data.tmp "zscore" 1> hopach_stdout.log 2> hopach_stderr.log
	# save a copy of hopach output
	cp hopach_results.out hopach_results.out-rnaseq-${genelist_number}
	# use col2 "UID" (geneid) to add the rank order from col7 "Final.Level.Order" to both ${genelist_number}-cluster_data.tmp
	#	for each sample, each gene should have the same rank order value
	#	in awk, add zero padding up to 5 digits, very unlikely to have >99999 genes in a list
	awk -F'\t' -v listnumber=$list_counter '{if(NR==FNR){rank_order[$2]=$7}else{if(rank_order[$5]!=""){printf("%s\t%s%05d\n",$0,listnumber,rank_order[$5])}}}' <(tail -n+2 hopach_results.out | sed 's/"//g') ${genelist_number}-cluster_data.tmp >> output_rna-seq-cluster_data.tmp
	((list_counter++))
done





# setup heatmap file variables
data_rnaseq="output_rna-seq-cluster_data.tmp"			# expression rank col12 (for vst z-score)
cp output_rna-seq-cluster_data.tmp output_rna-seq-forheatmap.tsv

#	HEATMAP no percentile normalization
printf "\n\n"
printf "make gct files for vst z-score heatmap generation\n"

# merge outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window		<- (REMOVING FOR EXP-ONLY)
#		row9 - data_type (vst z-score)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (vst z-score)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv


# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-3colors.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF
mv heatmap.html heatmap_vst_zscore.html












#########################################################################################
#					HEATMAPS FROM TOTALREADS & RPKM COUNT DATA							#
#					SCALED EXPRESSION VALUES 0-99										#
#########################################################################################

printf "generating master samplesheet from all inputs using TotalReads and RPKM counts\n"
#	start counter for genelist arrays (needs to start at 0 for getting correct values of names and annotations array indices)
list_counter=0
# 	set "internal field separator" to use comma on array string 'itemSeparator: ","'inputs from cwl (default is IFS=$' \t\n')
IFS=$','
#	turn genelist annotation files and genelist names into arrays
annotations_array=($(echo "$GENELIST_ANNOTATION_FILES"))
names_array=($(echo "$GENELIST_NAMES"))
for f in ${GENELIST_FILTERED_FILES[@]}; do
	# loop through RNA names
	#	start counter for name array
	name_counter=0
	#	turn expression tsv files into an array
	exp_array=($(echo "$FILES_RNASEQ_EXP"))
	for n in ${NAMES_RNASEQ[@]}; do
		# print formatted samplesheet row
		if [[ ${exp_array[name_counter]} != "" ]]; then
			# 20231102 - ensure sample name uniqueness, possible fix for duplicate rows > acast > aggregate length issue
			printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "genelist_${list_counter}" ${names_array[list_counter]} $f ${annotations_array[list_counter]} "rna-seq" ${name_counter}_${n} ${exp_array[name_counter]}
		fi
		((name_counter++))
	done
	((list_counter++))
done > master_samplesheet_scaled.tsv
# reset IFS
IFS=$' \t\n'




# initialize output files with headers
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tTotalReads\tRpkm\n" > output_rna-seq.tmp
printf "\n\n"
printf "\trunning master_samplesheet_scaled.tsv through loop, organizing and filtering expression data...\n"
while read master; do
	printf "\t\tprocessing $master\n"
	genelist_number=$(printf "$master" | cut -f1)
	genelist_name=$(printf "$master" | cut -f2)
	genelist_file_tmp=$(printf "$master" | cut -f3)
	# ensure only 1 gene name per row (iaintersect can sometimes report 2+ comma separated genes on a single row)
	awk -F'\t' '{split($4,col4,","); for(i in col4){printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,col4[i],$5,$6)}}' $genelist_file_tmp > genelist.tmp
	# also ensure only a single geneid per chr (different chr can share a geneid)
	cut -f1,4 genelist.tmp | sort | uniq | while read chr_geneid; do
		chr=$(printf "$chr_geneid" | cut -f1)
		geneid=$(printf "$chr_geneid" | cut -f2)
		grep "$chr" genelist.tmp | grep -m1 "$geneid"
	done > genelist.tsv
	genelist_file="genelist.tsv"
	genelist_annotation_file=$(printf "$master" | cut -f4)
	experiment_type=$(printf "$master" | cut -f5)	# determines data extraction method of file at $6 ("na-binding" or "rna-seq" - all in this script should be "rna-seq")
	sample_name=$(printf "$master" | cut -f6)
	sample_data=$(printf "$master" | cut -f7)		# should all be *.genes.tsv files (expression data)

	if [[ $experiment_type == "rna-seq" ]]; then
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
			# search for each gene, only use $chr and $geneid due to peak annotations having multiple geneids per line (when split coordinates for each are not accurate)
			#	this also addresses the same geneid on multiple chromosomes without using exact coordinates, since the nearest peaks to multiple genes much be co-localized (on same chr)
			g=$(grep -P -m 1 "\t$geneid\t$chr\t" $exp)
			if [[ $g == "" ]]; then
				printf "\t\t\tWARNING, gene does not exist in expression data file: $geneid"
			else
				TotalReads=$(printf "$g" | cut -f7)
				Rpkm=$(printf "$g" | cut -f8)
				# print formatted output for rna-seq experiments
				printf "$genelist_number\t$genelist_name\t$experiment_type\t$sample_name\t$geneid\t$chr\t$txStart\t$txEnd\t$strand\tNA\t$TotalReads\t$Rpkm\n" >> output_rna-seq.tmp
			fi
		done
	fi
done < master_samplesheet_scaled.tsv




# ensure rows are unique
head -1 output_rna-seq.tmp > output_rna-seq_raw.tsv
tail -n+2 output_rna-seq.tmp | sort | uniq >> output_rna-seq_raw.tsv




# 	CLUSTERING by HOPACH method for Rpkm expression data, for each genelist
printf "\n\n"
printf "\trunning hopach clustering per genelist on RPKM values\n"
#		apply same rank orders (of each list and data type) to row metadata file

# Rpkm clustering
#	run cluster per gene list
printf "\t\tdata type: rna-seq\n"
#	print header for new rna-seq data file for heatmaps
printf "genelist_number\tgenelist_name\texperiment_type\tsample_name\tgeneid\tchr\ttxStart\ttxEnd\tstrand\ttss_window\tTotalReads\tRpkm\thopach_rank\n" > output_rna-seq-cluster_data.tmp
#	make another genelist counter, but start at 1, not 0, as leading zeros will be stripped in R after prepending to the cluster rank
list_counter=1
tail -n+2 output_rna-seq_raw.tsv | cut -f1 | sort | uniq | while read genelist_number; do
	head -1 output_rna-seq_raw.tsv > ${genelist_number}-cluster_data.tmp
	grep "$genelist_number" output_rna-seq_raw.tsv >> ${genelist_number}-cluster_data.tmp
	printf "\t\t\tclustering for $genelist_number\n"
	R CMD /usr/local/bin/run_hopach_clustering.R ${genelist_number}-cluster_data.tmp "Rpkm" 1> hopach_stdout.log 2> hopach_stderr.log
	# save a copy of hopach output
	cp hopach_results.out hopach_results.out-rnaseq-${genelist_number}
	# use col2 "UID" (geneid) to add the rank order from col7 "Final.Level.Order" to both ${genelist_number}-cluster_data.tmp
	#	for each sample, each gene should have the same rank order value
	#	in awk, add zero padding up to 5 digits, very unlikely to have >99999 genes in a list
	awk -F'\t' -v listnumber=$list_counter '{if(NR==FNR){rank_order[$2]=$7}else{if(rank_order[$5]!=""){printf("%s\t%s%05d\n",$0,listnumber,rank_order[$5])}}}' <(tail -n+2 hopach_results.out | sed 's/"//g') ${genelist_number}-cluster_data.tmp >> output_rna-seq-cluster_data.tmp
	((list_counter++))
done










###		TotalReads heatmap
awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13)}' $data_rnaseq > output_rna-seq.tsv
cp output_rna-seq.tsv output_rna-seq-forheatmap.tsv

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-2colors.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
# make a copy of R inputs
cp output_row_metadata.tsv output_row_metadata-TotalReads.tsv
cp output_col_metadata.tsv output_col_metadata-TotalReads.tsv
cp output_counts.tsv output_counts-TotalReads.tsv

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF
mv heatmap.html heatmap_TotalReads.html




###		RPKM heatmap
awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$13)}' $data_rnaseq > output_rna-seq.tsv
cp output_rna-seq.tsv output_rna-seq-forheatmap.tsv

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv

# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-2colors.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
# make a copy of R inputs
cp output_row_metadata.tsv output_row_metadata-Rpkm.tsv
cp output_col_metadata.tsv output_col_metadata-Rpkm.tsv
cp output_counts.tsv output_counts-Rpkm.tsv

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF
mv heatmap.html heatmap_Rpkm.html







###		SCALING (100, 99, and 95 percentile heatmaps)

# setup heatmap file variables
data_rnaseq="output_rna-seq-cluster_data.tmp"			# Rpkm hopach rank col13



#	HEATMAP no percentile normalization
# still scaling among ALL samples
printf "\n\n"
printf "scaling expression data, no percentile normalization\n"

# scale RPKM data between 0-99 for better visualization
max=$(tail -n+2 $data_rnaseq | cut -f12 | awk 'BEGIN{max=0};{if ($1 > max) max=$1}END{print max}')
head -1 $data_rnaseq | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$13)}' > output_rna-seq-forheatmap.tsv
tail -n+2 $data_rnaseq  | awk -F'\t' -v max=$max '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,99*($12/max),$13)}' >> output_rna-seq-forheatmap.tsv
cp output_rna-seq-forheatmap.tsv output_rna-seq-100.tsv


# merge outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window		<- (REMOVING FOR EXP-ONLY)
#		row9 - data_type (avg_depth, TotalReads, or Rpkm)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (avg_depth, TotalReads, Rpkm)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv


# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-scaled.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./

# inject javascript to configure the heatmap
ed heatmap.html <<EOF
/^<script>(function(global){"use strict";var morpheus=typeof morpheus!=="undefined"?morpheus:{}
-2
a
setTimeout( function() { let groupByBtn = document.querySelector('div.btn-group.bootstrap-select.show-tick.form-control button[data-toggle="dropdown"]'); groupByBtn.click(); let groupRowSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let geneListOption = groupRowSelectionOptions.filter(function (el) { return el.textContent === 'genelist_name' })[0]; geneListOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); setTimeout( () => { let toolsBtn = document.getElementById("morpheus4"); toolsBtn.click(); let sortGroupBtn = document.querySelector('[data-action="Sort/Group"]'); sortGroupBtn.click(); document.querySelector('input[name="rowsOrColumns"][value="columns"]').click(); groupByBtn.click(); let groupColSelectionOptions = Array.from(document.querySelectorAll('ul.dropdown-menu.inner li a[role="option"] span.text')); let samplenameOption = groupColSelectionOptions.filter(function (el) { return el.textContent === 'sample_name' })[0]; samplenameOption.click(); let toolConfirmationBtn = document.getElementsByClassName('modal-footer')[0].querySelector('[name="ok"]'); toolConfirmationBtn.click(); }, 250); }, 1000);
.
wq
EOF
mv heatmap.html heatmap_scaled100.html





#	HEATMAP 99th percentile

# normalize expression data within each sample and scale from 0-99 per sample (for better visualization)
#	for each sample, find the 99th percentile average depth (pad), and apply normalization by changing:
#		values >= pad to pad value
#		values < 0-pad remain unchanged
printf "\n\n"
printf "running each rna-seq data sample through normalization to 99th percentile...\n"
# loop for each sample, grep sample rows, calc percentile, and apply normalization
#	in header also dropping Total_reads [$11] column (was not in nabinding data)
head -1 $data_rnaseq | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$13)}' > output_rna-seq.tsv
tail -n+2 $data_rnaseq | cut -f4 | sort | uniq | while read sample_name; do
	awk -F'\t' -v sample_name="$sample_name" '{if($4==sample_name){print($0)}}' $data_rnaseq > rna_norm.tmp
	#	awk explanation:
	#		Sort the file numerically (on Rpkm col12)
	#		drop the top 1%
	#		pick the next value
	pad=$(cut -f12 rna_norm.tmp | sort -n | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p=(c/100*1); p=p%1?int(p)+1:p; print a[c-p-1]}')
	#	apply pad normalization and scale from 0-99 (also dropping Total_reads column here)
	awk -F'\t' -v pad=$pad '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10); if($12 >= pad){printf("%s\t",99)}else{printf("%s\t",99*($12/pad))}; printf("%s\n",$13)}' rna_norm.tmp
done >> output_rna-seq.tsv
cp output_rna-seq.tsv output_rna-seq-forheatmap.tsv


# merge outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window		<- (REMOVING FOR EXP-ONLY)
#		row9 - data_type (avg_depth, TotalReads, or Rpkm)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (avg_depth, TotalReads, Rpkm)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv


# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-scaled.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
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
mv heatmap.html heatmap_scaled99.html





#	HEATMAP 95th percentile data

# normalize expression data within each sample and scale from 0-99 per sample (for better visualization)
#	for each sample, find the 95th percentile average depth (pad) [or percentile rpkm], and apply normalization by changing:
#		values >= pad to pad value
#		values < 0-pad remain unchanged
printf "\n\n"
printf "running each rna-seq data sample through normalization to 95th percentile...\n"
# loop for each sample, grep sample rows, calc percentile, and apply normalization
#	in header also dropping Total_reads [$11] column (was not in nabinding data)
head -1 $data_rnaseq | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$13)}' > output_rna-seq.tsv
tail -n+2 $data_rnaseq | cut -f4 | sort | uniq | while read sample_name; do
	awk -F'\t' -v sample_name="$sample_name" '{if($4==sample_name){print($0)}}' $data_rnaseq > rna_norm.tmp
	#	awk explanation:
	#		Sort the file numerically (on Rpkm col12)
	#		drop the top 5%
	#		pick the next value
	pad=$(cut -f12 rna_norm.tmp | sort -n | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p=(c/100*5); p=p%1?int(p)+1:p; print a[c-p-1]}')
	#	apply pad normalization and scale from 0-99 (also dropping Total_reads column here)
	awk -F'\t' -v pad=$pad '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10); if($12 >= pad){printf("%s\t",99)}else{printf("%s\t",99*($12/pad))}; printf("%s\n",$13)}' rna_norm.tmp
done >> output_rna-seq.tsv
cp output_rna-seq.tsv output_rna-seq-95.tsv


# merge outputs into GCT formatted file for morpheus heatmap compatibility
#	needs 3 files: counts matrix, row metadata, column metadata

# GCT format:
# line 1, version string "#1.3"
# line 2, col1 (# rows in matrix), col2 (# cols in matrix), col3 (# columns of row metadata), col4 (# rows of col metadata)
# col metadata (row4-a), sample metadata
#		row4 - genelist_number
#		row5 - genelist_name
#		row6 - experiment_type
#		row7 - sample_name
#		row8 - tss_window		<- (REMOVING FOR EXP-ONLY)
#		row9 - data_type (avg_depth, TotalReads, or Rpkm)
# row metadata (col2-b), gene metadata
#		col2 - chr (chr[n+])
#		col3 - txStart (int)
#		col4 - txEnd (int)
#		col5 - strand (+ or -)

# make unique row and column names with "value" as (avg_depth, TotalReads, Rpkm)

# count matrix data (rows (a+1)-x, cols (b+1)-y)
printf "rid\tcid\tvalue\n" > output_counts.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s:%s:%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$3,$4,"expression",$11)}}' output_rna-seq-forheatmap.tsv >> output_counts.tsv


# row metadata file
printf "rid\tgenelist_number\tgenelist_name\tgene\tchr\ttxStart\ttxEnd\tstrand\thopach_rank\n" > output_row_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s:%s:%s:%s:%s:%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$5,$6,$7,$8,$9,$12,$1,$2,$5,$6,$7,$8,$9,$12)}}' output_rna-seq-forheatmap.tsv > output_row_metadata.tmp
# ensure data rows are unique
sort output_row_metadata.tmp | uniq >> output_row_metadata.tsv

# col metadata file
printf "cid\texperiment_type\tsample_name\tdata_type\n" > output_col_metadata.tsv
awk -F'\t' '{if(NR!=1){printf("%s:%s:%s\t%s\t%s\t%s\n",$3,$4,"expression",$3,$4,"expression")}}' output_rna-seq-forheatmap.tsv > output_col_metadata.tmp
# get unique rows of col metadata (many are repeated per gene, since that's part of row metadata and is omitted here)
sort output_col_metadata.tmp | uniq >> output_col_metadata.tsv

# run r script to generate gct data file and morpheus heatmap
Rscript /usr/local/bin/run_genelists_heatmap_rnaseq-scaled.R output_row_metadata.tsv output_col_metadata.tsv output_counts.tsv ./
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
mv heatmap.html heatmap_scaled95.html





printf "\n\nWorkflow script run_genelists_rnaseq.sh complete!\n"