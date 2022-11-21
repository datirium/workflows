#!/bin/bash

#   Shell wrapper for run_rnbeads_diff.R
#
##########################################################################################
#
# v0.0.1
# - format input lists for rnbeads in Rscript
#
##########################################################################################
printf "$(date)\nLog file for run_rnbeads_diff.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_rnbeads_diff.sh\`:
	Wrapper for RnBeads differential methylation pipeline.
    Output reports directory in container at '/tmp/reports/', includes:
    reports/
    ├── data_import.html
    ├── differential_methylation.html
    ├── preprocessing.html
    ├── quality_control.html
    ├── tracks_and_tables.html

OPTIONS:
 -h  help	show this message
 -g  STRING   Sample genome, available options: hg19, hg38, mm9, mm10, rn5
 -t  INT	number of threads
 -a  STRING     name of condition1
 -b  STRING     name of condition2
 -c  LIST	comma separated list of absolute filepaths to all condition1 bed files (BismarkCov format)
 -d  LIST	comma separated list of absolute filepaths to all condition2 bed files (BismarkCov format)
OPTIONAL:
 -o  DIR        absolute path to output directory for reports, default '/tmp/[reports]'

BismarkCov formatted bed:
    https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
    The genome-wide cytosine report (optional) is tab-delimited in the following format (1-based coords):
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

____________________________________________________________________________________________________
References:
	https://rnbeads.org/materials/example_3/differential_methylation.html
        Makambi, K. (2003) Weighted inverse chi-square method for correlated significance tests. Journal of Applied Statistics, 30(2), 225234
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4216143/
        Assenov Y, Müller F, Lutsik P, Walter J, Lengauer T, Bock C. Comprehensive analysis of DNA methylation data with RnBeads. Nat Methods. 2014 Nov;11(11):1138-1140. doi: 10.1038/nmeth.3115. Epub 2014 Sep 28. PMID: 25262207; PMCID: PMC4216143.
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:g:a:b:c:d:o:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
        g) GENOME=$OPTARG ;;
		a) CONDITION1_NAME=$OPTARG ;;
		b) CONDITION2_NAME=$OPTARG ;;
		c) CONDITION1_BED_FILEPATHS=$OPTARG ;;
		d) CONDITION2_BED_FILEPATHS=$OPTARG ;;
        o) OUTDIR=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# defaults
workdir=$PWD/tmp
mkdir -p "$workdir"
# if -o param not specified, use current working dir
if [[ ! -v "$OUTDIR" ]]; then OUTDIR=$PWD; fi
#in dockercontainer, run from:
printf "List of defaults:\n"
printf "\tPWD - $PWD\n"
printf "\tworkdir - $workdir\n\n"
printf "List of inputs:\n"
printf "\tTHREADS - $THREADS\n"
printf "\tGENOME - $GENOME\n"
printf "\tCONDITION1_NAME - $CONDITION1_NAME\n"
printf "\tCONDITION2_NAME - $CONDITION2_NAME\n"
printf "\tCONDITION1_BED_FILEPATHS - $CONDITION1_BED_FILEPATHS\n"
printf "\tCONDITION2_BED_FILEPATHS - $CONDITION2_BED_FILEPATHS\n"
printf "\tOUTDIR - $OUTDIR\n"


#	MAIN
#===============================================================================
# copy all data into workdir, and use absolute path as data.source[1] input vector
#   also generate sample_annotation.csv
printf "Sample_ID,condition\n" > $workdir/sample_annotation.csv
echo "$CONDITION1_BED_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    cp "$filepath" "$workdir/"
    printf "%s,%s\n" $(basename "$filepath") "$CONDITION1_NAME";
done >> $workdir/sample_annotation.csv
echo "$CONDITION2_BED_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    cp "$filepath" "$workdir/"
    printf "%s,%s\n" $(basename "$filepath") "$CONDITION2_NAME";
done >> $workdir/sample_annotation.csv
# check sample sheet csv
cat $workdir/sample_annotation.csv

# run Rscript wrapper for rnbeads with formatted inputs
Rscript /usr/local/bin/run_rnbeads_diff.R $workdir/sample_annotation.csv "$OUTDIR" "$workdir" "$GENOME" "$THREADS"

#	OUTPUTS
#===============================================================================
# move sample sheet for output
cp $workdir/sample_annotation.csv ./
# format for Overview tab
head -1 sample_annotation.csv | awk -F',' '{printf("| %s\t| %s\t|\n",$1,$2)}' > sample_annotation.md
printf "| --------\t| --------\t|\n" >> sample_annotation.md
tail -n+2 sample_annotation.csv | awk -F',' '{printf("| %s\t| %s\t|\n",$1,$2)}' >> sample_annotation.md
# package report dir
tar -cf reports.tar ./reports
gzip reports.tar

# clean up
rm -rf "$workdir"