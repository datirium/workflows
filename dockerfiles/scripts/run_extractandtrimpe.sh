#!/bin/bash
printf "extractandtrim-pe.cwl\n$(date)\n"          # prints tool name/date to stdout (error_report.txt)
printf "extractandtrim-pe.cwl\n$(date)\n" 1>&2     # prints tool name/date to stderr (error_msg.txt)

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_extractandtrimpe.sh\`:

PARAMS:
 -h     help    show this message
 -a     FILE    paired-end fastq R1 file
 -b     FILE    paired-end fastq R2 file

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ha:b:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
        a) R1=$OPTARG ;;
        b) R2=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$R1" == "" ]]; then
    echo "Missing required read1 input param (-a). Exiting."
    exit
fi
if [[ "$R2" == "" ]]; then
    echo "Missing required read2 input param (-b). Exiting."
    exit
fi

#	MAIN
#===============================================================================
printf "\tExtracting fastq files if compressed, and concatenate if more than 1 per pair\n"
function extract {
    FILE=$1
    COMBINED=$2
    T=`file -b "${FILE}" | awk '{print $1}'`
    case "${T}" in
        "bzip2"|"gzip"|"Zip")
        7z e -so "${FILE}" >> "${COMBINED}"
        ;;
        "ASCII")
        cat "${FILE}" >> "${COMBINED}" || true
        ;;
        *)
        echo "Error: file type '${T}' unknown" > error_report.txt
        rm -f "${COMBINED}"
        exit 1
    esac
}
# run extracts for r1 and r2 in parallel
COMBINED="extracted_combined_R1.fastq"
for FILE in "$R1"; do
    echo "Extracting:" $FILE;
    extract "${FILE}" "${COMBINED}"
done
COMBINED="extracted_combined_R2.fastq"
for FILE in "$R2"; do
    echo "Extracting:" $FILE;
    extract "${FILE}" "${COMBINED}"
done
printf "\n\n"



printf "\tTrim the extracted and concatentated paired-end reads with Trimgalore\n"
# expected outputs are:
#   extracted_combined_R1_val_1.fq
#   extracted_combined_R2_val_2.fq
#   extracted_combined_R1.fastq_trimming_report.txt
#   extracted_combined_R2.fastq_trimming_report.txt
trim_galore "extracted_combined_R1.fastq" "extracted_combined_R2.fastq" --dont_gzip --length 30 --paired --cores 8
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from trim_galore execution.\n\nPlease check error_msg.txt file for details."
    exit
fi
if [[ $(grep "Premature end of file encountered" error_msg.txt) != "" ]]; then
    echo "Premature end of file encountered, fastq is likely truncated or misformatted."
fi



printf "\tRunning fastx_quality_stats from Fastx Toolkit on trimmed reads\n"
fastx_quality_stats -i extracted_combined_R1_val_1.fq -o fastx_quality_stats_r1.tsv
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from fastx_quality_stats execution on trimmed r1 file.\n\nPlease check error_msg.txt file for details."
    exit
fi
fastx_quality_stats -i extracted_combined_R2_val_2.fq -o fastx_quality_stats_r2.tsv
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from fastx_quality_stats execution on trimmed r2 file.\n\nPlease check error_msg.txt file for details."
    exit
fi









