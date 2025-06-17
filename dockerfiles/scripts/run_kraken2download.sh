#!/bin/bash
printf "k2-download-db.cwl\n$(date)\n"          # prints tool name/date to stdout (error_report.txt)
printf "k2-download-db.cwl\n$(date)\n" 1>&2     # prints tool name/date to stderr (error_msg.txt)

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_kraken2download.sh\`:

PARAMS:
 -h  help	show this message
 -d  STRING   kraken2 database enum (acceptable strings are: Viral, Standard, Standard-16, MinusB, PlusPFP-16, EuPathDB46, 16S_Greengenes, 16S_Silva_138)

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "hd:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
        d) DATABASE=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$DATABASE" == "" ]]; then
    echo "Missing required Database input param (-d). Exiting."
    exit
fi

#	MAIN
#===============================================================================
if [[ "$DATABASE" == "Viral" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240904.tar.gz"
db="k2_viral_20240904"
elif [[ "$DATABASE" == "Standard" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz"
db="k2_standard_20240904"
elif [[ "$DATABASE" == "Standard-16" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240904.tar.gz"
db="k2_standard_16gb_20240904"
elif [[ "$DATABASE" == "MinusB" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20240904.tar.gz"
db="k2_minusb_20240904"
elif [[ "$DATABASE" == "PlusPFP-16" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240904.tar.gz"
db="k2_pluspfp_16gb_20240904"
elif [[ "$DATABASE" == "EuPathDB46" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20230407.tar.gz"
db="k2_eupathdb48_20230107"
elif [[ "$DATABASE" == "16S_Greengenes" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz"
db="k2_16S_Greengenes_20200326"
elif [[ "$DATABASE" == "16S_Silva_138" ]]; then
url="https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz"
db="k2_16S_Silva_138_20200326"
fi
printf "Downloading Kraken2 $DATABASE database from: $url\n\n"
wget $url
if [[ $? == 1 ]]; then
    printf "Error from wget, possible URL (tried: '$url') issue.\n\nPlease check error_msg.txt file for details."
    exit
fi
mkdir ./k2db
tar -xf *.tar.gz -C ./k2db