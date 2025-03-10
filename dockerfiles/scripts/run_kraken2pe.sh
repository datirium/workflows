#!/bin/bash
printf "k2-classify-pe.cwl\n$(date)\n"          # prints tool name/date to stdout (error_report.txt)
printf "k2-classify-pe.cwl\n$(date)\n" 1>&2     # prints tool name/date to stderr (error_msg.txt)

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_kraken2download.sh\`:

PARAMS:
 -h     help    show this message
 -d     STRING  kraken2 database 
 -a     FILE    paired-end fastq R1 file
 -b     FILE    paired-end fastq R2 file


Kraken-report format:
    1. Percentage of fragments covered by the clade rooted at this taxon
    2. Number of fragments covered by the clade rooted at this taxon
    3. Number of fragments assigned directly to this taxon
    4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass,
        (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a
        rank code that is formed by using the rank code of the closest ancestor rank with a
        number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon
        is between genus and species and the grandparent taxon is at the genus rank.
    5. NCBI taxonomic ID number
    6. Indented scientific name

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "hd:a:b:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
        d) DATABASE=$OPTARG ;;
        a) R1=$OPTARG ;;
        b) R2=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$DATABASE" == "" ]]; then
    echo "Missing required Database input param (-d). Exiting."
    exit
fi
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
printf "\trun classification for PE reads\n"
kraken2 --db $DATABASE --threads 20 --paired --classified-out k2_classified_reads#.fastq --unclassified-out k2_unclassified_reads#.fastq --output k2.output --report k2.report $R1 $R2 2> k2.stderr
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from kraken2 execution.\n\nPlease check error_msg.txt file for details."
    exit
fi
# check that k2 output is not empty
if [[ $(awk 'END{print(NR)}' k2.report) == "0" ]]; then
    echo "k2.report file is empty; exiting"
    exit 1
fi
printf "\tformatting outputs\n"
# make stderr output markdown compatible for overview tab view
head -1 k2.stderr > parsed.stderr
tail -n+2 k2.stderr | sed 's/^ *//' | awk '{printf(" - %s\n",$0)}' >> parsed.stderr
# format report into tsv for table tab view
printf "percent_classified\treads_assigned_at_and_below_taxid\treads_assigned_directly_to_taxid\ttaxonomic_rank\ttaxid\tname\n" > k2_report.tsv
sed 's/^ *//' k2.report | sed 's/\t  */\t/' >> k2_report.tsv
printf "\tgenerate krona plot\n"
python3 /usr/local/src/KrakenTools/kreport2krona.py -r k2.report -o k2.krona --no-intermediate-ranks
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from kreport2krona.py execution.\n\nPlease check error_msg.txt file for details."
    exit
fi
# check that kreport2krona didn't produce no or empty output
if [[ ! -f k2.krona ]]; then
    echo "k2.krona file was not produced; exiting"
    exit 1
fi
if [[ $(awk 'END{print(NR)}' k2.krona) == "0" ]]; then
    echo "k2.krona file is empty; exiting"
    exit 1
fi
perl /usr/local/src/Krona/KronaTools/scripts/ImportText.pl -o krona.html k2.krona
if [[ $? == 1 ]]; then
    printf "Non-zero exit status from ImportText.pl execution.\n\nPlease check error_msg.txt file for details."
    exit
fi
# check that transform to html didn't produce no or empty output
if [[ ! -f krona.html ]]; then
    echo "krona.html file was not produced; exiting"
    exit 1
fi
if [[ $(awk 'END{print(NR)}' krona.html) == "0" ]]; then
    echo "krona.html file is empty; exiting"
    exit 1
fi