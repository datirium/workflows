#!/bin/bash

#   NOT CURRENTLY IN USE (rolled up into the cwl tool)
#   Shell wrapper for run_kraken2pe.sh
#
##########################################################################################
printf "$(date)\nLog file for run_kraken2pe.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_kraken2pe.sh\`:
    Wrapper for kraken2 taxonomic sequence classification system and associated report and visualization (viz pending)

    Primary Output files:
     - kraken2.tsv
     - kraken2.report

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

    Databases (various combinations of RefSeq databases):
        Viral (0.5 GB)          https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz
            All refseq viral genomes
        MinusB (8.7 GB)         https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20221209.tar.gz 
            Standard minus bacteria (archaea, viral, plasmid, human1, UniVec_Core)
        PlusPFP-16 (15.0 GB)    https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz
            Standard (archaea, bacteria, viral, plasmid, human1, UniVec_Core) + (protozoa, fungi & plant) capped at 16 GB (shrunk via random kmer downselect)
        EuPathDB46 (34.1 GB)    https://genome-idx.s3.amazonaws.com/kraken/k2_eupathdb48_20201113.tar.gz
            Eukaryotic pathogen genomes with contaminants removed
            reference: https://veupathdb.org/veupathdb/app

PARAMS:
 -h  help	show this message
 -t  INT	number of threads
 -d  STRING   kraken2 database enum (Viral, MinusB, PlusPFP-16, EuPathDB46) - this will be determined by the cwl tool calling this script


____________________________________________________________________________________________________
References:
    Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:d:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
        d) DATABASE=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# defaults
workdir=$PWD/tmp
mkdir -p "$workdir"
# use current working dir as outdir
OUTDIR=$PWD
#in dockercontainer, run from:
printf "List of defaults:\n"
printf "\tPWD - $PWD\n"
printf "\tworkdir - $workdir\n\n"
printf "List of inputs:\n"
printf "\tTHREADS - $THREADS\n"
printf "\tGENOME - $GENOME\n"
printf "\tOUTDIR - $OUTDIR\n"


#	MAIN
#===============================================================================
# 

# Database check
if [[ "$DATABASE" == "Viral" ]]; then dbpath="/data/k2_viral_20221209"
elif [[ "$DATABASE" == "MinusB" ]]; then dbpath="/data/k2_minusb_20221209"
elif [[ "$DATABASE" == "PlusPFP-16" ]]; then dbpath="/data/k2_pluspfp_16gb_20221209"
elif [[ "$DATABASE" == "EuPathDB46" ]]; then dbpath="/data/k2_eupathdb48_20201113"; fi




#   PRINT TEST READS
# test for human (human mitochondria)
printf ">hgmito-9606-read1_1\nGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTC\n" > metatest_R1.fa
printf ">hgmito-9606-read1_2\nATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC\n" > metatest_R2.fa
# test for k2_pluspfp_16gb_20221209 (e coli)
printf ">ecoli-562-read1_1\nAGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA\n" >> metatest_R1.fa
printf ">ecoli-562-read1_2\nTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC\n" >> metatest_R2.fa
# test for k2_viral_20221209 (phiX bacteriophage)
printf ">phix-2886930-read1_1\nTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC\n" >> metatest_R1.fa
printf ">phix-2886930-read1_2\nTGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGC\n" >> metatest_R2.fa
R1="metatest_R1.fa"
R2="metatest_R2.fa"
dbpath="/data/k2_pluspfp_16gb_20221209"
dbpath="/data/k2_viral_20221209"
THREADS=4

# run classification for PE reads
#   get read file extension
read_extension=$(basename "$R1" | awk -F'.' '{print($NF)}')
kraken2 --db $dbpath --threads $THREADS --paired --classified-out classified_reads#.$read_extension --output k2.output --report k2.report $R1 $R2 2> k2.summary

# for SE read classification (will use in different script)
#kraken2 --db /data/k2db/k2_pluspfp_16gb_20221209/ --threads 2 $FASTQ --output kraken.tsv --report kraken.report


# clean up
rm -rf "$workdir"