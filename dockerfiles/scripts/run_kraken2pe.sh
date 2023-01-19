#!/bin/bash

#   NOT CURRENTLY IN USE
#   Shell wrapper for run_kraken2pe.sh
#
##########################################################################################
#
# v0.0.1
# - format input lists for rnbeads in Rscript
# - run rnbeads Rscript
##########################################################################################
printf "$(date)\nLog file for run_rnbeads_diff.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_rnbeads_diff.sh\`:
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

 -a  STRING     name of condition1


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
while getopts "ht:d:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
        d) DATABASE=$OPTARG ;;

		a) CONDITION1_NAME=$OPTARG ;;
		b) CONDITION2_NAME=$OPTARG ;;

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
printf "\tCONDITION1_NAME - $CONDITION1_NAME\n"

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


#	OUTPUTS
#===============================================================================
# 

#   format kraken report for table





# clean up
rm -rf "$workdir"