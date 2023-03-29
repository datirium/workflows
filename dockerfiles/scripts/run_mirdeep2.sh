#!/bin/bash

#   Shell wrapper for miRNA id and discovery
#
##########################################################################################
#
# v1.0.0
# - tbd
#
##########################################################################################
printf "$(date)\nLog file for run_mirdeep2.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_mirdeep2.sh\`:
Shell wrapper for miRDeep2 tool, that discovers known or novel miRNAs from deep sequencing data.
To be used in cwl mirdeep2.cwl tool and workflow.

Primary Output files:
 - mirs_known.tsv, known mature miRNAs detected by mirdeep2, "Known miRNAs" tab
 - mirs_novel.tsv, known novel miRNAs detected by mirdeep2, "Novel miRNAs" tab
Secondary Output files:
 - mirs_known_exocarta_deepmirs.tsv, list of detected miRNA also in ExoCarta's exosome database, "Detected Exosome miRNAs" tab
 - mirs_known_gene_targets.tsv, pre-computed gene targets of known mature mirs, downloadable
 - known_mirs_mature.fa, known mature mir sequences, downloadable
 - known_mirs_precursor.fa, known precursor mir sequences, downloadable
 - novel_mirs_mature.fa, novel precursor mir sequences, downloadable
 - novel_mirs_precursor.fa, novel precursor mir sequences, downloadable
Reports:
 - overview.md (input list, alignment & mir metrics), "Overview" tab
 - mirdeep2_result.html, summary of mirdeep2 results, "miRDeep2 Results" tab

PARAMS:
    SECTION 1: general
    -h	help	show this message
    -t  INT		number of threads
    -b  DIR		path to bowtie2 indices folder of genome reference
    -r  FILE	path to genome reference FASTA file
    -g  STRING	genome short name (hg19, hg38, mm10, rn7, dm3)
    -a  STRING	sequencing adapter for clipping reads (default: TCGTAT)
    -f  FILE	path to sequence read FASTQ file

GENOME:
For filtering mirbase by organism.
	genome	#organism   #division   #name   				#tree   																		#NCBI-taxid
	hg19	hsa     	HSA     	Homo sapiens    		Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606
	hg38	hsa     	HSA     	Homo sapiens    		Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Primates;Hominidae; 9606
	mm10	mmu     	MMU     	Mus musculus    		Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   		10090
	rn7		rno     	RNO     	Rattus norvegicus   	Metazoa;Bilateria;Deuterostoma;Chordata;Vertebrata;Mammalia;Rodentia;   		10116
	dm3		dme     	DME     	Drosophila melanogaster Metazoa;Bilateria;Ecdysozoa;Arthropoda;Hexapoda;        						7227

NOTES:
For the identification of novel miRNA candidates, the following may be used as a filtering guideline:
1. miRDeep score > 4 (but also some authors use 1 sometimes)
2. not present a match with rfam
3. should present a significant RNAfold ("yes")
4. a number of mature reads > 10
5. (optional) novel mir must be expressed in multiple samples

____________________________________________________________________________________________________
References:
 - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245920
 - https://github.com/rajewsky-lab/mirdeep2
 - https://biocontainers.pro/tools/mirdeep2
 - https://www.mirbase.org/
 - http://exocarta.org/index.html
 - https://www.targetscan.org/vert_80/
    
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:b:r:g:a:f:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		b) BT2DIR=$OPTARG ;;
		r) REFGENOME=$OPTARG ;;
		g) GENOME=$OPTARG ;;
		a) ADAPTER=$OPTARG ;;
		f) FASTQ=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$BT2DIR" == "" ]]; then echo "error: required param missing (-b)"; exit; fi
if [[ ! -d "$BT2DIR" ]]; then echo "error: dir does not exist (-b)"; exit; fi
if [[ "$REFGENOME" == "" ]]; then echo "error: required param missing (-r)"; exit; fi
if [[ ! -f "$REFGENOME" ]]; then echo "error: dir does not exist (-r)"; exit; fi
if [[ "$GENOME" == "" ]]; then echo "error: required param missing (-g)"; exit; fi
if [[ "$ADAPTER" == "" ]]; then echo "error: required param missing (-a)"; exit; fi
if [[ ! -f "$FASTQ" ]]; then echo "error: file does not exist (-f)"; exit; fi


# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-b, \$BT2DIR, $BT2DIR\n"
printf "\t-r, \$REFGENOME, $REFGENOME\n"
printf "\t-g, \$GENOME, $GENOME\n"
printf "\t-a, \$ADAPTER, $ADAPTER\n"
printf "\t-f, \$FASTQ, $FASTQ\n"



#	MAIN
#===============================================================================
printf "\n\nStep 1 - the Mapper module (processes Illumina output and maps it to the reference genome)\n"
#   ‑c	input file is FASTA format
#   ‑e	input file is FASTQ format
#   ‑h	parse to FASTA format                           <- required for miRDeep2.pl
#   ‑m	collapse reads                                  <- required for miRDeep2.pl
#   ‑p  <genome>	map to genome (must be indexed by bowtie-build). The genome string must be the prefix of the bowtie index. For instance, if the first indexed file is called h_sapiens_37_asm.1.ebwt then the prefix is h_sapiens_37_asm.
#   ‑s  file	print processed reads to this file      <- use this output as input to miRDeep2.pl
#   -t  file    print read mappings to this file
# get bt2 index basename for input
bt2name=$(basename $BT2DIR)
bt2index=$BT2DIR/$bt2name
# mapper
# 	note: converting fastq to fasta separately takes longer than simply using the FASTQ with the -e and -h params
time mapper.pl $FASTQ -e -h -i -j -k $ADAPTER -l 18 -m -p $BT2DIR/$bt2name -s collapsed.fasta -t mir.arf -o $THREADS


printf "\n\nStep 2 - the miRDeep2 module (identifies known and novel miRNAs in high-throughput sequencing data)\n"
# miRDeep2.pl reads genome mappings miRNAs_ref/none miRNAs_other/none precursors/none [options] 2>report.log
#   keep randfold, takes a while (~20 mins) but gives probabilities for TP mir detection
# grep out mir seqs from mirbase fasta:
#   - using genome 'short name' (miRNAs_ref)
#	- genome short name (hg19, hg38, mm10, rn7, dm3)
if [[ "$GENOME" == "hg19" ]]; then organism="hsa"; fi
if [[ "$GENOME" == "hg38" ]]; then organism="hsa"; fi
if [[ "$GENOME" == "mm10" ]]; then organism="mmu"; fi
if [[ "$GENOME" == "rn7" ]]; then organism="rno"; fi
if [[ "$GENOME" == "dm3" ]]; then organism="dme"; fi
mature="/dockerdata/refs/mature.fa"
precursors="/dockerdata/refs/hairpin.fa"
extract_miRNAs.pl $mature $organism > mirbase_mature.fasta
extract_miRNAs.pl $precursors $organism > mirbase_precursors.fasta
time miRDeep2.pl collapsed.fasta $REFGENOME mir.arf mirbase_mature.fasta none mirbase_precursors.fasta
# make mir results standardized output filenames for cwl
cp result_*.html mirdeep2_result.html
mv mirna_results_* mirna_results
cp mirna_results/known_mature_*.fa known_mirs_mature.fa
cp mirna_results/known_pres_*.fa known_mirs_precursor.fa
cp mirna_results/novel_mature_*.fa novel_mirs_mature.fa
cp mirna_results/novel_pres_*.fa novel_mirs_precursor.fa


#	Exocarta and TargetScan setup
if [[ $organism == "hsa" ]]; then exo_org="Homo sapiens"; taxid="9606"; fi
if [[ $organism == "mmu" ]]; then exo_org="Mus musculus"; taxid="10090"; fi
if [[ $organism == "dm3" ]]; then taxid="7227"; fi
# mirdeep2 list of detected mirs for overlapping with exocarta and targetscan
grep -A1000000 "^mature miRBase miRNAs detected by miRDeep2" result_*.csv | grep -B1000000 "^#miRBase miRNAs not detected by miRDeep2" | tail -n+2 | head -n-2 | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$10,$14,$2,$4,$9,$6)}' | sed 's/ /_/g' > mirs_known.tsv
# make input file for DESeq
#		get total read count for pseudoRPKM values
mature_read_total=$(tail -n+2 mirs_known.tsv | awk -F'\t' '{total+=$7}END{print(total)}')
tail -n+2 mirs_known.tsv | awk -F'\t' -v total="$mature_read_total" 'BEGIN{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand", "TotalReads", "Rpkm")};{split($1,chr_pos,"_"); printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.0f\n", $1, $2, chr_pos[1], chr_pos[2], chr_pos[2]+length($3), ".", $7, ($7/total)*1000000)}' > deseq_input.tsv
# mirdeep2 list of novel mirs (POSSIBLE DOWNSTREAM ANALYSIS INPUT - mature sequence used in a sequence-based target prediction tool)
grep -A1000000 "^novel miRNAs predicted by miRDeep2" result_*.csv | grep -B1000000 "^mature miRBase miRNAs detected by miRDeep2" | tail -n+2 | head -n-4 | awk -F'\t' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$10,$14,$2,$4,$9,$6)}' | sed 's/ /_/g' > mirs_novel.tsv
# trim down mir name to the part that will match mirs in other lists
tail -n+2 mirs_known.tsv | cut -f2 | sed 's/'"$organism"'-//' | sort | uniq > mirs_known_names_for_overlap.tsv

printf "\n\nStep 3 - exocarta (find overlap between mirs_known.tsv and exosome associated miRNAs)\n"
# exosome mir summary (only for hg19, hg38, mm10)
#   wget http://exocarta.org/Archive/EXOCARTA_MIRNA_DETAILS_5.txt
exocarta="/dockerdata/refs/EXOCARTA_MIRNA_DETAILS_5.txt"
# exocarta only has mirs for human and mouse
if [[ "$exo_org" != "" ]]; then
	# get list of organism specific exosome mirs from exocarta
	grep "$exo_org" $exocarta > exo_$organism.tsv
	# print exo lines that match trimmed mir names, only output "Clear hit to Entrez gene ID"
	head -1 $exocarta > mirs_known_exocarta.tsv
	awk -F'\t' '{if(FNR==NR){exo[$3]=$0}else{print(exo[$0])}}' exo_$organism.tsv mirs_known_names_for_overlap.tsv | sort | uniq | grep "Clear hit to Entrez gene ID" >> mirs_known_exocarta.tsv
	# make file with mirdeep2 stats for detected exocarta mirs
	head -1 mirs_known.tsv > mirs_known_exocarta_deepmirs.tsv
	cut -f3 mirs_known_exocarta.tsv | while read mir; do grep "$organism-$mir" mirs_known.tsv; done | sort | uniq >> mirs_known_exocarta_deepmirs.tsv
else
	printf "\tOrganism $organism does not have miRNAs in Exocarta DB for miRNA found in exosomes, skipping step."
fi

printf "\n\nStep 4 - TargetScan (find overlap between mirs_known.tsv and predicted miRNA targets per organism)\n"
# targetscan only has mirs for human, mouse, and fly
if [[ "$taxid" != "" ]]; then
	# set organism-specific target scan file
	targetscan="/dockerdata/refs/predicted_targets_${organism}_${taxid}.txt"
	# print targetscan lines that match trimmed mir names
	head -1 $targetscan > mirs_known_gene_targets.tsv
	while read mir; do grep "^$mir" $targetscan; done < mirs_known_names_for_overlap.tsv >> mirs_known_gene_targets.tsv
else
	printf "\tOrganism $organism does not have miRNAs in TargetScan DB for miRNA target predictions, skipping step."
fi



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
printf "2. not present a match with rfam\n" >> overview.md
printf "3. should present a significant RNAfold (yes)\n" >> overview.md
printf "4. a number of mature reads (>10)\n" >> overview.md
printf "5. if applicable, novel mir must be in >1 sample\n" >> overview.md
printf "\n" >> overview.md

printf "## INPUTS\n" >> overview.md
printf "-" >> overview.md
printf " \$THREADS, $THREADS\n" >> overview.md
printf "-" >> overview.md
printf " \$GENOME, $GENOME\n" >> overview.md
printf "-" >> overview.md
printf " \$ADAPTER, $ADAPTER\n" >> overview.md
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