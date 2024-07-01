#!/bin/bash

#   Shell wrapper for 16S sample quantitation
#
##########################################################################################
#
# v1.0.0
# - 
#
##########################################################################################
printf "$(date)\nLog file for run_qiime2_sample_pe.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_qiime2_sample_pe.sh\`:
Shell wrapper for import and quantitation of one paired-end 16S sequencing sample using DADA2. Alpha rarefaction and taxonomic classification plots are also output.
Taxonomy classification is performed using a Naive Bayes classifier trained on the Greengenes2 database "gg_2022_10_backbone_full_length.nb.qza".
Generally, this workflow follows the "moving-pictures" turorial: https://docs.qiime2.org/2023.5/tutorials/moving-pictures/

    Output files:
    - overview.md, list of inputs
    - demux.qzv, summary visualizations of imported data
    - alpha-rarefaction.qzv, plot of OTU rarefaction
    - taxa-bar-plots.qzv, relative frequency of taxomonies barplot


PARAMS:
    SECTION 1: general
    -h  help   show this message
    -t  INT    number of threads
    -s  STR    sample name
    -a  FILE   path to read1 fastq file
    -b  FILE   path to read2 fastq file
    -j  J      trims the first J bases from the 5' end of each forward read
    -k  K      trims the first K bases from the 5' end of each reverse read
    -m  M      clips the forward read starting M bases from the 5' end (before trimming)
    -n  N      clips the reverse read starting N bases from the 5' end (before trimming)

____________________________________________________________________________________________________
References:
    Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
    
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:s:a:b:j:k:m:n:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		s) SAMPLENAME=$OPTARG ;;
		a) READ1=$OPTARG ;;
		b) READ2=$OPTARG ;;
    j) trimLeftF=$OPTARG ;;
    k) trimLeftR=$OPTARG ;;
    m) truncLenF=$OPTARG ;;
    n) truncLenR=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$SAMPLENAME" == "" ]]; then echo "error: required param missing (-s)"; exit; fi
if [[ "$READ1" == "" ]]; then echo "error: required param missing (-a)"; exit; fi
if [[ ! -f "$READ1" ]]; then echo "error: file does not exist (-a)"; exit; fi
if [[ "$READ2" == "" ]]; then echo "error: required param missing (-b)"; exit; fi
if [[ ! -f "$READ2" ]]; then echo "error: file does not exist (-b)"; exit; fi
if [[ "$trimLeftF" == "" ]]; then trimLeftF=0; fi
if [[ "$trimLeftR" == "" ]]; then trimLeftR=0; fi
if [[ "$truncLenF" == "" ]]; then truncLenF=120; fi
if [[ "$truncLenR" == "" ]]; then truncLenR=120; fi

# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-p, \$SAMPLENAME, $SAMPLENAME\n"
printf "\t-a, \$READ1, $READ1\n"
printf "\t-b, \$READ2, $READ2\n"
printf "\t-j, \$trimLeftF, $trimLeftF\n"
printf "\t-k, \$trimLeftR, $trimLeftR\n"
printf "\t-m, \$truncLenF, $truncLenF\n"
printf "\t-n, \$truncLenR, $truncLenR\n"



#	MAIN
#===============================================================================
printf "\n\nStep 1 - Create sample manifest (for single sample)\n"
printf "%s\t%s\t%s\n" "sample-id" "forward-absolute-filepath" "reverse-absolute-filepath" > pe-manifest.tsv
printf "%s\t" "$SAMPLENAME" >> pe-manifest.tsv
printf "%s\t" "$READ1" >> pe-manifest.tsv
printf "%s" "$READ2" >> pe-manifest.tsv

printf "\n\nStep 2 - Import the data, make a qiime artifact (demux.qza), and summary visualization\n"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-manifest.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
# summary visualization of raw read stats
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

printf "\n\nStep 3 - Denoise with DADA2"
# detect and correct (where possible) Illumina amplicon sequence data. This quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-max-ee-f 5 \
  --p-max-ee-r 5 \
  --p-trim-left-f $trimLeftF \
  --p-trim-left-r $trimLeftR \
  --p-trunc-len-f $truncLenF \
  --p-trunc-len-r $truncLenR \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza \
  --p-n-threads $THREADS
# summarize denoising
qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

printf "\n\nStep 4 - Generate a tree for phylogenetic diversity analyses"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

printf "\n\nStep 5 - Alpha rarefaction plotting"
# spoof metadata file
printf "%s\t%s\n" "sample-id" "condition" > sample-metadata.tsv
printf "%s\t%s\n" "$SAMPLENAME" "single-sample-processing" >> sample-metadata.tsv
# run rarefaction
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 1000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

printf "\n\nStep 6 - Taxonomic analysis"
qiime feature-classifier classify-sklearn \
  --i-classifier /dockerdata/gg_2022_10_backbone_full_length.nb.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
# summarize taxonomic analysis
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
# visualize taxonomic analysis
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv






#	OUTPUTS
#===============================================================================
printf "\tformatting...\n"

printf "## INPUTS\n" > overview.md
printf "\n" >> overview.md
printf "#### SECTION 1: general\n" >> overview.md
printf "-" >> overview.md
printf " \$THREADS, $THREADS\n" >> overview.md
printf "-" >> overview.md
printf " \$SAMPLENAME, $SAMPLENAME\n" >> overview.md
printf "-" >> overview.md
printf " \$READ1, $READ1\n" >> overview.md
printf "-" >> overview.md
printf " \$READ2, $READ2\n" >> overview.md
printf "-" >> overview.md
printf " \$trimLeftF, $trimLeftF\n" >> overview.md
printf "-" >> overview.md
printf " \$trimLeftR, $trimLeftR\n" >> overview.md
printf "-" >> overview.md
printf " \$truncLenF, $truncLenF\n" >> overview.md
printf "-" >> overview.md
printf " \$truncLenR, $truncLenR\n" >> overview.md
printf "\n" >> overview.md

printf "\n\nWorkflow script run_qiime2_sample_pe.sh complete!\n"