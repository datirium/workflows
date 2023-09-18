#!/bin/bash

#   Shell wrapper for 16S sample quantitation
#
##########################################################################################
#
# v1.0.0
# - 
#
##########################################################################################
printf "$(date)\nLog file for run_qiime2_aggregate.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_qiime2_aggregate.sh\`:
Shell wrapper for import and quantitation of 2+ paired-end 16S sequencing samples that (each) have already been processed in the SciDAP platform (using the qiime2-sample-pe.cwl workflow).
Alpha rarefaction and taxonomic classification plots are also output for the aggregated samples.
Taxonomy classification is performed using a Naive Bayes classifier trained on the Greengenes2 database "gg_2022_10_backbone_full_length.nb.qza".
Generally, this workflow follows the "moving-pictures" turorial: https://docs.qiime2.org/2023.5/tutorials/moving-pictures/
If "Metadata header name for PCoA axis label" is provided, principle coordinates analysis (PCoA) will be performed using the unweighted unifrac and bray curtis methods. 3D plots are produced with PCo1, PCo2, and the provided axis label on the x, y, and z axes.
If the sampling depth and metadata header for differential analysis are provided, differential abundance analysis will be performed using Gneiss and ANCOM methods at the family, genus, and species taxonomic levels. A unsupervised hierarchical clustering heatmap (Gneiss) and volcano plot (ANCOM) are produced at the taxonomic level between the specified group.

    Primary output files:
    - overview.md, list of inputs
    - demux.qzv, summary visualizations of imported data
    - alpha-rarefaction.qzv, plot of OTU rarefaction
    - taxa-bar-plots.qzv, relative frequency of taxomonies barplot
    - table.qza, table containing how many sequences are associated with each sample and with each feature (OTU)
    Optional output files:
    - pcoa-unweighted-unifrac-emperor.qzv, PCoA using unweighted unifrac method
    - pcoa-bray-curtis-emperor.qzv, PCoA using bray curtis method
    - heatmap.qzv, output from gneiss differential abundance analysis using unsupervised correlation-clustering method (this will define the partitions of microbes that commonly co-occur with each other using Ward hierarchical clustering)
    - ancom-\$LEVEL.qzv, output from ANCOM differential abundance analysis at family, genus, and species taxonomic levels (includes volcano plot)

PARAMS:
    SECTION 1: general
    -h  help         show this message
    -t  INT          number of threads
    -r  FILE         sample metadata file (sample-id [col1] should be identical to sample names being aggregated)
    -a  FILE ARRAY   array (csv) of paths to read1 fastq files
    -b  FILE ARRAY   array (csv) of paths to read2 fastq files
    -j  INT ARRAY    trims the first J bases from the 5' end of each forward read (first int in csv array)
    -k  INT ARRAY    trims the first K bases from the 5' end of each reverse read (first int in csv array)
    -m  INT ARRAY    clips the forward read starting M bases from the 5' end (before trimming) (first int in csv array)
    -n  INT ARRAY    clips the reverse read starting N bases from the 5' end (before trimming) (first int in csv array)
    -c  STR          custom axis label
Must be identical to one of the headers of the sample-metadata file. The corresponding column may only contain INT data.
    -d  INT          rarefaction sampling depth (required for differential abundance execution)
This step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of INT. If the total count for any sample(s) are smaller than this value, those samples will be dropped from the diversity analysis. It's recommend making your choice by reviewing the rarefaction plot. Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.
    -g  STR          group or experimental condition column name from sample metadata file (required for differential abundance execution)
Must be identical to one of the headers of the sample-metadata file. The corresponding column should only have two groups/conditions.

NOTES:
  Example sample metadata file (-r):
  Spacing should be single-tab separated (below things are lined up for clarity)

sample-id       sample-name     condition       time
SRR25508255     sample1 c1      1
SRR25508256     sample2 c1      1
SRR25508257     sample3 c2      2
SRR25508258     sample4 c2      2
 
____________________________________________________________________________________________________
References:
    Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9

EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:r:a:b:j:k:m:n:c:d:g:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		r) METADATAFILE=$OPTARG ;;
		a) READ1ARRAY=$OPTARG ;;
		b) READ2ARRAY=$OPTARG ;;
		j) trimLeftFarr=$OPTARG ;;
		k) trimLeftRarr=$OPTARG ;;
		m) truncLenFarr=$OPTARG ;;
		n) truncLenRarr=$OPTARG ;;
		c) CUSTOMLABEL=$OPTARG ;;
		d) SAMPLINGDEPTH=$OPTARG ;;
		g) GROUP=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$METADATAFILE" == "" ]]; then echo "error: required param missing (-r)"; exit; fi
if [[ ! -f "$METADATAFILE" ]]; then echo "error: file does not exist (-r)"; exit; fi
if [[ "$READ1ARRAY" == "" ]]; then echo "error: required param missing (-a)"; exit; fi
if [[ "$READ2ARRAY" == "" ]]; then echo "error: required param missing (-b)"; exit; fi
if [[ "$trimLeftFarr" == "" ]]; then echo "error: required param missing (-j)"; exit; fi
if [[ "$trimLeftRarr" == "" ]]; then echo "error: required param missing (-k)"; exit;  fi
if [[ "$truncLenFarr" == "" ]]; then echo "error: required param missing (-m)"; exit;  fi
if [[ "$truncLenRarr" == "" ]]; then echo "error: required param missing (-n)"; exit;  fi

# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-r, \$METADATAFILE, $METADATAFILE\n"
printf "\t-a, \$READ1ARRAY, $READ1ARRAY\n"
printf "\t-b, \$READ2ARRAY, $READ2ARRAY\n"
printf "\t-j, \$trimLeftF, $trimLeftFarr\n"
printf "\t-k, \$trimLeftR, $trimLeftRarr\n"
printf "\t-m, \$truncLenF, $truncLenFarr\n"
printf "\t-n, \$truncLenR, $truncLenRarr\n"
printf "\t-c, \$CUSTOMLABEL, $CUSTOMLABEL\n"
printf "\t-d, \$SAMPLINGDEPTH, $SAMPLINGDEPTH\n"
printf "\t-g, \$GROUP, $GROUP\n"


#	MAIN
#===============================================================================
printf "\n\nStep 1 - Create sample manifest\n"
# cannot import fastq in bz2 compressed format, need to change to gz
printf "$READ1ARRAY\n" | sed 's/,/\n/g' | while read f; do
  bn=$(printf "$f" | sed -e 's/\/.*\///' -e 's/\..*//') # get basename for each file
  bunzip2 --stdout $f | gzip > $bn.fastq.gz             # decompress from bz2, recompress with gz
  printf "$PWD/$bn.fastq.gz," >> read1path.array.csv    # store absolute file paths in csv file
done
printf "$READ2ARRAY\n" | sed 's/,/\n/g' | while read f; do
  bn=$(printf "$f" | sed -e 's/\/.*\///' -e 's/\..*//') # get basename for each file
  bunzip2 --stdout $f | gzip > $bn.fastq.gz             # decompress from bz2, recompress with gz
  printf "$PWD/$bn.fastq.gz," >> read2path.array.csv    # store absolute file paths in csv file
done
# make manifest using gz files
printf "%s\t" "sample-id" > pe-manifest.transpose
cut -f1 "$METADATAFILE" | tail -n+2 | awk '{printf("%s\t",$0)}END{printf("\n")}' >> pe-manifest.transpose
printf "%s\t" "forward-absolute-filepath" >> pe-manifest.transpose
sed 's/,$/\n/' read1path.array.csv | sed 's/,/\t/g' >> pe-manifest.transpose
printf "%s\t" "reverse-absolute-filepath" >> pe-manifest.transpose
sed 's/,$/\n/' read2path.array.csv | sed 's/,/\t/g' >> pe-manifest.transpose
# transpose into expected manifest format
transpose_tsv.awk pe-manifest.transpose > pe-manifest.tsv

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

printf "\n\nStep 3 - Denoise with DADA2\n"
# get the first int in each trim/trunc array to use for parameters
trimLeftF=$(printf "$trimLeftFarr\n" | sed 's/,.*//')
trimLeftR=$(printf "$trimLeftRarr\n" | sed 's/,.*//')
truncLenF=$(printf "$truncLenFarr\n" | sed 's/,.*//')
truncLenR=$(printf "$truncLenRarr\n" | sed 's/,.*//')
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

printf "\n\nStep 4 - Generate a tree for phylogenetic diversity analyses\n"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

printf "\n\nStep 5 - Alpha rarefaction plotting\n"
# rename input metadata file
cp $METADATAFILE sample-metadata.tsv
# run rarefaction
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 10000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

printf "\n\nStep 6 - Taxonomic analysis\n"
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



# if a sampling depth is provided, run these
printf "\n\nStep 7 - Checking input params for execution of diversity, pcoa, and differential abundance analysis\n"
if [[ "$SAMPLINGDEPTH" != "" && "$GROUP" != "" ]]; then
  printf "\n\nStep 7 - Input params SAMPLINGDEPTH (-d) and GROUP (-g) are not empty. Proceeding...\n"
  printf "\n\nStep 7a - Diversity analysis\n"
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table.qza \
    --p-sampling-depth $SAMPLINGDEPTH \
    --m-metadata-file sample-metadata.tsv \
    --output-dir core-metrics-results

  if [[ "$CUSTOMLABEL" != "" ]]; then
    printf "\n\nStep 7b - PCoA analysis\n"
    # PCoA of unweighted unifrac
    qiime emperor plot \
      --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
      --m-metadata-file sample-metadata.tsv \
      --p-custom-axes $CUSTOMLABEL \
      --o-visualization core-metrics-results/pcoa-unweighted-unifrac-emperor.qzv
    # PCoA of bray curtis
    qiime emperor plot \
      --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
      --m-metadata-file sample-metadata.tsv \
      --p-custom-axes $CUSTOMLABEL \
      --o-visualization core-metrics-results/pcoa-bray-curtis-emperor.qzv
  fi

  printf "\n\nStep 7c1 - Differential abundance analysis with gneiss\n"
  qiime gneiss correlation-clustering \
    --i-table core-metrics-results/rarefied_table.qza \
    --o-clustering hierarchy.qza
  # generate heatmap to see which groups of taxa they represent
  qiime gneiss dendrogram-heatmap \
    --i-table core-metrics-results/rarefied_table.qza \
    --i-tree hierarchy.qza \
    --m-metadata-file sample-metadata.tsv \
    --m-metadata-column $GROUP \
    --p-color-map seismic \
    --o-visualization heatmap.qzv

  printf "\n\nStep 7c2 - Differential abundance analysis with ANCOM at family, genus, and species levels between $GROUP groups/conditions\n"
  for LEVEL in "Family" "Genus" "Species"; do
    if [[ "$LEVEL" == "Kingdom" ]]; then LEVELNUMBER=1; fi
    if [[ "$LEVEL" == "Phylum" ]]; then LEVELNUMBER=2; fi
    if [[ "$LEVEL" == "Class" ]]; then LEVELNUMBER=3; fi
    if [[ "$LEVEL" == "Order" ]]; then LEVELNUMBER=4; fi
    if [[ "$LEVEL" == "Family" ]]; then LEVELNUMBER=5; fi
    if [[ "$LEVEL" == "Genus" ]]; then LEVELNUMBER=6; fi
    if [[ "$LEVEL" == "Species" ]]; then LEVELNUMBER=7; fi
    # collapse table to user-specified taxa level
    qiime taxa collapse \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --p-level $LEVELNUMBER \
      --o-collapsed-table collapsed-table-$LEVEL.qza
    # add pseudocounts of 1 where value == 0
    qiime composition add-pseudocount \
      --i-table collapsed-table-$LEVEL.qza \
      --o-composition-table comp-collapsed-table-$LEVEL.qza
    # run ancom
    qiime composition ancom \
      --i-table comp-collapsed-table-$LEVEL.qza \
      --m-metadata-file sample-metadata.tsv \
      --m-metadata-column $GROUP \
      --o-visualization ancom-$LEVEL.qzv
  done

fi


#	OUTPUTS
#===============================================================================
printf "\n\n\tformatting overview.md file...\n"

printf "## INPUTS\n" > overview.md
printf "\n" >> overview.md

printf "#### SECTION 1: general\n" >> overview.md

printf "-" >> overview.md
printf " \$THREADS, $THREADS\n" >> overview.md
printf "-" >> overview.md
printf " \$METADATAFILE, $METADATAFILE\n" >> overview.md
printf "-" >> overview.md
printf " \$READ1ARRAY, $READ1ARRAY\n" >> overview.md
printf "-" >> overview.md
printf " \$READ2ARRAY, $READ2ARRAY\n" >> overview.md
printf "-" >> overview.md
printf " \$trimLeftF, $trimLeftF\n" >> overview.md
printf "-" >> overview.md
printf " \$trimLeftR, $trimLeftR\n" >> overview.md
printf "-" >> overview.md
printf " \$truncLenF, $truncLenF\n" >> overview.md
printf "-" >> overview.md
printf " \$truncLenR, $truncLenR\n" >> overview.md
printf "-" >> overview.md
printf " \$CUSTOMLABEL, $CUSTOMLABEL\n" >> overview.md
printf "-" >> overview.md
printf " \$SAMPLINGDEPTH, $SAMPLINGDEPTH\n" >> overview.md
printf "-" >> overview.md
printf " \$GROUP, $GROUP\n" >> overview.md
printf "-" >> overview.md
printf " \$LEVEL, $LEVEL\n" >> overview.md

printf "\n" >> overview.md

printf "\n\nWorkflow script run_qiime2_aggregate.sh complete!\n"