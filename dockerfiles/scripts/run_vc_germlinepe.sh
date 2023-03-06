#!/bin/bash

#   Shell wrapper for germline variant calling
#
##########################################################################################
#
# v0.0.1
# - 
##########################################################################################
printf "$(date)\nLog file for run_vc_germline.sh\n\n"


#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_vc_germline.sh\`:
Shell wrapper for the Broad Institute's best practices gatk4 germline variant calling pipeline.

    Primary Output files:
    - bqsr2_indels.vcf, filtered and recalibrated indels
    - bqsr2_snps.ann.vcf, filtered and recalibrated snps with effect annotations
    Secondary Output files:
    - raw_indels.vcf, first pass indel calls
    - raw_snps.vcf, first pass snp calls
    Reports:
    - overview.md (input list, alignment metrics, variant counts)
    - insert_size_histogram.pdf
    - recalibration_plots.pdf
    - snpEff_summary.html

PARAMS:
    SECTION 1: general
    -h  help   show this message
    -t  INT    number of threads
    -r  DIR    path to bwa indices folder of genome reference
    -p  INT    ploidy, number of copies per chromosome (default should be 2)
    -s  DIR    path to SNPEFF database (see SNPEFFDB section below)
    -a  FILE   path to read1 fastq file
    -b  FILE   path to read2 fastq file
    SECTION 2: SNP filters (see Step 6 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
    -c  FLOAT  QD
    -d  FLOAT  FS
    -e  FLOAT  MQ
    -f  FLOAT  SOR
    -g  FLOAT  MQRankSum
    -i  FLOAT  ReadPosRankSum
    SECTION 3: Indel filters (see Step 7 Notes: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
    -j  FLOAT  QD
    -k  FLOAT  FS
    -l  FLOAT  SOR

SNPEFFDB:
    # snpeff databases
    docker run --rm -ti gatk4-dev /bin/bash
    java -jar $SNPEFF_JAR databases
    # use first column as SNPEFF -v input (e.g. "hg38")
    hg38, Homo_sapiens (USCS), http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg38.zip
    mm10, Mus_musculus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_mm10.zip
    dm6.03, Drosophila_melanogaster, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_dm6.03.zip
    Rnor_6.0.86, Rattus_norvegicus, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_Rnor_6.0.86.zip
    R64-1-1.86, Saccharomyces_cerevisiae, http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_R64-1-1.86.zip


____________________________________________________________________________________________________
References:
        https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:r:p:s:a:b:c:d:e:f:g:i:j:k:l:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		r) GENOMEDIR=$OPTARG ;;
		p) PLOIDY=$OPTARG ;;
		s) SNPEFFDB=$OPTARG ;;
		a) READ1=$OPTARG ;;
		b) READ2=$OPTARG ;;
		c) snp_QD=$OPTARG ;;
		d) snp_FS=$OPTARG ;;
		e) snp_MQ=$OPTARG ;;
		f) snp_SOR=$OPTARG ;;
		g) snp_MQRankSum=$OPTARG ;;
		i) snp_ReadPosRankSum=$OPTARG ;;
		j) indel_QD=$OPTARG ;;
		k) indel_FS=$OPTARG ;;
		l) indel_SOR=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$GENOMEDIR" == "" ]]; then echo "error: required param missing (-r)"; exit; fi
if [[ ! -d "$GENOMEDIR" ]]; then echo "error: dir does not exist (-r)"; exit; fi
if [[ "$SNPEFFDB" == "" ]]; then echo "error: required param missing (-s)"; exit; fi
if [[ "$READ1" == "" ]]; then echo "error: required param missing (-a)"; exit; fi
if [[ ! -f "$READ1" ]]; then echo "error: file does not exist (-a)"; exit; fi
if [[ "$READ2" == "" ]]; then echo "error: required param missing (-b)"; exit; fi
if [[ ! -f "$READ2" ]]; then echo "error: file does not exist (-b)"; exit; fi

# defaults
printf "List of defaults:\n"
printf "\t\$PWD - $PWD\n"
printf "List of inputs:\n"
printf "  SECTION 1: general\n"
printf "\t-t, \$THREADS, $THREADS\n"
printf "\t-p, \$PLOIDY, $PLOIDY\n"
printf "\t-r, \$GENOMEDIR, $GENOMEDIR\n"
printf "\t-s, \$SNPEFFDB, $SNPEFFDB\n"
printf "\t-a, \$READ1, $READ1\n"
printf "\t-b, \$READ2, $READ2\n"
printf "  SECTION 2: snp filters\n"
printf "\t-c, \$snp_QD, $snp_QD\n"
printf "\t-d, \$snp_FS, $snp_FS\n"
printf "\t-e, \$snp_MQ, $snp_MQ\n"
printf "\t-f, \$snp_SOR, $snp_SOR\n"
printf "\t-g, \$snp_MQRankSum, $snp_MQRankSum\n"
printf "\t-i, \$snp_ReadPosRankSum, $snp_ReadPosRankSum\n"
printf "  SECTION 3: indel filters\n"
printf "\t-c, \$indel_QD, $indel_QD\n"
printf "\t-d, \$indel_FS, $indel_FS\n"
printf "\t-f, \$indel_SOR, $indel_SOR\n"


#	MAIN
#===============================================================================
#printf "\n\nStep 0 - Build bwa mem index for '$GENOME'\n"
#bwa index $GENOME/genome.fasta
#   COMMENTED Step 0 OUT, TAKEN CARE OF WITH 'bwa-index.cwl' workflow
GENOME=$GENOMEDIR/genome.fasta

printf "\n\nStep 1 - Alignment to Reference\n"
bwa mem \
    -t $THREADS \
    -K 100000000 \
    -Y \
    -R '@RG\tID:id\tLB:library\tPL:platform\tPM:machine\tSM:sample' \
    $GENOME \
    $READ1 \
    $READ2 \
    > aligned_reads.sam

printf "\n\nStep 2 - Sort + Mark Duplicates\n"
samtools sort -@$THREADS -o sorted_reads.bam aligned_reads.sam
gatk MarkDuplicates \
    -I sorted_reads.bam \
    -M dedup_metrics.txt \
    -O sorted_dedup_reads.bam
samtools index -@$THREADS sorted_dedup_reads.bam
# generate coverage metrics file (percent coverage, depth mean, depth stdev)
/usr/local/bin/samtools-1.17/bin/samtools depth -@$THREADS -a sorted_dedup_reads.bam | cut -f3 | awk '{x+=$0;y+=$0^2;if($0==0){z++}}END{printf("%.2f\t%.0f\t%.0f\n",100*((NR-z)/NR),x/NR,sqrt(y/NR-(x/NR)^2))}' > coverage_metrics.tsv

printf "\n\nStep 3 - Collect Alignment & Insert Size Metrics\n"
java -jar $PICARD_JAR \
    CollectAlignmentSummaryMetrics \
    R=$GENOME \
    I=sorted_dedup_reads.bam \
    O=alignment_metrics.txt
java -jar $PICARD_JAR \
    CollectInsertSizeMetrics \
    INPUT=sorted_dedup_reads.bam \
    OUTPUT=insert_metrics.txt \
    HISTOGRAM_FILE=insert_size_histogram.pdf    # VISUALIZATION OPTION!
#samtools depth -a sorted_dedup_reads.bam > depth_out.txt   # not required downstream... (also, could take out positions with depth=0?)

printf "\n\nStep 4 - Call Variants\n"
gatk HaplotypeCaller \
    -R $GENOME \
    -ploidy $PLOIDY \
    -I sorted_dedup_reads.bam \
    -O raw_variants.vcf

printf "\n\nStep 5 - Extract SNPs & Indels\n"
gatk SelectVariants \
	-R $GENOME \
	-V raw_variants.vcf \
	--select-type-to-include SNP \
	--output raw_snps.vcf
gatk SelectVariants \
	-R $GENOME \
	-V raw_variants.vcf \
	--select-type-to-include INDEL \
	--output raw_indels.vcf

printf "\n\nStep 6 - Filter SNPs\n"
#   Note: SNPs which are ‘filtered out’ at this step will remain in the
#       filtered_snps.vcf file, however they will be marked as ‘_filter’,
#       while SNPs which passed the filter will be marked as ‘PASS’. We need to
#       extract and provide only the passing SNPs to the BQSR tool, we do this
#       in the next step (step 9). 
gatk VariantFiltration \
	-R $GENOME \
	-V raw_snps.vcf \
	--output filtered_snps.vcf \
	--filter-name "QD_filter" -filter "QD < $snp_QD" \
	--filter-name "FS_filter" -filter "FS > $snp_FS" \
	--filter-name "MQ_filter" -filter "MQ < $snp_MQ" \
	--filter-name "SOR_filter" -filter "SOR > $snp_SOR" \
	--filter-name "MQRankSum_filter" -filter "MQRankSum < $snp_MQRankSum" \
	--filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < $snp_ReadPosRankSum"

printf "\n\nStep 7 - Filter Indels\n"
gatk VariantFiltration \
	-R $GENOME \
	-V raw_indels.vcf \
	--output filtered_indels.vcf \
	--filter-name "QD_filter" -filter "QD < $indel_QD" \
	--filter-name "FS_filter" -filter "FS > $indel_FS" \
	--filter-name "SOR_filter" -filter "SOR > $indel_SOR"

printf "\n\nStep 8 - Exclude Filtered Variants\n"
gatk SelectVariants \
	--exclude-filtered \
	-V filtered_snps.vcf \
	--output bqsr_snps.vcf
gatk SelectVariants \
	--exclude-filtered \
	-V filtered_indels.vcf \
	--output bqsr_indels.vcf

printf "\n\nStep 9 - Base Quality Score Recalibration (BQSR) #1\n"
gatk BaseRecalibrator \
	-R $GENOME \
	-I sorted_dedup_reads.bam \
	--known-sites bqsr_snps.vcf \
	--known-sites bqsr_indels.vcf \
	--output recal_data.table

printf "\n\nStep 10 - Apply BQSR\n"
gatk ApplyBQSR \
	-R $GENOME \
	-I sorted_dedup_reads.bam \
	-bqsr recal_data.table \
	--output recal_reads.bam

printf "\n\nStep 11 - Base Quality Score Recalibration (BQSR) #2\n"
gatk BaseRecalibrator \
	-R $GENOME \
	-I recal_reads.bam \
	--known-sites bqsr_snps.vcf \
	--known-sites bqsr_indels.vcf \
	--output post_recal_data.table

printf "\n\nStep 12 - Analyze Covariates\n"
gatk AnalyzeCovariates \
	-before recal_data.table \
	-after post_recal_data.table \
	-plots recalibration_plots.pdf      # VISUALIZATION OPTION!

printf "\n\nStep 13-16  - Call Variants from recalibrated bam, extract, filter, and exclude\n"
gatk HaplotypeCaller \
	-R $GENOME \
	-ploidy 2 \
	-I recal_reads.bam \
	-O raw_variants_recal.vcf
gatk SelectVariants \
	-R $GENOME \
	-V raw_variants_recal.vcf \
	--select-type-to-include SNP \
	--output raw_snps_recal.vcf
gatk SelectVariants \
	-R $GENOME \
	-V raw_variants_recal.vcf \
	--select-type-to-include INDEL \
	--output raw_indels_recal.vcf
gatk VariantFiltration \
	-R $GENOME \
	-V raw_snps_recal.vcf \
	--output filtered_snps_final.vcf \
	-filter-name "QD_filter" -filter "QD < $snp_QD" \
	-filter-name "FS_filter" -filter "FS > $snp_FS" \
	-filter-name "MQ_filter" -filter "MQ < $snp_MQ" \
	-filter-name "SOR_filter" -filter "SOR > $snp_SOR" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < $snp_MQRankSum" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < $snp_ReadPosRankSum"
gatk VariantFiltration \
	-R $GENOME \
	-V raw_indels_recal.vcf \
	--output filtered_indels_final.vcf \
	-filter-name "QD_filter" -filter "QD < $indel_QD" \
	-filter-name "FS_filter" -filter "FS > $indel_FS" \
	-filter-name "SOR_filter" -filter "SOR > $indel_SOR"
gatk SelectVariants \
	--exclude-filtered \
	-V filtered_snps_final.vcf \
	--output bqsr2_snps.vcf
gatk SelectVariants \
	--exclude-filtered \
	-V filtered_indels_final.vcf \
	--output bqsr2_indels.vcf

printf "\n\nStep 17 - Annotate SNPs and Predict Effects\n"
# check available databases: java -jar $SNPEFF_JAR databases | grep -i <genus|species>
# download db first
java -jar $SNPEFF_JAR ann -v \
	$SNPEFFDB \
	bqsr2_snps.vcf > bqsr2_snps.ann.vcf



# setting up outputs for other steps:
#	chrom_length file for bam-bedgraph-bigwig.cwl
#		  chrom_length_file:
#		    type: File
#		    doc: "Tab delimited chromosome length file: <chromName><TAB><chromSize>"
fasta_rmlinebreaks.awk $GENOME | fasta_2col.sh | sed -e 's/^>//' -e 's/ .*\t/\t/' | awk -F'\t' '{printf("%s\t%s\n",$1,length($2))}' > chrom_length.tsv



#	OUTPUTS
#===============================================================================
printf "\n\nGenerating metrics and formatting overview.md file\n"
printf "\talignment metrics...\n"
read_count_r1=$(cut -f2 alignment_metrics.txt | tail -5 | head -1)
read_count_r2=$(cut -f2 alignment_metrics.txt | tail -4 | head -1)
read_count_total=$(cut -f2 alignment_metrics.txt | tail -3 | head -1)
mapped_count_r1=$(cut -f6 alignment_metrics.txt | tail -5 | head -1)
mapped_count_r2=$(cut -f6 alignment_metrics.txt | tail -4 | head -1)
mapped_count_total=$(cut -f6 alignment_metrics.txt | tail -3 | head -1)
mapped_perc_r1=$(cut -f7 alignment_metrics.txt | tail -5 | head -1)
mapped_perc_r2=$(cut -f7 alignment_metrics.txt | tail -4 | head -1)
mapped_perc_total=$(cut -f7 alignment_metrics.txt | tail -3 | head -1)

perc_coverage=$(cut -f1 coverage_metrics.tsv)
depth_mean=$(cut -f2 coverage_metrics.tsv)
depth_stdev=$(cut -f3 coverage_metrics.tsv)

reads_dedup=$(samtools view -c -F 0x4 -F 0x100 -F 0x800 sorted_dedup_reads.bam)
reads_recal=$(samtools view -c -F 0x4 -F 0x100 -F 0x800 recal_reads.bam)
printf "\tvariant counts...\n"
variants_raw_total=$(grep -c -v "^#" raw_variants.vcf)
variants_raw_snps=$(grep -c -v "^#" raw_snps.vcf)
variants_raw_indels=$(grep -c -v "^#" raw_indels.vcf)
variants_bqsr1_snps=$(grep -c -v "^#" bqsr_snps.vcf)     # bqsr contain "PASS" only snps after 1st recal/filtering/exclusion
variants_bqsr1_indels=$(grep -c -v "^#" bqsr_indels.vcf)     # bqsr contain "PASS" only indels after 1st recal/filtering/exclusion
variants_bqsr2_snps=$(grep -c -v "^#" bqsr2_snps.vcf)     # bqsr2 contain "PASS" only snps after 2nd recal/filtering/exclusion
variants_bqsr2_indels=$(grep -c -v "^#" bqsr2_indels.vcf)     # bqsr2 contain "PASS" only indels after 2nd recal/filtering/exclusion


printf "\tformatting...\n"

printf "## INPUTS\n" > overview.md
printf "\n" >> overview.md
printf "#### SECTION 1: general\n" >> overview.md
printf "-" >> overview.md
printf " \$THREADS, $THREADS\n" >> overview.md
printf "-" >> overview.md
printf " \$PLOIDY, $PLOIDY\n" >> overview.md
printf "-" >> overview.md
printf " \$GENOMEDIR, $GENOMEDIR\n" >> overview.md
printf "-" >> overview.md
printf " \$SNPEFFDB, $SNPEFFDB\n" >> overview.md
printf "-" >> overview.md
printf " \$READ1, $READ1\n" >> overview.md
printf "-" >> overview.md
printf " \$READ2, $READ2\n" >> overview.md
printf "\n" >> overview.md
printf "#### SECTION 2: snp filters\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_QD, $snp_QD\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_FS, $snp_FS\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_MQ, $snp_MQ\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_SOR, $snp_SOR\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_MQRankSum, $snp_MQRankSum\n" >> overview.md
printf "-" >> overview.md
printf " \$snp_ReadPosRankSum, $snp_ReadPosRankSum\n" >> overview.md
printf "\n" >> overview.md
printf "#### SECTION 3: indel filters\n" >> overview.md
printf "-" >> overview.md
printf " \$indel_QD, $indel_QD\n" >> overview.md
printf "-" >> overview.md
printf " \$indel_FS, $indel_FS\n" >> overview.md
printf "-" >> overview.md
printf " \$indel_SOR, $indel_SOR\n" >> overview.md
printf "\n" >> overview.md

printf "## READ & ALIGNMENT METRICS\n" >> overview.md
printf "-" >> overview.md
printf " Read 1 (R1) read count: $read_count_r1\n" >> overview.md
printf "-" >> overview.md
printf " Read 2 (R2) read count: $read_count_r2\n" >> overview.md
printf "-" >> overview.md
printf " Total read count: $read_count_total\n" >> overview.md
printf "-" >> overview.md
printf " R1 mapped reads: $mapped_count_r1\n" >> overview.md
printf "-" >> overview.md
printf " R2 mapped reads: $mapped_count_r2\n" >> overview.md
printf "-" >> overview.md
printf " Total mapped reads: $mapped_count_total\n" >> overview.md
printf "-" >> overview.md
printf " R1 percent mapped reads: $mapped_perc_r1\n" >> overview.md
printf "-" >> overview.md
printf " R2 percent mapped reads: $mapped_perc_r2\n" >> overview.md
printf "-" >> overview.md
printf " Total percent mapped reads: $mapped_perc_total\n" >> overview.md
printf "-" >> overview.md
printf "%s%s\n" "Percent Coverage: $perc_coverage " "%" >> overview.md
printf "-" >> overview.md
printf " Mean depth: $depth_mean\n" >> overview.md
printf "-" >> overview.md
printf " Stdev depth: $depth_stdev\n" >> overview.md
printf "-" >> overview.md
printf " Deduplicated reads: $reads_dedup\n" >> overview.md
printf "-" >> overview.md
printf " Base Quality Score Recalibrated (BQSR) reads: $reads_recal\n" >> overview.md
printf "\n" >> overview.md

printf "## VARIANT COUNTS\n" >> overview.md
printf "-" >> overview.md
printf " Total raw variants: $variants_raw_total\n" >> overview.md
printf "-" >> overview.md
printf " Raw SNPs: $variants_raw_snps\n" >> overview.md
printf "-" >> overview.md
printf " Raw Indels: $variants_raw_indels\n" >> overview.md
printf "-" >> overview.md
printf " BQSR1 SNPs: $variants_bqsr1_snps\n" >> overview.md
printf "-" >> overview.md
printf " BQSR1 Indels: $variants_bqsr1_indels\n" >> overview.md
printf "-" >> overview.md
printf " BQSR2 SNPs: $variants_bqsr2_snps\n" >> overview.md
printf "-" >> overview.md
printf " BQSR2 Indels: $variants_bqsr2_indels\n" >> overview.md



# compress output vcf files (better for IGV)
gzip raw_indels.vcf
gzip raw_snps.vcf
gzip bqsr2_indels.vcf
gzip bqsr2_snps.vcf
gzip bqsr2_snps.ann.vcf

printf "\n\nWorkflow script run_vc_germline.sh complete!\n"