#!/bin/bash

#   Shell wrapper for running GSEApy and associated scripts
#
##########################################################################################
#
# v1.0.0
# - 20240910, initial release
#
##########################################################################################
printf "$(date)\nLog file for run_gseapy.sh\n\n"

#	FUNCTIONS
#===============================================================================
usage()
{
cat << EOF
Help message for \`run_genelists.sh\`:
Shell wrapper for GSEApy and associated scripts


Primary Output files:
 - gseapy.gsea.gene_set.report.csv, table of gene set enrichment scores and pvalues in CSV format
 - enrichment plots and heatmaps per gene set, depending on user-selected dataset

Secondary Output files:
 - reportsummary.md, markdown with general filtering and analysis statistics

PARAMS:
    SECTION 1: general
	-h	help		show this message
	-t  INT			number of threads
	-c  FILE		Input class vector (phenotype) file in CLS format (from DESeq workflow)
	-d	FILE		Input gene expression dataset file in txt or gct format (from DESeq workflow)
	-g	FILE		Gene set database, from prefilled dropdown or use-provided .gmt file
	-r	STRING		Permutation type. Default: "gene_set"
	-n	INT			Number of random permutations. For calculating esnulls. Default: 1000
	-w	INT			Min size of input genes presented in Gene Sets. Default: 15
	-x	INT			Max size of input genes presented in Gene Sets. Default: 500
	-m	STRING		Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes
	-s	INT			Number of random seed. Default: None
	-p	FLOAT		Output only graphs from gene sets with less than thsi set p-value. Default: 1.00
    
EOF
}


#	INPUTS & CHECKS & DEFAULTS
#===============================================================================
# parse args
while getopts "ht:c:d:g:r:n:w:x:m:s:p:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
		c) PHENOFILE=$OPTARG ;;
		d) EXPDATA=$OPTARG ;;
		g) DATASETIN=$OPTARG ;;
		r) PERMTYPE=$OPTARG ;;
		n) PERMNUM=$OPTARG ;;
		w) GENEFILTERMIN=$OPTARG ;;
		x) GENEFILTERMAX=$OPTARG ;;
		m) CORRMETRICS=$OPTARG ;;
		s) RANDSEED=$OPTARG ;;
		p) PVALUE=$OPTARG ;;
		?) usage; exit ;;
	esac
done
# arg checks
if [[ "$THREADS" == "" ]]; then THREADS=1; fi
if [[ "$PHENOFILE" == "" ]]; then echo "error: required param missing (-c), FILE, Input class vector (phenotype) file in CLS format (from DESeq workflow)"; exit; fi
if [[ "$EXPDATA" == "" ]]; then echo "error: required param missing (-d), FILE, Input gene expression dataset file in txt or gct format (from DESeq workflow)"; exit; fi
if [[ "$DATASETIN" == "" ]]; then echo "error: required param missing (-g), FILE, Gene set database, from prefilled dropdown or use-provided .gmt file"; exit; fi
if [[ "$PERMTYPE" == "" ]]; then echo "warning: optional param missing (-r), STRING, Permutation type. Default: gene_set"; fi
if [[ "$PERMNUM" == "" ]]; then echo "warning: optional param missing (-n), INT, Number of random permutations. For calculating esnulls. Default: 1000"; fi
if [[ "$GENEFILTERMIN" == "" ]]; then echo "warning: optional param missing (-w), INT, Min size of input genes presented in Gene Sets. Default: 15"; fi
if [[ "$GENEFILTERMAX" == "" ]]; then echo "warning: optional param missing (-x), INT, Max size of input genes presented in Gene Sets. Default: 500"; fi
if [[ "$CORRMETRICS" == "" ]]; then echo "warning: optional param missing (-m), STRING, Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes"; fi
if [[ "$RANDSEED" == "" ]]; then echo "warning: optional param missing (-s), INT, Number of random seed. Default: None"; fi
if [[ "$PVALUE" == "" ]]; then echo "warning: optional param missing (-p), FLOAT, Output only graphs from gene sets with less than thsi set p-value. Default: 0.05"; fi


#	MAIN
#===============================================================================
printf "Selecting dataset...\n"
if [[ "$DATASETIN" == "H_hallmark_gene_sets" ]]; then
	DATASET="/opt/gseapy/H_hallmark_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C1_positional_gene_sets" ]]; then
	DATASET="/opt/gseapy/C1_positional_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C2_curated_gene_sets" ]]; then
	DATASET="/opt/gseapy/C2_curated_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C3_regulatory_target_gene_sets" ]]; then
	DATASET="/opt/gseapy/C3_regulatory_target_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C4_computational_gene_sets" ]]; then
	DATASET="/opt/gseapy/C4_computational_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C5_ontology_gene_sets" ]]; then
	DATASET="/opt/gseapy/C5_ontology_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C6_oncogenic_signature_gene_sets" ]]; then
	DATASET="/opt/gseapy/C6_oncogenic_signature_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C7_immunologic_signature_gene_sets" ]]; then
	DATASET="/opt/gseapy/C7_immunologic_signature_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "C8_cell_type_signature_gene_sets" ]]; then
	DATASET="/opt/gseapy/C8_cell_type_signature_gene_sets.v2024.1.gmt";
elif [[ "$DATASETIN" == "KEGG_2021_Human" ]]; then
	DATASET="/opt/gseapy/enrichr.KEGG_2021_Human.gmt";
elif [[ "$DATASETIN" == "Reactome_2022" ]]; then
	DATASET="/opt/gseapy/enrichr.Reactome_2022.gmt";
elif [[ "$DATASETIN" == "WikiPathway_2023_Human" ]]; then
	DATASET="/opt/gseapy/enrichr.WikiPathway_2023_Human.gmt";
else
	DATASET=$DATASETIN;
fi



printf "Running GSEApy program...\n"
gseapy gsea \
	--threads $THREADS \
	--cls "$PHENOFILE" \
	--data "$EXPDATA" \
	--gmt "$DATASET" \
	--permu-type "$PERMTYPE" \
	--permu-num "$PERMNUM" \
	--min-size $GENEFILTERMIN \
	--max-size $GENEFILTERMAX \
	--method "$CORRMETRICS" \
	--seed $RANDSEED \
	--graph 1000000


printf "filtering plots on p-value..."
report="GSEApy_reports/gseapy.gsea.gene_set.report.csv"
cut -d',' -f 1-4 $report | tail -n+2 | awk -F"," -v pval=$PVALUE '{if($4<pval){print()}}' > gseapy_report_filtered.csv
# aggregate and compress filtered plots
mkdir filtered_enrichment_plots filtered_heatmap_plots
cut -d',' -f1 gseapy_report_filtered.csv | while read term; do find GSEApy_reports/ -name "*$term*gsea.pdf"; done | while read enrich; do mv $enrich filtered_enrichment_plots/; done
cut -d',' -f1 gseapy_report_filtered.csv | while read term; do find GSEApy_reports/ -name "*$term*heatmap.pdf"; done | while read enrich; do mv $enrich filtered_heatmap_plots/; done
tar -cf filtered_enrichment_plots.tar filtered_enrichment_plots/
tar -cf filtered_heatmap_plots.tar filtered_heatmap_plots/
gzip filtered_enrichment_plots.tar filtered_heatmap_plots.tar

# aggregate and compress all plots
tar -cf all_gseapy_enrichment_plots.tar GSEApy_reports/*.gsea.pdf
tar -cf all_gseapy_heatmap_plots.tar GSEApy_reports/*.heatmap.pdf
gzip all_gseapy_enrichment_plots.tar all_gseapy_heatmap_plots.tar



# generate report summary TSV file for kableExtra
report="GSEApy_reports/gseapy.gsea.gene_set.report.csv"
log=GSEApy_reports/gseapy.*.log

#   DATASET DETAILS
# total unique genes in rna-seq analysis (remove first 3 rows of header info)
total_unique_genes=$(cut -f1 $EXPDATA | tail -n+4 | sort | uniq | awk 'END{print(NR)}')
# total unique genes among all leading edge genes of all marker datasets
ledge_genes=$(cut -d$',' -f9 $report | tail -n+2 | sed 's/;/\n/g' | sort | uniq | awk 'END{print(NR)}')
# get total number of gene sets in selected dataset
total_gene_sets=$(grep -i "gene_sets" $log | sed -e 's/.*] \+//' -e 's/ .*//' | awk '{x+=$0+0}END{print x}')
# get number of filtered gene sets
total_filtered=$(grep -i "gene_sets" $log | sed -e 's/.*] \+//' -e 's/ .*//' | awk '{x=$0+0;print x}' | head -1)
# get number of gene sets used for statistical testing
total_analyzed=$(grep -i "gene_sets" $log | sed -e 's/.*] \+//' -e 's/ .*//' | awk '{x=$0+0;print x}' | tail -1)
# total unique genes in remaining marker gene datasets
feature_genes=$(cut -d$',' -f8 $report | tail -n+2 | sed 's/;/\n/g' | sort | uniq | awk 'END{print(NR)}')

#   FORMATTING OUTPUT
# header for kable
printf "%s\t%s\n" "GSEApy Result" "Value" > reportsummary.tsv

# totals
dataset_name=$(basename "$DATASET" | sed 's/\.gmt//')
printf "Total unique genes from input experiment data\t$total_unique_genes\n" >> reportsummary.tsv
printf "Leading edge subset genes\t$ledge_genes\n" >> reportsummary.tsv
printf "Total gene sets in $dataset_name\t$total_gene_sets\n" >> reportsummary.tsv
printf "Total gene sets remaining after min/max filtering\t$total_analyzed\n" >> reportsummary.tsv
printf "Total features (genes) analyzed in $dataset_name\t$feature_genes\n" >> reportsummary.tsv

#   PHENOTYPE AND ANALYSIS DETAILS
# get phenotype names
p1=$(tail -1 $PHENOFILE | awk -F' ' '{print($1)}')
p2=$(tail -1 $PHENOFILE | awk -F' ' '{print($NF)}')
#echo $p1 $p2
# phenotype 1
enriched=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1>0){x++}}END{print(x)}')
sigsets_fdr25=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1>0){if($4<0.25){x++}}}END{print(x)}')
sigsets_pv5=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1>0){if($3<0.05){x++}}}END{print(x)}')
sigsets_pv1=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1>0){if($3<0.01){x++}}}END{print(x)}')
printf "Gene sets enriched in phenotype $p1\t$enriched\n" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets at FDR<25%" "$sigsets_fdr25" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets w/ p-value<5%" "$sigsets_pv5" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets w/ p-value <1%" "$sigsets_pv1" >> reportsummary.tsv
# phenotype 2
enriched=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1<0){x++}}END{print(x)}')
sigsets_fdr25=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1<0){if($4<0.25){x++}}}END{print(x)}')
sigsets_pv5=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1<0){if($3<0.05){x++}}}END{print(x)}')
sigsets_pv1=$(cut -d$',' -f2-5 $report | tail -n+2 | awk -F',' 'BEGIN{x=0};{if($1<0){if($3<0.01){x++}}}END{print(x)}')
printf "Gene sets enriched in phenotype $p2\t$enriched\n" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets at FDR<25%" "$sigsets_fdr25" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets w/ p-value<5%" "$sigsets_pv5" >> reportsummary.tsv
printf "%s\t%s\n" "Gene sets w/ p-value <1%" "$sigsets_pv1" >> reportsummary.tsv

# generate formatted markdown table
Rscript /usr/local/bin/run_tsv_to_kable.R reportsummary.tsv reportsummary.md



printf "\n\nWorkflow script run_gseapy.sh complete!\n"