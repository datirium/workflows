#!/bin/bash

#   Shell wrapper for run_rnbeads_diff.R
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
    Wrapper for RnBeads differential methylation pipeline, with downstream processing for tables and IGV

    Primary Output files:
     - sig_dm_sites.bed (bed for IGV; sig diff meth sites)
     - sig_dm_sites_annotated.tsv (tsv for TABLE; for each site above, closest single gene annotation)

    Output rnbeads reports directory in container at '/tmp/reports/', includes:
    reports/
    ├── data_import.html
    ├── differential_methylation.html
    ├── preprocessing.html
    ├── quality_control.html
    ├── tracks_and_tables.html

PARAMS:
 -h  help	show this message
 -g  STRING   Sample genome, available options: hg19, hg38, mm9, mm10, rn5
 -t  INT	number of threads
 -a  STRING     name of condition1
 -b  STRING     name of condition2
 -c  LIST	comma separated list of absolute filepaths to all condition1 bed files (BismarkCov format)
 -d  LIST	comma separated list of absolute filepaths to all condition2 bed files (BismarkCov format)
 -j  LIST   comma separated list of sample names in condition1
 -k  LIST   comma separated list of sample names in condition2
 -m  FILE   refGene.txt file for annotating DM sites with gene information
 -p  LIST   comma separated list of absolute filepaths to all condition1 bigWig files (from Bismark upstream)
 -q  LIST   comma separated list of absolute filepaths to all condition2 bigWig files (from Bismark upstream)

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
while getopts "ht:g:a:b:c:d:j:k:m:p:q:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		t) THREADS=$OPTARG ;;
        g) GENOME=$OPTARG ;;
		a) CONDITION1_NAME=$OPTARG ;;
		b) CONDITION2_NAME=$OPTARG ;;
		c) CONDITION1_BED_FILEPATHS=$OPTARG ;;
		d) CONDITION2_BED_FILEPATHS=$OPTARG ;;
		j) CONDITION1_ALIASES=$OPTARG ;;
        k) CONDITION2_ALIASES=$OPTARG ;;
        m) REFGENETXT=$OPTARG ;;
        p) CONDITION1_BIGWIG_FILEPATHS=$OPTARG ;;
		q) CONDITION2_BIGWIG_FILEPATHS=$OPTARG ;;
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
printf "\tCONDITION2_NAME - $CONDITION2_NAME\n"
printf "\tCONDITION1_BED_FILEPATHS - $CONDITION1_BED_FILEPATHS\n"
printf "\tCONDITION2_BED_FILEPATHS - $CONDITION2_BED_FILEPATHS\n"
printf "\tCONDITION1_ALIASES - $CONDITION1_ALIASES\n"
printf "\tCONDITION2_ALIASES - $CONDITION2_ALIASES\n"
printf "\REFGENETXT - $REFGENETXT\n"
printf "\tCONDITION1_BIGWIG_FILEPATHS - $CONDITION1_BIGWIG_FILEPATHS\n"
printf "\tCONDITION2_BIGWIG_FILEPATHS - $CONDITION2_BIGWIG_FILEPATHS\n"
printf "\tOUTDIR - $OUTDIR\n"


#	MAIN
#===============================================================================
# copy all data into workdir, and use absolute path as rnbeads data.source[1] input vector
#   also generate sample_annotation.csv
printf "Sample_ID,condition\n" > $workdir/sample_annotation.csv
echo "$CONDITION1_BED_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    cp "$filepath" "$workdir/"
    printf "%s,%s\n" $(basename "$filepath") "$CONDITION1_NAME";
done >> $workdir/sample_annotation.csv
echo "$CONDITION2_BED_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    cp "$filepath" "$workdir/"
    printf "%s,%s\n" $(basename "$filepath") "$CONDITION2_NAME";
done >> $workdir/sample_annotation.csv
# check sample sheet csv
cat $workdir/sample_annotation.csv

# copy and rename (force ".bigWig" file ext) bismark bigwig upstreams for IGV visualization
rep=1; echo "$CONDITION1_BIGWIG_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    bn=$(basename "$filepath" | sed 's/\.bigWig//')
    cp "$filepath" "c1_rep${rep}_${bn}.bigWig"
    ((rep++))
done
rep=1; echo "$CONDITION2_BIGWIG_FILEPATHS" | sed 's/,/\n/g' | while read filepath; do
    bn=$(basename "$filepath" | sed 's/\.bigWig//')
    cp "$filepath" "c2_rep${rep}_${bn}.bigWig"
    ((rep++))
done

# run Rscript wrapper for rnbeads with formatted inputs
Rscript /usr/local/bin/run_rnbeads_diff.R $workdir/sample_annotation.csv "$OUTDIR" "$workdir" "$GENOME" "$THREADS"

#	OUTPUTS
#===============================================================================
# move sample sheet for output
cp $workdir/sample_annotation.csv ./
# format for Overview tab (use sample alias instead of filename)
head -1 sample_annotation.csv | awk -F',' '{printf("| %s | %s |\n",$1,$2)}' > sample_annotation.md
printf "| -- | -- |\n" >> sample_annotation.md
echo "$CONDITION1_ALIASES" | sed 's/,/\n/g' | while read alias; do
    echo "$alias" | awk -v x="$CONDITION1_NAME" '{printf("| %s | %s |\n",$0,x)}'
done >> sample_annotation.md
echo "$CONDITION2_ALIASES" | sed 's/,/\n/g' | while read alias; do
    echo "$alias" | awk -v x="$CONDITION2_NAME" '{printf("| %s | %s |\n",$0,x)}'
done >> sample_annotation.md

# get annotationed DMs (differential methylation) data (https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/RnBeads/inst/doc/RnBeads_Annotations.pdf)
#   2.1 Sites ($sites)
#Currently, every data package examines cytosines in the context of CpG and contains an
#annnotation table of all CpGs in the respective genome. CpG density and GC content are
#also computed for the neighborhood of length 100 base pairs centered on each locus. The
#total number of dinucleotides annotated in HG19 is 28,217,009 represented both on the
#forward and reverse DNA strands.
#   2.4 Regions ($tiling, $cpg, $genes)
#Every data package defines the following sets of regions for the dedicated assembly:
# - GpG islands The CpG island track is downloaded from the dedicated FTP directory of the UCSC Genome Browser.
# - Tiling regions Tiling regions with a window size of 5 kilobases are defined over the whole genome.
# - Genes and promoters Ensembl3 gene definitions are downloaded using the biomaRt package. A promoter is defined as the region spanning 1,500 bases upstream and 500 bases downstream of the transcription start site of the corresponding gene.
#CpG density and GC content are computed for all region types listed above.
dm="./reports/differential_methylation_data"
genes="$dm/diffMethTable_region_cmp1_genes.csv"
promoters="$dm/diffMethTable_region_cmp1_promoters.csv"
cpg="$dm/diffMethTable_region_cmp1_cpgislands.csv"
tiling="$dm/diffMethTable_region_cmp1_tiling.csv"
sites="$dm/diffMethTable_site_cmp1.csv"

# FOR DOWNLOAD
#   quot.log2 - log2 of the quotient in methylation: log2(mean.g1+ε/mean.g2+ε), where ε = 0.01. In case of paired analysis, it is the mean of the pairwise quotients
#   FOR TABLE VIEW (& DOWNLOAD): id,Chromosome,Start,End,symbol,entrezID,comb.p.adj.fdr,num.sites,mean.mean.covg.condition1,mean.mean.covg.condition2
sed 's/,/\t/g' $genes | cut -f1-6,10,12,14,17,18 > dm_genes.tsv &
sed 's/,/\t/g' $promoters | cut -f1-6,10,12,14,17,18 > dm_promoters.tsv &
#   $cpg, $tiling: id,Chromosome,Start,End,mean.quot.log2,comb.p.adj.fdr,num.sites,mean.mean.covg.condition1,mean.mean.covg.condition2
sed 's/,/\t/g' $cpg | cut -f1-4,8,10,12,15,16 > dm_cpg.tsv &
sed 's/,/\t/g' $tiling | cut -f1-4,8,10,12,15,16 > dm_tiling.tsv &
#   $sites: id,Chromosome,Start,Strand,mean.quot.log2,diffmeth.p.adj.fdr,mean.covg.condition1,mean.covg.condition2
sed 's/,/\t/g' $sites | cut -f1-4,8,17,21,22 > dm_sites.tsv &
wait

# FOR IGV (sig diff meth with quotient of meth log2, FDR, and score)
#    reduce to those with dm p.adj<0.10, and output as bed format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#       1-3 chr,start,start+1
#       4   quot.meth.log2=x;FDR=y (semicolon separated)
#       5   score (transform 0.000-0.1000 to value between 0-1000, a lower padj will be a higher score - closer to 1000)
#           `score=1000-(1000/($7[cov1]-$8[cov2])^2)`, if diff==0, score=0, elif score<1==1, else score
#           function: https://www.wolframalpha.com/input?i=solve+1000-%281000%2F%28x%5E2%29%29
#       6   strand
#   Sites
awk -F'\t' '{if($6<0.10){if(($7-$8)==0){score=0}else{y=1000-(1000/(($7-$8)^2)); if(y<1){score=1}else{score=y}}; printf("%s\t%s\t%s\tquot.meth.log2=%.2f;FDR=%.2f\t%.0f\t%s\n",$2,$3,$3+1,$5,$6,score,$4)}}' dm_sites.tsv > sig_dm_sites.bed &
#   CpG Islands
awk -F'\t' '{if($6<0.10){if(($8-$9)==0){score=0}else{y=1000-(1000/(($8-$9)^2)); if(y<1){score=1}else{score=y}}; printf("%s\t%s\t%s\tquot.meth.log2=%.2f;FDR=%.2f\t%.0f\t%s\n",$2,$3,$4,$5,$6,score,"=")}}' dm_cpg.tsv > sig_dm_cpg.bed &


# FOR TABLES (annotated closest gene to each DM Site or CpG Island - separate files)
#    reformat refGene.tsv into bed
#    get headers from ucsc: wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.sql
awk -F'\t' '{printf("%s\t%s\t%s\t%s\n",$3,$5,$6,$2)}' <(tail -n+2 $REFGENETXT) > refGene.bed &
wait
#    prep gene annotations for matching
#         reduce refGene.txt to single row per Refseq_id
#         some ids have more than 1 location (e.g. `grep -P "\tNR_125730\t" refGene.txt | less` shows 9 different chr locations)
#         split these ids and annotate `loc_multi` genes with "multi" for the table
grep -v "_alt" $REFGENETXT | grep -v "_fix" | cut -f2 | sort | uniq -c | sort -rnk1 | sed 's/^ \+//' | awk -F' ' '{if($1==1){print($2)}}' > genes.loc_single &
grep -v "_alt" $REFGENETXT | grep -v "_fix" | cut -f2 | sort | uniq -c | sort -rnk1 | sed 's/^ \+//' | awk -F' ' '{if($1>1){print($2)}}' > genes.loc_multi &
wait
#         output separate annotation files for single (all) and multi (refSeq_id and Gene_id [symbol] only)
awk -F'\t' '{if(FNR==NR){gene_id[$2]=$13; chr[$2]=$3; start[$2]=$5; end[$2]=$6; strand[$2]=$4}else{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$0,gene_id[$0],chr[$0],start[$0],end[$0],strand[$0])}}' $REFGENETXT genes.loc_single > refGene.reduced.txt
awk -F'\t' '{if(FNR==NR){gene_id[$2]=$13}else{printf("%s\t%s\n",$0,gene_id[$0])}}' $REFGENETXT genes.loc_multi >> refGene.reduced.txt
# find closest single gene to each DM (don't forget to sort beds!), NOTE: .cfs extension, closest-feature shortest format
closest-features --shortest <(sort-bed sig_dm_sites.bed) <(sort-bed refGene.bed) > sig_dm_sites.refGene.cfs
closest-features --shortest <(sort-bed sig_dm_cpg.bed) <(sort-bed refGene.bed) > sig_dm_cpg.refGene.cfs
# merge gene annotations with dms for table output
#    output table headers
#         dm_<sites|cpg>.tsv - $1,$2,$3,$4,$5,$6,$7,$8 (all columns)
#         refGene.txt - $2,$13,$3,$5,$6,$4
printf "Site_id\tSite_Chr\tSite_position\tSite_strand\tLog2_Meth_Quotient\tFDR\tCoverage_score\tmeancov_${CONDITION1_NAME}\tmeancov_${CONDITION2_NAME}\trefSeq_id\tGene_id\tChr\ttxStart\ttsEnd\ttxStrand\n" > sig_dm_sites_annotated.tsv
printf "CpGisland_id\tCpGisland_Chr\tCpGisland_position\tCpGisland_strand\tLog2_Meth_Quotient\tFDR\tCoverage_score\tmeancov_${CONDITION1_NAME}\tmeancov_${CONDITION2_NAME}\trefSeq_id\tGene_id\tChr\ttxStart\ttsEnd\ttxStrand\n" > sig_dm_cpg_annotated.tsv
#    reformat sig dm files for pretty tabling and annotation matching
sed -e 's/quot.meth.log2=//' -e 's/;FDR=/\t/' -e 's/|/\t/' sig_dm_sites.refGene.cfs > tmp1.sites
sed -e 's/quot.meth.log2=//' -e 's/;FDR=/\t/' -e 's/|/\t/' sig_dm_cpg.refGene.cfs > tmp1.cpg
#    match all sig_dm with statistics (file1 match chr and start (col2,3) with same in file 2 (col1,2)), link to refseq_id (file2 $6), reorder per "table headers" (above)
awk -F'\t' '{if(FNR==NR){site_id[$2,$3]=$1; cov1[$1]=$7; cov2[$1]=$8}else{id=site_id[$1,$2]; printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",id,$1,$2,$7,$4,$5,$6,cov1[id],cov2[id],$11)}}' dm_sites.tsv tmp1.sites > tmp2.sites
awk -F'\t' '{if(FNR==NR){cpg_id[$2,$3]=$1; cov1[$1]=$8; cov2[$1]=$9}else{id=cpg_id[$1,$2]; printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",id,$1,$2,$7,$4,$5,$6,cov1[id],cov2[id],$11)}}' dm_cpg.tsv tmp1.cpg > tmp2.cpg
#    store reduced annotations, match to sig dm refseq_ids annotation (single location reported from bedops) with associated dm data
awk -F'\t' '{if(FNR==NR){anno[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6}else{printf("%s\t%s\n",$0,anno[$10])}}' refGene.reduced.txt tmp2.sites >> sig_dm_sites_annotated.tsv
awk -F'\t' '{if(FNR==NR){anno[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6}else{printf("%s\t%s\n",$0,anno[$10])}}' refGene.reduced.txt tmp2.cpg >> sig_dm_cpg_annotated.tsv


# package full report dir
tar -cf reports.tar ./reports
gzip reports.tar

# minimize report dir (for CWL output dir fix)
rm -r ./reports/tracks_and_tables*
for minimizer in "data_import" "differential_methylation" "preprocessing" "quality_control"; do
    rm -r ./reports/${minimizer}_data
    rm -r ./reports/${minimizer}_pdfs
done
find ./reports -name "*_high_resolution.png" -exec rm {} +

# clean up
rm -rf "$workdir"