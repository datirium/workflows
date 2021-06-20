# Available workflows
## RNA-Seq

- [RNA-Seq single-read](trim-rnaseq-se.md)
- [RNA-Seq paired-end](trim-rnaseq-pe.md)
- [RNA-Seq single-read strand specific](trim-rnaseq-se-dutp.md)
- [RNA-Seq paired-end strand specific](trim-rnaseq-pe-dutp.md)
- [RNA-Seq single-read (mitochondrial)](rnaseq-se-dutp-mitochondrial.md)
- [RNA-Seq paired-end (mitochondrial)](rnaseq-pe-dutp-mitochondrial.md)
- [RNA-Seq single-read strand specific (smarter)](trim-rnaseq-pe-smarter-dutp.md)

## ChIP-Seq

- [ChIP-Seq single-read](trim-chipseq-se.md)
- [ChIP-Seq paired-end](trim-chipseq-pe.md)

## ATAC-Seq

- [ATAC-Seq single-read](trim-atacseq-se.md)
- [ATAC-Seq paired-end](trim-atacseq-pe.md)

## Quant-Seq

- [QuantSeq 3' mRNA-Seq single-read](trim-quantseq-mrnaseq-se.md)
- [QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data](trim-quantseq-mrnaseq-se-strand-specific.md)

## Single Cell

- [Single-Cell Preprocessing Cell Ranger](single-cell-preprocess-cellranger.md)
- [Single-Cell Preprocessing](single-cell-preprocess.md)
- [Cellranger aggr - aggregates data from multiple Cellranger runs](cellranger-aggr.md)
- [Cell Ranger Build Reference Indices](cellranger-mkref.md)
- [Cellranger reanalyze - reruns secondary analysis performed on the feature-barcode matrix](cellranger-reanalyze.md)
- [Seurat for comparative scRNA-seq analysis of across experimental conditions](seurat-cluster.md)
- [AltAnalyze cell-level matching and comparison of single-cell transcriptomes](altanalyze-cellharmony.md)
- [AltAnalyze Iterative Clustering and Guide-gene Selection](altanalyze-icgs.md)
- [AltAnalyze Prepare Genome](altanalyze-prepare-genome.md)
- [SoupX - an R package for the estimation and removal of cell free mRNA contamination](soupx.md)

## Clip-Seq

- [CLIP-Seq single-read NNNNG](clipseq-se.md)

## Methylation

- [Bismark Methylation single-read](bismark-methylation-se.md)

## Cut-n-Run

- [Cut-n-Run paired-end](trim-chipseq-pe-cut-n-run.md)

## Differential expression

- [DESeq - differential gene expression analysis](deseq.md)
- [DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test](deseq-lrt.md)

## Differential binding analyses

- [DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data](diffbind.md)
- [THOR - differential peak calling of ChIP-seq signals with replicates](rgt-thor.md)
- [MAnorm SE - quantitative comparison of ChIP-Seq single-read data](manorm-se.md)
- [MAnorm PE - quantitative comparison of ChIP-Seq paired-end data](manorm-pe.md)

## Motif analyses

- [Motif Finding with HOMER with custom background regions](homer-motif-analysis-bg.md)
- [Motif Finding with HOMER with target and background regions from peaks](homer-motif-analysis-peak.md)
- [Motif Finding with HOMER with random background regions](homer-motif-analysis.md)

## Quality control

- [FastQC - a quality control tool for high throughput sequence data](fastqc.md)

## Gene set enrichment

- [GSEApy - Gene Set Enrichment Analysis in Python](gseapy.md)

## Enchancers analyses

- [ROSE: rank ordering of super-enhancers](super-enhancer.md)

## Dimensionality reduction, clustering, heatmaps

- [PCA - Principal Component Analysis](pca.md)
- [HOPACH - Hierarchical Ordered Partitioning and Collapsing Hybrid](hopach.md)
- [Tag density profile around regions of interest](heatmap.md)

## Genomic regions and lists manipulation

- [Genomic regions intersection and visualization](intervene.md)
- [Pairwise genomic regions intersection](peak-intersect.md)
- [Feature expression merge - combines feature expression from several experiments](feature-merge.md)
- [GAT - Genomic Association Tester](gat-run.md)
- [Interval overlapping alignments counts](bedtools-multicov.md)
- [Filter differentially expressed genes from DESeq for Tag Density Profile Analyses](filter-deseq-for-heatmap.md)
- [Filter ChIP/ATAC peaks for Tag Density Profile or Motif Enrichment analyses](filter-peaks-for-heatmap.md)

## Genome indices

- [Generate genome indices for STAR & bowtie](genome-indices.md)
- [Build STAR indices](star-index.md)
- [Build Bowtie indices](bowtie-index.md)
- [Build Bismark indices](bismark-index.md)

# Deprecated workflows

- rnaseq-se.cwl
- rnaseq-pe.cwl
- rnaseq-se-dutp.cwl
- rnaseq-pe-dutp.cwl
- chipseq-se.cwl
- chipseq-pe.cwl