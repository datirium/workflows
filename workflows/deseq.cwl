cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  rnaseq_cond_1:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se.cwl"
  rnaseq_cond_2:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  rpkm_isoforms_cond_1:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, e.g. 'treatment')"
    doc: "CSV/TSV input files grouped by isoforms (condition 1, e.g. 'treatment')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_isoforms"
    'sd:localLabel': true

  rpkm_genes_cond_1:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, e.g. 'treatment')"
    doc: "CSV/TSV input files grouped by genes (condition 1, e.g. 'treatment')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_genes"

  rpkm_common_tss_cond_1:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, e.g. 'treatment')"
    doc: "CSV/TSV input files grouped by common TSS (condition 1, e.g. 'treatment')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_common_tss"

  rpkm_isoforms_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, e.g. 'control')"
    doc: "CSV/TSV input files grouped by isoforms (condition 2, e.g. 'control')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_isoforms"
    'sd:localLabel': true

  rpkm_genes_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, e.g. 'control')"
    doc: "CSV/TSV input files grouped by genes (condition 2, e.g. 'control')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_genes"

  rpkm_common_tss_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, e.g. 'control')"
    doc: "CSV/TSV input files grouped by common TSS (condition 2, e.g. 'control')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_common_tss"

  group_by:
    type:
      - "null"
      - type: enum
        symbols: ["isoforms", "genes", "common tss"]
    default: "genes"
    label: "Group by"
    doc: "Grouping method for features: isoforms, genes or common tss"

  rpkm_cutoff:
    type: float?
    default: 0
    label: "Minimum rpkm cutoff. Applied before running DEseq"
    doc: "Minimum threshold for rpkm filtering. Default: 5"
    'sd:layout':
      advanced: true

  batch_file:
    type: File?
    default: null
    label: "Headerless TSV/CSV file for multi-factor analysis. First column - experiments' names from condition 1 and 2, second column - batch name"
    format: "http://edamontology.org/format_2330"
    doc: |
      Metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --ua and --ta, second column - batch name.
      Default: None

  alias_cond_1:
    type: string?
    default: "treatment"
    label: "Alias for condition 1, e.g. 'treatment' (letters and numbers only)"
    doc: "Name to be displayed for condition 1, e.g. 'treatment' (letters and numbers only)"
    'sd:layout':
      advanced: true

  alias_cond_2:
    type: string?
    default: "control"
    label: "Alias for condition 2, e.g. 'control' (letters and numbers only)"
    doc: "Name to be displayed for condition 2, e.g. 'control' (letters and numbers only)"
    'sd:layout':
      advanced: true

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
      - "none"
    default: "none"
    label: "Hopach clustering method to be run on normalized read counts"
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization analysis. Default: do not run clustering
    'sd:layout':
      advanced: true

  row_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "cosangle"
    label: "Distance metric for HOPACH row clustering"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided. Default: cosangle
    'sd:layout':
      advanced: true

  column_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "euclid"
    label: "Distance metric for HOPACH column clustering"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid
    'sd:layout':
      advanced: true

  center_row:
    type: boolean?
    default: false
    label: "Apply mean centering for feature expression prior to running clustering by row"
    doc: |
      Apply mean centering for feature expression prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered
    'sd:layout':
      advanced: true

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum P-adjusted to show features in the exploratory visualization analysis"
    doc: |
      In the exploratory visualization analysis output only features with
      adjusted P-value not bigger than this value. Default: 0.05
    'sd:layout':
      advanced: true

  sample_names_cond_1:
    type:
      - "null"
      - string[]
    default: null
    label: "Sample names for RNA-Seq experiments (condition 1, e.g. 'treatment')"
    doc: |
      Aliases for RNA-Seq experiments (condition 1, e.g. 'treatment') to make the
      legend for generated plots. Order corresponds to the rpkm_isoforms_cond_1
    'sd:upstreamSource': "rnaseq_cond_1/alias"

  sample_names_cond_2:
    type:
      - "null"
      - string[]
    default: null 
    label: "Sample names for RNA-Seq experiments (condition 2, e.g. 'control')"
    doc: |
      Aliases for RNA-Seq experiments (condition 2, e.g. 'control') to make the
      legend for generated plots. Order corresponds to the rpkm_isoforms_cond_2
    'sd:upstreamSource': "rnaseq_cond_2/alias"

  threads:
    type: int?
    default: 1
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  diff_expr_file:
    type: File
    label: "Differentially expressed features grouped by isoforms, genes or common TSS"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq generated file of differentially expressed features grouped by isoforms, genes or common TSS in TSV format"
    outputSource: deseq/diff_expr_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differential Expression Analysis'
        Title: 'Combined DESeq results'
    - scatter:
        tab: 'Volcano Plot'
        Title: 'Volcano'
        xAxisTitle: 'log fold change'
        yAxisTitle: '-log10(pAdj)'
        colors: ["#b3de69"]
        height: 600
        data: [$2, $9, $13]

  read_counts_file:
    type: File
    label: "Normalized read counts in GCT format. Compatible with GSEA"
    format: "http://edamontology.org/format_3709"
    doc: "DESeq generated file of with normalized read counts in GCT format. Compatible with GSEA"
    outputSource: deseq/read_counts_file

  phenotypes_file:
    type: File
    label: "Phenotype data file in CLS format. Compatible with GSEA"
    format: "http://edamontology.org/format_2330"
    doc: "DESeq generated file with phenotypes in CLS format. Compatible with GSEA"
    outputSource: deseq/phenotypes_file

  mds_plot_html:
    type: File?
    outputSource: deseq/mds_plot_html
    label: "MDS plot of normalized counts"
    doc: |
      MDS plot of normalized counts
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  plot_lfc_vs_mean:
    type: File?
    label: "Plot of normalised mean versus log2 fold change"
    format: "http://edamontology.org/format_3603"
    doc: |
      Plot of the log2 fold changes attributable to a given variable
      over the mean of normalized counts for all the samples
    outputSource: deseq/plot_lfc_vs_mean
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'LFC vs mean'

  gene_expr_heatmap:
    type: File?
    label: "Heatmap of the 30 most highly expressed features"
    format: "http://edamontology.org/format_3603"
    doc: |
      Heatmap showing the expression data of the 30 most highly expressed features grouped by
      isoforms, genes or common TSS, based on the variance stabilisation transformed data
    outputSource: deseq/gene_expr_heatmap
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'The 30 most highly expressed features'

  plot_pca:
    type: File?
    label: "PCA plot for variance stabilized count data"
    format: "http://edamontology.org/format_3603"
    doc: |
      PCA plot for variance stabilized count data. Values are now approximately
      homoskedastic (have constant variance along the range of mean values)
    outputSource: deseq/plot_pca
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'PCA plot for variance stabilized count data'

  plot_lfc_vs_mean_pdf:
    type: File?
    label: "Plot of normalised mean versus log2 fold change"
    format: "http://edamontology.org/format_3508"
    doc: |
      Plot of the log2 fold changes attributable to a given variable
      over the mean of normalized counts for all the samples
    outputSource: deseq/plot_lfc_vs_mean_pdf

  gene_expr_heatmap_pdf:
    type: File?
    label: "Heatmap of the 30 most highly expressed features"
    format: "http://edamontology.org/format_3508"
    doc: |
      Heatmap showing the expression data of the 30 most highly expressed features grouped by
      isoforms, genes or common TSS, based on the variance stabilisation transformed data
    outputSource: deseq/gene_expr_heatmap_pdf

  plot_pca_pdf:
    type: File?
    label: "PCA plot for variance stabilized count data"
    format: "http://edamontology.org/format_3508"
    doc: |
      PCA plot for variance stabilized count data. Values are now approximately
      homoskedastic (have constant variance along the range of mean values)
    outputSource: deseq/plot_pca_pdf

  volcano_plot_html_file:
    type: File
    outputSource: make_volcano_plot/html_file
    label: "Volcano Plot"
    doc: |
      HTML index file with volcano plot data.
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_css_file:
    type: File
    outputSource: make_volcano_plot/css_file
    label: "Volcano Plot CSS"
    doc: |
      CSS index file with volcano plot data.

  volcano_plot_js_file:
    type: File
    outputSource: make_volcano_plot/js_file
    label: "Volcano Plot JS"
    doc: |
      JS index file with volcano plot data.

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    label: "Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  deseq_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq stdout log"
    doc: "DESeq stdout log"
    outputSource: deseq/stdout_log

  deseq_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq stderr log"
    doc: "DESeq stderr log"
    outputSource: deseq/stderr_log

  morpheus_stdout_log:
    type: File
    outputSource: morpheus_heatmap/stdout_log
    label: "Morpheus heatmap stdout log"
    doc: "Morpheus heatmap stdout log"

  morpheus_stderr_log:
    type: File
    outputSource: morpheus_heatmap/stderr_log
    label: "Morpheus heatmap stderr log"
    doc: "Morpheus heatmap stderr log"


steps:

  deseq:
    run: ../tools/deseq-advanced.cwl
    in:
      untreated_files:
        source: [group_by, rpkm_isoforms_cond_1, rpkm_genes_cond_1, rpkm_common_tss_cond_1]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      treated_files:
        source: [group_by, rpkm_isoforms_cond_2, rpkm_genes_cond_2, rpkm_common_tss_cond_2]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      untreated_name: alias_cond_1
      treated_name: alias_cond_2
      untreated_sample_names: sample_names_cond_1
      treated_sample_names: sample_names_cond_2
      rpkm_cutoff: rpkm_cutoff
      batch_file: batch_file
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row: center_row
      maximum_padj: maximum_padj
      threads: threads
    out:
      - diff_expr_file
      - plot_lfc_vs_mean
      - gene_expr_heatmap
      - plot_pca
      - plot_lfc_vs_mean_pdf
      - gene_expr_heatmap_pdf
      - plot_pca_pdf
      - read_counts_file
      - phenotypes_file
      - mds_plot_html
      - stdout_log
      - stderr_log

  make_volcano_plot:
    run: ../tools/volcanot-plot.cwl
    in:
      diff_expr_file: deseq/diff_expr_file
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        source: group_by
        valueFrom: |
          ${
              if (self == "isoforms") {
                return "RefseqId";
              } else if (self == "genes") {
                return "GeneId";
              } else {
                return "GeneId";
              }
          }
    out:
      - html_file
      - css_file
      - js_file

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: deseq/read_counts_file
    out:
    - heatmap_html
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq - differential gene expression analysis"
label: "DESeq - differential gene expression analysis"
s:alternateName: "Differential gene expression analysis based on the negative binomial distribution"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/deseq.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
  - class: s:Organization
    s:legalName: "Datirium, LLC"
    s:member:
      - class: s:Person
        s:name: Artem BArski
        s:email: mailto:Artem.Barski@datirum.com
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@datirium.com
        s:sameAs:
          - id: http://orcid.org/0000-0001-9102-5681


doc: |
  # Differential gene expression analysis

  This differential gene expression (DGE) analysis takes as input samples from two experimental conditions that have been processed with an RNA-Seq workflow (see list of "Upstream workflows" below).

  DESeq estimates variance-mean dependence in count data from high-throughput sequencing assays, then tests for DGE based on a model which assumes a negative binomial distribution of gene expression (aligned read count per gene).

  ### Experimental Setup and Results Interpretation

  The workflow design uses as its fold change (FC) calculation: condition 1 (c1, e.g. treatment) over condition 2 (c2, e.g. control).
  
  In other words: `FC == (c1/c2)`

  Therefore:

  - if FC<1 the log2(FC) is <0 (negative), meaning expression in condition1<condition2 (gene is downregulated in c1)
  - if FC>1 the log2(FC) is >0 (positive), meaning expression in condition1>condition2 (gene is upregulated in c1)

  In other words, if you have input TREATMENT samples as condition 1, and CONTROL samples as condition 2, a positive L2FC for a gene indicates that expression of the gene in TREATMENT is greater (or upregulated) compared to CONTROL.
  
  Next, threshold the p-adjusted values with your FDR (false discovery rate) cutoff to determine if the change may be considered significant or not.

  It is important to note when DESeq1 or DESeq2 is used in our DGE analysis workflow. If a user inputs only a single sample per condition DESeq1 is used for calculating DGE. In this experimental setup, there are no repeated measurements per gene per condition, therefore biological variability in each condition cannot be captured so the output p-values are assumed to be purely "technical". On the other hand, if >1 sample(s) are input per condition DESeq2 is used. In this case, biological variability per gene within each condition is available to be incorporated into the model, and resulting p-values are assumed to be "biological". Additionally, DESeq2 fold change is "shrunk" to account for sample variability, and as Michael Love (DESeq maintainer) puts it, "it looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution. So the large fold changes from genes with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk. This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior".

  In either case, the null hypothesis (H0) tested is that there are no significantly differentially expressed genes between conditions, therefore a smaller p-value indicates a lower probability of the H0 occurring by random chance and therefore, below a certain threshold (traditionally <0.05), H0 should be rejected. Additionally, due to the many thousands of independent hypotheses being tested (each gene representing an independent test), the p-values attained by the Wald test are adjusted using the Benjamini and Hochberg method by default. These "padj" values should be used for determination of significance (a reasonable value here would be <0.10, i.e. below a 10% FDR).

  Further Analysis: Output from the DESeq workflow may be used as input to the GSEA (Gene Set Enrichment Analysis) workflow for identifying enriched marker gene sets between conditions.

  ### DESeq1

  High-throughput sequencing assays such as RNA-Seq, ChIP-Seq or barcode counting provide quantitative readouts in the form of count data. To infer differential signal in such data correctly and with good statistical power, estimation of data variability throughout the dynamic range and a suitable error model are required. Simon Anders and Wolfgang Huber propose a method based on the negative binomial distribution, with variance and mean linked by local regression and present an implementation, [DESeq](http://www.bioconductor.org/packages/3.8/bioc/html/DESeq.html), as an R/Bioconductor package.

  ### DESeq2

  In comparative high-throughput sequencing assays, a fundamental task is the analysis of count data, such as read counts per gene in RNA-seq, for evidence of systematic changes across experimental conditions. Small replicate numbers, discreteness, large dynamic range and the presence of outliers require a suitable statistical approach. [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html), a method for differential analysis of count data, using shrinkage estimation for dispersions and fold changes to improve stability and interpretability of estimates. This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.

  ### __References__
    - Anders S, Huber W (2010). “Differential expression analysis for sequence count data.” Genome Biology, 11, R106. doi: 10.1186/gb-2010-11-10-r106, http://genomebiology.com/2010/11/10/R106/.
    - Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
