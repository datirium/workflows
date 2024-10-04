cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  rnaseq_cond_1:
  - "trim-rnaseq-pe.cwl"
  - "trim-rnaseq-se.cwl"
  - "trim-rnaseq-pe-dutp.cwl"
  - "trim-rnaseq-pe-smarter-dutp.cwl"
  - "trim-rnaseq-se-dutp.cwl"
  - "trim-quantseq-mrnaseq-se-strand-specific.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-smarter-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-quantseq-mrnaseq-se-strand-specific.cwl"
  rnaseq_cond_2:
  - "trim-rnaseq-pe.cwl"
  - "trim-rnaseq-se.cwl"
  - "trim-rnaseq-pe-dutp.cwl"
  - "trim-rnaseq-pe-smarter-dutp.cwl"
  - "trim-rnaseq-se-dutp.cwl"
  - "trim-quantseq-mrnaseq-se-strand-specific.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-smarter-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se-dutp.cwl"
  - "https://github.com/datirium/workflows/workflows/trim-quantseq-mrnaseq-se-strand-specific.cwl"


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
    label: "RNA-Seq experiments (condition 1, aka 'untreated')"
    doc: "CSV/TSV input files grouped by isoforms (condition 1, aka 'untreated')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_isoforms"
    'sd:localLabel': true

  rpkm_genes_cond_1:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, aka 'untreated')"
    doc: "CSV/TSV input files grouped by genes (condition 1, aka 'untreated')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_genes"

  rpkm_common_tss_cond_1:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 1, aka 'untreated')"
    doc: "CSV/TSV input files grouped by common TSS (condition 1, aka 'untreated')"
    'sd:upstreamSource': "rnaseq_cond_1/rpkm_common_tss"

  rpkm_isoforms_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, aka 'treated')"
    doc: "CSV/TSV input files grouped by isoforms (condition 2, aka 'treated')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_isoforms"
    'sd:localLabel': true

  rpkm_genes_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, aka 'treated')"
    doc: "CSV/TSV input files grouped by genes (condition 2, aka 'treated')"
    'sd:upstreamSource': "rnaseq_cond_2/rpkm_genes"

  rpkm_common_tss_cond_2:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments (condition 2, aka 'treated')"
    doc: "CSV/TSV input files grouped by common TSS (condition 2, aka 'treated')"
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
    default: "untreated"
    label: "Alias for condition 1, aka 'untreated' (letters and numbers only)"
    doc: "Name to be displayed for condition 1, aka 'untreated' (letters and numbers only)"
    'sd:layout':
      advanced: true

  alias_cond_2:
    type: string?
    default: "treated"
    label: "Alias for condition 2, aka 'treated' (letters and numbers only)"
    doc: "Name to be displayed for condition 2, aka 'treated' (letters and numbers only)"
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
    label: "Sample names for RNA-Seq experiments (condition 1, aka 'untreated')"
    doc: |
      Aliases for RNA-Seq experiments (condition 1, aka 'untreated') to make the
      legend for generated plots. Order corresponds to the rpkm_isoforms_cond_1
    'sd:upstreamSource': "rnaseq_cond_1/alias"

  sample_names_cond_2:
    type:
      - "null"
      - string[]
    default: null 
    label: "Sample names for RNA-Seq experiments (condition 2, aka 'treated')"
    doc: |
      Aliases for RNA-Seq experiments (condition 2, aka 'treated') to make the
      legend for generated plots. Order corresponds to the rpkm_isoforms_cond_2
    'sd:upstreamSource': "rnaseq_cond_2/alias"

  threads:
    type: int?
    default: 6
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
      HTML index file for Volcano Plot
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    outputSource: make_volcano_plot/html_data
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for Volcano Plot

  ma_plot_html_file:
    type: File
    outputSource: make_ma_plot/html_file
    label: "MA-plot"
    doc: |
      HTML index file for MA-plot
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  ma_plot_html_data:
    type: Directory
    outputSource: make_ma_plot/html_data
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for MA-plot

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
    run: ../tools/deseq-advanced-deprecated.cwl
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
    run: ../tools/volcano-plot.cwl
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
    - html_data
    - html_file

  make_ma_plot:
    run: ../tools/ma-plot.cwl
    in:
      diff_expr_file: deseq/diff_expr_file
      x_axis_column:
        default: "baseMean"
      y_axis_column:
        default: "log2FoldChange"
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
    - html_data
    - html_file

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

s:name: "Deprecated. DESeq - differential gene expression analysis"
label: "Deprecated. DESeq - differential gene expression analysis"
s:alternateName: "Deprecated. Differential gene expression analysis based on the negative binomial distribution"

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


# doc:
#   $include: ../descriptions/deseq.md


doc: |
  Deprecated. Differential gene expression analysis
  =====================================

  Differential gene expression analysis based on the negative binomial distribution

  Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.

  DESeq1
  ------

  High-throughput sequencing assays such as RNA-Seq, ChIP-Seq or barcode counting provide quantitative readouts
  in the form of count data. To infer differential signal in such data correctly and with good statistical power,
  estimation of data variability throughout the dynamic range and a suitable error model are required.
  Simon Anders and Wolfgang Huber propose a method based on the negative binomial distribution, with variance and mean
  linked by local regression and present an implementation, [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html),
  as an R/Bioconductor package 

  DESeq2
  ------

  In comparative high-throughput sequencing assays, a fundamental task is the analysis of count data,
  such as read counts per gene in RNA-seq, for evidence of systematic changes across experimental conditions.
  Small replicate numbers, discreteness, large dynamic range and the presence of outliers require a
  suitable statistical approach. [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html),
  a method for differential analysis of count data,
  using shrinkage estimation for dispersions and fold changes to improve stability and interpretability of estimates.
  This enables a more quantitative analysis focused on the strength rather than the mere presence of differential expression.