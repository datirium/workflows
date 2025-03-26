cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.42


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This
      file should include genes expression information
      stored in the RNA assay. The dimensionality reductions
      selected in the --reduction parameter should be
      present in the loaded Seurat object.

  reduction:
    type: string
    inputBinding:
      prefix: "--reduction"
    doc: |
      Dimensionality reduction to be
      used for generating UMAP plots.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - library_id should
      correspond to all unique values from the new.ident
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied. Default: no extra metadata is added

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata by selected barcodes.
      First column should be named as barcode. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present. Default: all cells used, no extra
      metadata is added

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to group cells
      for optional subsetting when combined with --subset
      parameter. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Ignored if --subset is not set. Default: do not
      subset, include all cells into analysis.

  subset:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--subset"
    doc: |
      Values from the column set with --groupby parameter to
      subset cells before running differential expression
      analysis. Ignored if --groupby is not provided.
      Default: do not subset cells, include all of them

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells
      into two groups to run --second vs --first
      differential expression analysis. If --test parameter
      is set to deseq or deseq-lrt, the --splitby shouldn't
      put cells from the same dataset into the different
      comparison groups. May be one of the extra metadata
      columns added with --metadata or --barcodes
      parameters.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the first group of cells
      for differential expression analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the second group of
      cells for differential expression analysis.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcoxon"                   # (wilcox) Wilcoxon Rank Sum test
      - "likelihood-ratio"           # (bimod) Likelihood-ratio test
      - "t-test"                     # (t) Student's t-test
      - "negative-binomial"          # (negbinom) Negative Binomial Generalized Linear Model (supports --batchby)
      - "poisson"                    # (poisson) Poisson Generalized Linear Model (supports --batchby)
      - "logistic-regression"        # (LR) Logistic Regression (supports --batchby)
      - "mast"                       # (MAST) MAST package (supports --batchby)
      - "deseq"                      # DESeq2 Wald test on pseudobulk aggregated gene expression
      - "deseq-lrt"                  # DESeq2 LRT test on pseudobulk aggregated gene expression
    inputBinding:
      prefix: "--test"
    doc: |
      Test type to use in differential expression analysis.
      If set to deseq or deseq-lrt, gene expression will be
      aggregated to the pseudobulk level per dataset. For
      deseq, the pair-wise Wald test will be used. For
      deseq-lrt, the reduced formula will look like ~1 if
      --batchby parameter is omitted or will be set to
      ~batchby to exclude the criteria if interest (defined
      by --splitby). For all other values of the --test
      parameter the FindMarkers function will be used (genes
      will be prefiltered by minimum percentage >= 0.1 and
      by minimum log2FoldChange >= 0.25 before running
      differential expression analysis). Default: use
      FindMarkers with Wilcoxon Rank Sum test.

  batchby:
    type: string?
    inputBinding:
      prefix: "--batchby"
    doc: |
      Column from the Seurat object metadata to group cells
      into batches. If --test is set to deseq or deseq-lrt
      the --batchby parameter will be used in the design
      formula in the following way ~splitby+batchby.
      Additionally, the --batchby shouldn't put cells from
      the same dataset into the different batches. If --test
      is set to negative-binomial, poisson, logistic-
      regression, or mast it will be used as a latent
      variable in the FindMarkers function. Not supported
      for --test values equal to wilcoxon, likelihood-ratio,
      or t-test. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Default: do not model batch effect.

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially expressed genes with
      adjusted P-value not bigger than this value.
      Default: 0.05

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      Include only those genes that are detected in not lower than this
      fraction of cells in either of the two tested conditions.
      Default: 0.1

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      In the exploratory visualization part of
      the analysis output only differentially
      expressed genes with log2 Fold Change
      not smaller than this value.
      Default: 0.585 (1.5 folds)

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to label on the generated plots.
      Default: top 10 genes with the highest and the lowest
      log2FoldChange values.

  exclude_pattern:
    type: string?
    inputBinding:
      prefix: "--exclude"
    doc: |
      Regex pattern to identify and exclude specific genes
      from the differential expression analysis (not case-
      sensitive). If any of such genes are provided in the
      --genes parameter, they will be excluded from there as
      well. Default: use all genes

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
    inputBinding:
      prefix: "--cluster"
    doc: |
      Hopach clustering method to be run on the normalized
      read counts for the exploratory visualization part of
      the analysis. Clustering by column is supported only
      when --test is set to deseq or deseq-lrt. Default: do
      not run clustering

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
    inputBinding:
      prefix: "--rowdist"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if
      --cluster is set to column or not provided.
      Default: cosangle

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
    inputBinding:
      prefix: "--columndist"
    doc: |
      Distance metric for HOPACH column clustering. Ignored
      if --cluster is set to row or not provided.
      Default: euclid

  center_row:
    type: boolean?
    inputBinding:
      prefix: "--center"
    doc: |
      Apply mean centering for gene expression prior to
      running clustering by row. Ignored if --cluster is
      set to column or not provided. Default: do not
      center

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  color_theme:
    type:
    - "null"
    - type: enum
      symbols:
      - "gray"
      - "bw"
      - "linedraw"
      - "light"
      - "dark"
      - "minimal"
      - "classic"
      - "void"
    inputBinding:
      prefix: "--theme"
    doc: |
      Color theme for all generated plots. One of gray, bw,
      linedraw, light, dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_html_report:
    type: boolean?
    default: false
    doc: |
      Export tehcnical report. HTML format.
      Note, stdout will be less informative.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./sc

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between
      the workers when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  umap_spl_tst_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_tst.png"
    doc: |
      UMAP colored by selected for
      analysis cells. Split by tested
      condition; optionally subsetted
      to the specific group.
      PNG format.

  cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_cell_cnts.png"
    doc: |
      Number of cells per dataset or
      tested condition. Colored by tested
      condition; optionally subsetted to
      the specific group.
      PNG format.

  mds_plot_html:
    type: File?
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of pseudobulk aggregated
      normalized reads counts.
      HTML format.

  pca_1_2_plot_png:
    type: File?
    outputBinding:
      glob: "*_pca_1_2.png"
    doc: |
      Gene expression PCA (1,2).
      PNG format.

  pca_2_3_plot_png:
    type: File?
    outputBinding:
      glob: "*_pca_2_3.png"
    doc: |
      Gene expression PCA (2,3).
      PNG format

  dxpr_vlcn_plot_png:
    type: File?
    outputBinding:
      glob: "*_dxpr_vlcn.png"
    doc: |
      Differentially expressed genes.
      Volcano plot of differentially expressed genes.
      Highlighed genes are either provided by user or
      top 10 genes with the highest log2FoldChange
      values. The direction of comparison is defined
      as --second vs --first. Cells are optionally
      subsetted to the specific group and optionally
      coerced to the pseudobulk form.
      PNG format.

  xpr_avg_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_avg.png"
    doc: |
      Average gene expression plots split by dataset
      or tested condition for either user provided
      or top 10 differentially expressed genes with
      the highest log2FoldChange values.
      PNG format.

  xpr_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_dnst.png"
    doc: |
      Gene expression violin plots for either user
      provided or top 10 differentially expressed
      genes with the highest log2FoldChange values.
      The direction of comparison is defined as
      --second vs --first.
      PNG format.

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_*.png"
    doc: |
      UMAP colored by gene expression.
      Split by selected criteria; optionally
      subsetted to the specific group.
      PNG format.

  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.png"
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value and log2FC;
      optionally subsetted to the specific groups.
      PNG format.

  xpr_htmp_tsv:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.tsv"
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value and log2FC;
      optionally subsetted to the specific groups.
      TSV format.

  xpr_htmp_gct:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.gct"
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value and log2FC;
      optionally subsetted to the specific groups.
      GCT format.

  xpr_htmp_html:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.html"
    doc: |
      Gene expression heatmap.
      Filtered by adjusted p-value and log2FC;
      optionally subsetted to the specific groups.
      HTML format.

  diff_expr_genes:
    type: File?
    outputBinding:
      glob: "*_de_genes.tsv"
    doc: |
      Differentially expressed genes.
      Not filtered by adjusted p-value.
      TSV format.

  bulk_counts_gct:
    type: File?
    outputBinding:
      glob: "*_bulk_counts.gct"
    doc: |
      GSEA compatible not filtered normalized
      reads counts aggregated to pseudobulk
      form.
      GCT format.

  bulk_phntps_cls:
    type: File?
    outputBinding:
      glob: "*_bulk_phntps.cls"
    doc: |
      GSEA compatible phenotypes file defined
      based on --splitby, --first, and --second
      parameters.
      CLS format.

  all_plots_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.pdf"
    doc: |
      All generated plots.
      PDF format.

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

  human_log:
    type: File?
    outputBinding:
      glob: "error_report.txt"
    doc: |
      Human readable error log.
      TXT format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["Rscript"]
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_de_pseudobulk.R"]:"/usr/local/bin/sc_rna_de_pseudobulk.R")

stdout: sc_rna_de_pseudobulk_stdout.log
stderr: sc_rna_de_pseudobulk_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell RNA-Seq Differential Expression Analysis"
s:name: "Single-Cell RNA-Seq Differential Expression Analysis"
s:alternateName: "Single-Cell RNA-Seq Differential Expression Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-de-pseudobulk.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Single-Cell RNA-Seq Differential Expression Analysis

  Identifies differentially expressed genes between any
  two groups of cells, optionally aggregating gene
  expression data from single-cell to pseudobulk form.


s:about: |
  usage: /usr/local/bin/sc_rna_de_pseudobulk.R [-h] --query QUERY --reduction
                                              REDUCTION [--metadata METADATA]
                                              [--barcodes BARCODES]
                                              [--groupby GROUPBY]
                                              [--subset [SUBSET [SUBSET ...]]]
                                              --splitby SPLITBY --first FIRST
                                              --second SECOND
                                              [--test {wilcoxon,likelihood-ratio,t-test,negative-binomial,poisson,logistic-regression,mast,deseq,deseq-lrt}]
                                              [--batchby BATCHBY] [--padj PADJ]
                                              [--minpct MINPCT] [--logfc LOGFC]
                                              [--genes [GENES [GENES ...]]]
                                              [--exclude EXCLUDE]
                                              [--cluster {row,column,both}]
                                              [--rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
                                              [--columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
                                              [--center] [--pdf] [--verbose]
                                              [--output OUTPUT]
                                              [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                              [--cpus CPUS] [--memory MEMORY]
                                              [--seed SEED]

  Single-Cell RNA-Seq Differential Expression Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression information
                          stored in the RNA assay. The dimensionality reductions
                          selected in the --reduction parameter should be
                          present in the loaded Seurat object.
    --reduction REDUCTION
                          Dimensionality reduction to be used for generating
                          UMAP plots.
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - library_id should
                          correspond to all unique values from the new.ident
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          When combined with --barcodes parameter, first the
                          metadata will be extended, then barcode filtering will
                          be applied. Default: no extra metadata is added
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata by selected barcodes.
                          First column should be named as barcode. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --groupby GROUPBY     Column from the Seurat object metadata to group cells
                          for optional subsetting when combined with --subset
                          parameter. May be one of the extra metadata columns
                          added with --metadata or --barcodes parameters.
                          Ignored if --subset is not set. Default: do not
                          subset, include all cells into analysis.
    --subset [SUBSET [SUBSET ...]]
                          Values from the column set with --groupby parameter to
                          subset cells before running differential expression
                          analysis. Ignored if --groupby is not provided.
                          Default: do not subset cells, include all of them.
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first
                          differential expression analysis. If --test parameter
                          is set to deseq or deseq-lrt, the --splitby shouldn't
                          put cells from the same dataset into the different
                          comparison groups. May be one of the extra metadata
                          columns added with --metadata or --barcodes
                          parameters.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby parameter to define the first group of cells
                          for differential expression analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby parameter to define the second group of
                          cells for differential expression analysis.
    --test {wilcoxon,likelihood-ratio,t-test,negative-binomial,poisson,logistic-regression,mast,deseq,deseq-lrt}
                          Test type to use in differential expression analysis.
                          If set to deseq or deseq-lrt, gene expression will be
                          aggregated to the pseudobulk level per dataset. For
                          deseq, the pair-wise Wald test will be used. For
                          deseq-lrt, the reduced formula will look like ~1 if
                          --batchby parameter is omitted or will be set to
                          ~batchby to exclude the criteria if interest (defined
                          by --splitby). For all other values of the --test
                          parameter the FindMarkers function will be used (genes
                          will be prefiltered by minimum percentage >= 0.1 and
                          by minimum log2FoldChange >= 0.25 before running
                          differential expression analysis). Default: use
                          FindMarkers with Wilcoxon Rank Sum test.
    --batchby BATCHBY     Column from the Seurat object metadata to group cells
                          into batches. If --test is set to deseq or deseq-lrt
                          the --batchby parameter will be used in the design
                          formula in the following way ~splitby+batchby.
                          Additionally, the --batchby shouldn't put cells from
                          the same dataset into the different batches. If --test
                          is set to negative-binomial, poisson, logistic-
                          regression, or mast it will be used as a latent
                          variable in the FindMarkers function. Not supported
                          for --test values equal to wilcoxon, likelihood-ratio,
                          or t-test. May be one of the extra metadata columns
                          added with --metadata or --barcodes parameters.
                          Default: do not model batch effect.
    --padj PADJ           In the exploratory visualization part of the analysis
                          output only differentially expressed genes with
                          adjusted P-value not bigger than this value. Default:
                          0.05
    --minpct MINPCT       Include only those genes that are detected in not
                          lower than this fraction of cells in either of the two
                          tested conditions. Default: 0.1
    --logfc LOGFC         In the exploratory visualization part of the analysis
                          output only differentially expressed genes with log2
                          Fold Change not smaller than this value. Default:
                          0.585 (1.5 folds)
    --genes [GENES [GENES ...]]
                          Genes of interest to label on the generated plots.
                          Default: top 10 genes with the highest and the lowest
                          log2FoldChange values.
    --exclude EXCLUDE     Regex pattern to identify and exclude specific genes
                          from the differential expression analysis (not case-
                          sensitive). If any of such genes are provided in the
                          --genes parameter, they will be excluded from there as
                          well. Default: use all genes
    --cluster {row,column,both}
                          Hopach clustering method to be run on the normalized
                          read counts for the exploratory visualization part of
                          the analysis. Clustering by column is supported only
                          when --test is set to deseq or deseq-lrt. Default: do
                          not run clustering
    --rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH row clustering. Ignored if
                          --cluster is set to column or not provided. Default:
                          cosangle
    --columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH column clustering. Ignored
                          if --cluster is set to row or not provided. Default:
                          euclid
    --center              Apply mean centering for gene expression prior to
                          running clustering by row. Ignored if --cluster is set
                          to column or not provided. Default: do not centered
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42