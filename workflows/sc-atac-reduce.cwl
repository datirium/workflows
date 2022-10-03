cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };

'sd:upstream':
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-multiome-filter.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through Single-cell Multiome ATAC and RNA-Seq Filtering Analysis"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      chromatin accessibility information stored in the ATAC assay.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  datasets_metadata:
    type: File?
    label: "Path to the TSV/CSV file to optionally extend Seurat object metadata with categorical values"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - 'library_id' should
      correspond to all unique values from the 'new.ident'
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied.
      Default: no extra metadata is added

  barcodes_data:
    type: File?
    label: "Optional TSV/CSV file to prefilter and extend metadata be barcodes. First column should be named as 'barcode'"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

  dimensions:
    type: int?
    label: "Dimensionality to use for datasets integration and UMAP projection (from 2 to 50)"
    default: 40
    doc: |
      Dimensionality to use for datasets integration and
      UMAP projection (from 2 to 50). If single value N is
      provided, use from 2 to N LSI components. If multiple
      values are provided, subset to only selected LSI
      components. In combination with --ntgr set to harmony,
      selected principle components will be used in Harmony
      integration.
      Default: from 2 to 10

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "log-tfidf"
      - "tf-logidf"
      - "logtf-logidf"
      - "idf"
    label: "TF-IDF normalization method applied to chromatin accessibility counts"
    default: "log-tfidf"
    doc: |
      TF-IDF normalization method applied to chromatin
      accessibility counts. log-tfidf - Stuart & Butler et
      al. 2019, tf-logidf - Cusanovich & Hill et al. 2018,
      logtf-logidf - Andrew Hill, idf - 10x Genomics,
      Default: log-tfidf
    'sd:layout':
      advanced: true

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "signac"
      - "harmony"
      - "none"
    label: "Integration method used for joint analysis of multiple datasets"
    default: "signac"
    doc: |
      Integration method used for joint analysis of multiple
      datasets. Automatically set to 'none' if loaded Suerat
      object includes only one dataset. Default: signac
    'sd:layout':
      advanced: true

  integrate_by:
    type: string?
    label: "Variable(s) to be integrated out when running multiple integration with Harmony"
    default: "new.ident"
    doc: |
      Column(s) from the Seurat object metadata to define
      the variable(s) that should be integrated out when
      running multiple datasets integration with harmony.
      May include columns from the extra metadata added with
      --metadata parameter. Ignored if --ntgr is not set to
      harmony.
      Default: new.ident
    'sd:layout':
      advanced: true

  minimum_var_peaks_perc:
    type: int?
    label: "Minimum percentile for identifying the top most common peaks as highly variable"
    default: 0
    doc: |
      Minimum percentile for identifying the top most common peaks as highly variable.
      For example, setting to 5 will use the the top 95 percent most common among all cells
      peaks as highly variable. These peaks are used for datasets integration, scaling
      and dimensionality reduction.
      Default: 0 (use all available peaks)
    'sd:layout':
      advanced: true

  umap_spread:
    type: float?
    label: "UMAP Spread - the effective scale of embedded points (determines how clustered/clumped the embedded points are)"
    default: 1
    doc: |
      The effective scale of embedded points on UMAP. In combination with '--mindist'
      it determines how clustered/clumped the embedded points are.
      Default: 1
    'sd:layout':
      advanced: true

  umap_mindist:
    type: float?
    label: "UMAP Min. Dist. - controls how tightly the embedding is allowed compress points together"
    default: 0.3
    doc: |
      Controls how tightly the embedding is allowed compress points together on UMAP.
      Larger values ensure embedded points are moreevenly distributed, while smaller
      values allow the algorithm to optimise more accurately with regard to local structure.
      Sensible values are in the range 0.001 to 0.5.
      Default:  0.3
    'sd:layout':
      advanced: true

  umap_neighbors:
    type: int?
    label: "UMAP Neighbors Number - determines the number of neighboring points used"
    default: 30
    doc: |
      Determines the number of neighboring points used in UMAP. Larger values will result
      in more global structure being preserved at the loss of detailed local structure.
      In general this parameter should often be in the range 5 to 50.
      Default: 30
    'sd:layout':
      advanced: true

  umap_metric:
    type:
    - "null"
    - type: enum
      symbols:
      - "euclidean"
      - "cosine"
      - "correlation"
    label: "UMAP Dist. Metric - the metric to use to compute distances in high dimensional space"
    default: "cosine"
    doc: |
      The metric to use to compute distances in high dimensional space for UMAP.
      Default: cosine
    'sd:layout':
      advanced: true

  umap_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "uwot"
      - "uwot-learn"
      - "umap-learn"
    label: "UMAP implementation to run (if set to 'umap-learn' use 'correlation' distance metric)"
    default: "uwot"
    doc: |
      UMAP implementation to run. If set to 'umap-learn' use --umetric 'correlation'
      Default: uwot
    'sd:layout':
      advanced: true

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
    default: "classic"
    label: "Color theme for all generated plots"
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic
    'sd:layout':
      advanced: true

  parallel_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "32"
    default: "32"
    label: "Maximum memory in GB allowed to be shared between the workers when using multiple CPUs"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Forced to 32 GB
    'sd:layout':
      advanced: true

  vector_memory_limit:
    type:
    - "null"
    - type: enum
      symbols:
      - "96"
    default: "96"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Forced to 96 GB
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
      Forced to 1
    'sd:layout':
      advanced: true


outputs:

  qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_atac_reduce/qc_dim_corr_plot_png
    label: "Correlation plots between QC metrics and cells LSI dimensions"
    doc: |
      Correlation plots between QC metrics and cells LSI dimensions.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Correlation plots between QC metrics and cells LSI dimensions'

  umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_qc_mtrcs_plot_png
    label: "QC metrics on cells UMAP"
    doc: |
      QC metrics on cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'QC metrics on cells UMAP'

  umap_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_plot_png
    label: "Cells UMAP"
    doc: |
      Cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells UMAP'

  umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_idnt_plot_png
    label: "Split by dataset cells UMAP"
    doc: |
      Split by dataset cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset cells UMAP'

  umap_spl_cnd_plot_png:
    type: File?
    outputSource: sc_atac_reduce/umap_spl_cnd_plot_png
    label: "Split by grouping condition cells UMAP"
    doc: |
      Split by grouping condition cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition cells UMAP'

  seurat_data_rds:
    type: File
    outputSource: sc_atac_reduce/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  sc_atac_reduce_stdout_log:
    type: File
    outputSource: sc_atac_reduce/stdout_log
    label: "stdout log generated by sc_atac_reduce step"
    doc: |
      stdout log generated by sc_atac_reduce step

  sc_atac_reduce_stderr_log:
    type: File
    outputSource: sc_atac_reduce/stderr_log
    label: "stderr log generated by sc_atac_reduce step"
    doc: |
      stderr log generated by sc_atac_reduce step


steps:

  sc_atac_reduce:
    doc: |
      Integrates multiple single-cell ATAC-Seq datasets,
      reduces dimensionality using LSI
    run: ../tools/sc-atac-reduce.cwl
    in:
      query_data_rds: query_data_rds
      barcodes_data: barcodes_data
      datasets_metadata: datasets_metadata
      normalization_method: normalization_method
      integration_method: integration_method
      integrate_by:
        source: integrate_by
        valueFrom: $(split_features(self))
      minimum_var_peaks_perc: minimum_var_peaks_perc
      dimensions: dimensions
      umap_spread: umap_spread
      umap_mindist: umap_mindist
      umap_neighbors: umap_neighbors
      umap_metric: umap_metric
      umap_method: umap_method
      verbose:
        default: true
      export_ucsc_cb:
        default: false
      color_theme: color_theme
      parallel_memory_limit:
        source: parallel_memory_limit
        valueFrom: $(parseInt(self))
      vector_memory_limit:
        source: vector_memory_limit
        valueFrom: $(parseInt(self))
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell ATAC-Seq Dimensionality Reduction Analysis"
s:name: "Single-cell ATAC-Seq Dimensionality Reduction Analysis"
s:alternateName: "Integrates multiple single-cell ATAC-Seq datasets, reduces dimensionality using LSI"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-reduce.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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
  Single-cell ATAC-Seq Dimensionality Reduction Analysis

  Integrates multiple single-cell ATAC-Seq datasets,
  reduces dimensionality using LSI.