cwlVersion: v1.2
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement


'sd:upstream':
  sc_filter_sample:
  - "sc-rna-filter.cwl"
  - "sc-multiome-filter.cwl"
  - "sc-atac-reduce.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell RNA-Seq or Multiome ATAC and RNA-Seq Filtering or ATAC-Seq Dimensionality Reduction Analysis"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay.
    'sd:upstreamSource': "sc_filter_sample/seurat_data_rds"
    'sd:localLabel': true

  barcodes_data:
    type: File?
    label: "Optional headerless TSV/CSV file with the list of barcodes to select cells of interest (one barcode per line)"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters loaded Seurat object
      to include only specific set of cells.
      Default: use all cells.

  cell_cycle_data:
    type: File?
    label: "Optional TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'"
    doc: |
      Path to the TSV/CSV file with the information for cell cycle score assignment.
      First column - 'phase', second column 'gene_id'. If loaded Seurat object already
      includes cell cycle scores in 'S.Score' and 'G2M.Score' metatada columns they will
      be removed.
      Default: skip cell cycle score assignment.

  dimensions:
    type: int?
    label: "Dimensionality to use in UMAP projection (from 1 to 50)"
    default: 40
    doc: |
      Dimensionality to use in UMAP projection (from 1 to 50). If single value N
      is provided, use from 1 to N PCs. If multiple values are provided, subset to
      only selected PCs.
      Default: from 1 to 10
    'sd:layout':
      advanced: true

  highly_var_genes_count:
    type: int?
    label: "Number of highly variable genes used in datasets integration, scaling and dimensionality reduction"
    default: 3000
    doc: |
      Number of highly variable genes used in datasets integration, scaling and
      dimensionality reduction.
      Default: 3000
    'sd:layout':
      advanced: true

  regress_mito_perc:
    type: boolean?
    label: "Regress the percentage of transcripts mapped to mitochondrial genes as a confounding source of variation"
    default: false
    doc: |
      Regress the percentage of transcripts mapped to mitochondrial genes as a
      confounding source of variation.
      Default: false
    'sd:layout':
      advanced: true

  regress_rna_umi:
    type: boolean?
    label: "Regress UMI per cell counts as a confounding source of variation"
    default: false
    doc: |
      Regress UMI per cell counts as a confounding source of variation.
      Default: false
    'sd:layout':
      advanced: true

  regress_genes:
    type: boolean?
    label: "Regress genes per cell counts as a confounding source of variation"
    default: false
    doc: |
      Regress genes per cell counts as a confounding source of variation.
      Default: false
    'sd:layout':
      advanced: true

  regress_cellcycle:
    type: boolean?
    label: "Regress cell cycle scores as a confounding source of variation"
    default: false
    doc: |
      Regress cell cycle scores as a confounding source of variation.
      Ignored if --cellcycle is not provided.
      Default: false
    'sd:layout':
      advanced: true

  umap_spread:
    type: float?
    label: "The effective scale of embedded points on UMAP (determines how clustered/clumped the embedded points are)"
    default: 1
    doc: |
      The effective scale of embedded points on UMAP. In combination with '--mindist'
      it determines how clustered/clumped the embedded points are.
      Default: 1
    'sd:layout':
      advanced: true

  umap_mindist:
    type: float?
    label: "Controls how tightly the embedding is allowed compress points together on UMAP"
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
    label: "Determines the number of neighboring points used in UMAP"
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
    - type: enum
      symbols:
      - "euclidean"
      - "cosine"
      - "correlation"
    label: "The metric to use to compute distances in high dimensional space for UMAP"
    default: "cosine"
    doc: |
      The metric to use to compute distances in high dimensional space for UMAP.
      Default: cosine
    'sd:layout':
      advanced: true

  umap_method:
    type:
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

  export_ucsc_cb:
    type: boolean?
    label: "Export results to UCSC Cell Browser"
    default: false
    doc: |
      Export results to UCSC Cell Browser.
      Default: false
    'sd:layout':
      advanced: true

  parallel_memory_limit:
    type:
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

  elbow_plot_png:
    type: File?
    outputSource: sc_rna_reduce/elbow_plot_png
    label: "Elbow plot (from cells PCA)"
    doc: |
      Elbow plot (from cells PCA).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Elbow plot (from cells PCA)'

  qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_rna_reduce/qc_dim_corr_plot_png
    label: "Correlation plots between QC metrics and cells PCA components"
    doc: |
      Correlation plots between QC metrics and cells PCA components.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Correlation plots between QC metrics and cells PCA components'

  umap_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_qc_mtrcs_plot_png
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
    outputSource: sc_rna_reduce/umap_plot_png
    label: "Cells UMAP"
    doc: |
      Cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Cells UMAP'

  umap_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_ph_plot_png
    label: "Split by cell cycle phase cells UMAP"
    doc: |
      Split by cell cycle phase cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by cell cycle phase cells UMAP'

  umap_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_mito_plot_png
    label: "Split by the percentage of transcripts mapped to mitochondrial genes cells UMAP"
    doc: |
      Split by the percentage of transcripts mapped to mitochondrial genes cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by the percentage of transcripts mapped to mitochondrial genes cells UMAP'

  umap_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_umi_plot_png
    label: "Split by the UMI per cell counts cells UMAP"
    doc: |
      Split by the UMI per cell counts cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by the UMI per cell counts cells UMAP'

  umap_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_gene_plot_png
    label: "Split by the genes per cell counts cells UMAP"
    doc: |
      Split by the genes per cell counts cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by the genes per cell counts cells UMAP'

  umap_spl_idnt_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_spl_idnt_plot_png
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
    outputSource: sc_rna_reduce/umap_spl_cnd_plot_png
    label: "Split by grouping condition cells UMAP"
    doc: |
      Split by grouping condition cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition cells UMAP'

  umap_gr_cnd_spl_ph_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_ph_plot_png
    label: "Grouped by condition split by cell cycle cells UMAP"
    doc: |
      Grouped by condition split by cell cycle cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by condition split by cell cycle cells UMAP'

  umap_gr_cnd_spl_mito_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_mito_plot_png
    label: "Grouped by condition split by the percentage of transcripts mapped to mitochondrial genes cells UMAP"
    doc: |
      Grouped by condition split by the percentage of transcripts mapped to mitochondrial genes cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by condition split by the percentage of transcripts mapped to mitochondrial genes cells UMAP'

  umap_gr_cnd_spl_umi_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_umi_plot_png
    label: "Grouped by condition split by the UMI per cell counts cells UMAP"
    doc: |
      Grouped by condition split by the UMI per cell counts cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by condition split by the UMI per cell counts cells UMAP'

  umap_gr_cnd_spl_gene_plot_png:
    type: File?
    outputSource: sc_rna_reduce/umap_gr_cnd_spl_gene_plot_png
    label: "Grouped by condition split by the genes per cell counts cells UMAP"
    doc: |
      Grouped by condition split by the genes per cell counts cells UMAP.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Grouped by condition split by the genes per cell counts cells UMAP'

  ucsc_cb_config_data:
    type: File?
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputSource: sc_rna_reduce/ucsc_cb_html_data
    label: "Directory with UCSC Cellbrowser html data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: sc_rna_reduce/ucsc_cb_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_rna_reduce/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  sc_rna_reduce_stdout_log:
    type: File
    outputSource: sc_rna_reduce/stdout_log
    label: "stdout log generated by sc_rna_reduce step"
    doc: |
      stdout log generated by sc_rna_reduce step

  sc_rna_reduce_stderr_log:
    type: File
    outputSource: sc_rna_reduce/stderr_log
    label: "stderr log generated by sc_rna_reduce step"
    doc: |
      stderr log generated by sc_rna_reduce step


steps:

  sc_rna_reduce:
    doc: |
      Integrates multiple single-cell RNA-Seq datasets,
      reduces dimensionality using PCA
    run: ../tools/sc-rna-reduce.cwl
    in:
      query_data_rds: query_data_rds
      barcodes_data: barcodes_data
      cell_cycle_data: cell_cycle_data
      normalization_method:
        default: "sctglm"
      integration_method:
        default: "seurat"
      highly_var_genes_count: highly_var_genes_count
      regress_mito_perc: regress_mito_perc
      regress_rna_umi: regress_rna_umi
      regress_genes: regress_genes
      regress_cellcycle: regress_cellcycle
      dimensions: dimensions
      umap_spread: umap_spread
      umap_mindist: umap_mindist
      umap_neighbors: umap_neighbors
      umap_metric: umap_metric
      umap_method: umap_method
      verbose:
        default: true
      export_ucsc_cb: export_ucsc_cb
      low_memory:
        default: true
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
    - elbow_plot_png
    - qc_dim_corr_plot_png
    - umap_qc_mtrcs_plot_png
    - umap_plot_png
    - umap_spl_ph_plot_png
    - umap_spl_mito_plot_png
    - umap_spl_umi_plot_png
    - umap_spl_gene_plot_png
    - umap_spl_idnt_plot_png
    - umap_spl_cnd_plot_png
    - umap_gr_cnd_spl_ph_plot_png
    - umap_gr_cnd_spl_mito_plot_png
    - umap_gr_cnd_spl_umi_plot_png
    - umap_gr_cnd_spl_gene_plot_png
    - ucsc_cb_config_data
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    when: $(inputs.folder_to_compress != null)
    in:
      folder_to_compress: sc_rna_reduce/ucsc_cb_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell RNA-Seq Dimensionality Reduction Analysis"
s:name: "Single-cell RNA-Seq Dimensionality Reduction Analysis"
s:alternateName: "Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-reduce.cwl
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
  Single-cell RNA-Seq Dimensionality Reduction Analysis

  Integrates multiple single-cell RNA-Seq datasets, reduces dimensionality using PCA.