cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.21


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the query RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R or
      run_seurat_wnn.R scripts.

  ref_data_rds:
    type: File?
    inputBinding:
      prefix: "--ref"
    doc: |
      Path to the reference RDS file for automatic cell type
      assignment. Either --ref or --ctype parameters should be
      provided.

  cell_type_data:
    type: File?
    inputBinding:
      prefix: "--ctype"
    doc: |
      Path to the cell types metadata TSV/CSV file for manual
      cell type assignment. Either --ref or --ctype parameters
      should be provided.

  ref_source_column:
    type: string?
    inputBinding:
      prefix: "--refsource"
    doc: |
      Metadata column from the reference RDS file to select
      labels to transfer. Required if running with --ref

  ref_reduction:
    type: string?
    inputBinding:
      prefix: "--refreduction"
    doc: |
      Reduction name from the reference RDS file to transfer
      embeddings from. Required if running with --ref

  query_source_column:
    type: string?
    inputBinding:
      prefix: "--querysource"
    doc: |
      Metadata column from the query RDS file to select clusters
      for manual cell type assignment. Required if running with
      --ctype. When running with --ref and --downsample parameter
      was provided, defines clusters each of which will be dowsampled
      to a fixed number of cells.

  query_target_column:
    type: string
    inputBinding:
      prefix: "--querytarget"
    doc: |
      Column to be added to metadata of the query RDS file for
      automatically or manually assigned cell types.

  high_var_features_count:
    type: int?
    inputBinding:
      prefix: "--highvarcount"
    doc: |
      Number of highly variable features to detect. Used for anchors
      detection between the reference and query datasets.
      Default: 3000

  regress_mito_perc:
    type: boolean?
    inputBinding:
      prefix: "--regressmt"
    doc: |
      Regress the percentage of transcripts mapped to mitochondrial genes
      as a confounding source of variation. Applied for query dataset only.
      Default: false

  dimensionality:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--ndim"
    doc: |
      Dimensions from the specified with --refreduction reference dataset's
      reduction that should be used in label transfer. If single value N is
      provided, use from 1 to N dimensions. If multiple values are provided,
      subset to only selected dimensions. Ignored when running with --ctype.
      Default: from 1 to 30

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--features"
    doc: |
      Features of interest to highlight.
      Default: None
  
  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./seurat

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  downsample:
    type: int?
    inputBinding:
      prefix: "--downsample"
    doc: |
      Downsample query dataset to include a fixed number of
      cells per cluster defined by --querysource parameter.
      Applied only when running automatic cell type assignment
      with both --ref and --querysource parameters provided.
      Default: do not downsample

  memory:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  ref_elbow_plot_png:
    type: File?
    outputBinding:
      glob: "*_ref_elbow.png"
    doc: |
      Elbow plot for selected reduction of the reference dataset.
      PNG format

  ref_elbow_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ref_elbow.pdf"
    doc: |
      Elbow plot for selected reduction of the reference dataset.
      PDF format

  ref_depth_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_ref_depth_corr.png"
    doc: |
      Correlation plot between depth and selected reduction
      dimensions of the reference dataset.
      PNG format

  ref_depth_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ref_depth_corr.pdf"
    doc: |
      Correlation plot between depth and selected reduction
      dimensions of the reference dataset.
      PDF format

  query_umap_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_query_umap_ctype.png"
    doc: |
      Query UMAP grouped by cell type.
      PNG format

  query_umap_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_query_umap_ctype.pdf"
    doc: |
      Query UMAP grouped by cell type.
      PDF format

  query_umap_ctype_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_query_umap_ctype_spl_by_cond.png"
    doc: |
      Query UMAP split by condition grouped by cell type.
      PNG format

  query_umap_ctype_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_query_umap_ctype_spl_by_cond.pdf"
    doc: |
      Query UMAP split by condition grouped by cell type.
      PDF format

  query_expr_avg_per_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_query_expr_avg_per_ctype.png"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of query dataset.
      PNG format

  query_expr_avg_per_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_query_expr_avg_per_ctype.pdf"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets
      PDF format

  query_expr_per_ctype_cell_plot_png:
    type: File?
    outputBinding:
      glob: "*_query_expr_per_ctype_cell.png"
    doc: |
      Log normalized gene expression per cell of query dataset with predicted cell types.
      PNG format

  query_expr_per_ctype_cell_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_query_expr_per_ctype_cell.pdf"
    doc: |
      Log normalized gene expression per cell of query dataset with predicted cell types.
      PDF format

  query_expr_dnst_per_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_query_expr_dnst_per_ctype.png"
    doc: |
      Log normalized gene expression densities per predicted cell type of query dataset.
      PNG format

  query_expr_dnst_per_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_query_expr_dnst_per_ctype.pdf"
    doc: |
      Log normalized gene expression densities per predicted cell type of query dataset.
      PDF format

  seurat_ctype_data_rds:
    type: File?
    outputBinding:
      glob: "*_ctype_data.rds"
    doc: |
      Query Seurat data with assigned cell types.
      RDS format

  cellbrowser_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["assign_cell_types.R"]

stdout: assign_cell_types_stdout.log
stderr: assign_cell_types_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Assign Cell Types"
s:name: "Single-cell Assign Cell Types"
s:alternateName: "Assigns cell types to Seurat clusters"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-assign-cell-types.cwl
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
  Single-cell Assign Cell Types
  ================================
  Assigns cell types to Seurat clusters


s:about: |
  usage: assign_cell_types.R [-h] --query QUERY
                                            (--ref REF | --ctype CTYPE)
                                            [--refsource REFSOURCE]
                                            [--refreduction REFREDUCTION]
                                            [--querysource QUERYSOURCE]
                                            --querytarget QUERYTARGET
                                            [--highvarcount HIGHVARCOUNT]
                                            [--regressmt]
                                            [--ndim [NDIM [NDIM ...]]]
                                            [--features [FEATURES [FEATURES ...]]]
                                            [--verbose] [--output OUTPUT]
                                            [--pdf] [--cpus CPUS]
                                            [--downsample DOWNSAMPLE]
                                            [--memory MEMORY]

  Assigns cell types to clusters

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the query RDS file to load Seurat object from.
                          RDS file can be produced by run_seurat.R or
                          run_seurat_wnn.R scripts.
    --ref REF             Path to the reference RDS file for automatic cell type
                          assignment. Either --ref or --ctype parameters should
                          be provided.
    --ctype CTYPE         Path to the cell types metadata TSV/CSV file for
                          manual cell type assignment. Either --ref or --ctype
                          parameters should be provided.
    --refsource REFSOURCE
                          Metadata column from the reference RDS file to select
                          labels to transfer. Required if running with --ref
    --refreduction REFREDUCTION
                          Reduction name from the reference RDS file to transfer
                          embeddings from. Required if running with --ref
    --querysource QUERYSOURCE
                          Metadata column from the query RDS file to select
                          clusters for manual cell type assignment. Required if
                          running with --ctype. When running with --ref and
                          --downsample parameter was provided, defines clusters
                          each of which will be dowsampled to a fixed number of
                          cells.
    --querytarget QUERYTARGET
                          Column to be added to metadata of the query RDS file
                          for automatically or manually assigned cell types.
    --highvarcount HIGHVARCOUNT
                          Number of highly variable features to detect. Used for
                          anchors detection between the reference and query
                          datasets. Default: 3000
    --regressmt           Regress the percentage of transcripts mapped to
                          mitochondrial genes as a confounding source of
                          variation. Applied for query dataset only. Default:
                          false
    --ndim [NDIM [NDIM ...]]
                          Dimensions from the specified with --refreduction
                          reference dataset's reduction that should be used in
                          label transfer. If single value N is provided, use
                          from 1 to N dimensions. If multiple values are
                          provided, subset to only selected dimensions. Ignored
                          when running with --ctype. Default: from 1 to 30
    --features [FEATURES [FEATURES ...]]
                          Features of interest to highlight. Default: None
    --verbose             Print debug information. Default: false
    --output OUTPUT       Output prefix for all generated files. Default:
                          ./seurat
    --pdf                 Export plots in PDF format in additions to PNG.
                          Default: false
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --downsample DOWNSAMPLE
                          Downsample query dataset to include a fixed number of
                          cells per cluster defined by --querysource parameter.
                          Applied only when running automatic cell type
                          assignment with both --ref and --querysource
                          parameters provided. Default: do not downsample
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32

