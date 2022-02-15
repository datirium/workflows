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
          var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  seurat_cluster_sample:
  - "https://github.com/datirium/workflows/workflows/seurat-cluster.cwl"
  - "seurat-cluster.cwl"
  - "seurat-wnn-cluster.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Seurat Cluster Experiment"
    doc: |
      Path to the query RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R or
      run_seurat_wnn.R scripts.
    'sd:upstreamSource': "seurat_cluster_sample/seurat_clst_data_rds"
    'sd:localLabel': true

  ref_reduction:
    type: string?
    default: "pca"
    label: "Reduction name from the reference RDS file to transfer embeddings from"
    doc: |
      Reduction name from the reference RDS file to transfer
      embeddings from. Required if running with --ref

  dimensionality:
    type: int?
    default: 30
    label: "Number of dimensions to be transfered from the reference dataset"
    doc: |
      Dimensions from the specified with --refreduction reference dataset's
      reduction that should be used in label transfer. If single value N is
      provided, use from 1 to N dimensions. If multiple values are provided,
      subset to only selected dimensions. Ignored when running with --ctype.

  ref_source_column:
    type: string
    label: "Metadata column from the reference RDS file to select labels to transfer"
    doc: |
      Metadata column from the reference RDS file to select
      labels to transfer. Required if running with --ref

  ref_data_rds:
    type: File
    label: "Reference Seurat RDS file for automatic cell type assignment"
    doc: |
      Path to the reference RDS file for automatic cell type
      assignment. Either --ref or --ctype parameters should be
      provided.

  high_var_features_count:
    type: int?
    default: 3000
    label: "Number of highly variable features to use for anchors detection"
    doc: |
      Number of highly variable features to detect. Used for anchors
      detection between the reference and query datasets.
    'sd:layout':
      advanced: true

  regress_mito_perc:
    type: boolean?
    default: false
    label: "Regress the percentage of transcripts mapped to mitochondrial genes as a confounding source of variation"
    doc: |
      Regress the percentage of transcripts mapped to mitochondrial genes
      as a confounding source of variation. Applied for query dataset only.
    'sd:layout':
      advanced: true

  selected_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes of interest"
    doc: |
      Features of interest to evaluate expression.
    'sd:layout':
      advanced: true

  memory:
    type: int?
    default: 32
    label: "Maximum memory in GB allowed to be shared between the workers when using multiple CPUs"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 2
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
    'sd:layout':
      advanced: true


outputs:

  ref_elbow_plot_png:
    type: File?
    outputSource: assign_cell_types/ref_elbow_plot_png
    label: "Elbow plot for selected reduction of the reference dataset"
    doc: |
      Elbow plot for selected reduction of the reference dataset.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Reference QC'
        Caption: 'Elbow plot for selected reduction of the reference dataset'

  ref_depth_corr_plot_png:
    type: File?
    outputSource: assign_cell_types/ref_depth_corr_plot_png
    label: "Correlation plot between depth and selected reduction dimensions of the reference dataset"
    doc: |
      Correlation plot between depth and selected reduction
      dimensions of the reference dataset.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Reference QC'
        Caption: 'Correlation plot between depth and selected reduction dimensions of the reference dataset'

  query_umap_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/query_umap_ctype_plot_png
    label: "Query UMAP grouped by cell type"
    doc: |
      Query UMAP grouped by cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Cell Types'
        Caption: 'Query UMAP grouped by cell type'

  query_umap_ctype_spl_by_cond_plot_png:
    type: File?
    outputSource: assign_cell_types/query_umap_ctype_spl_by_cond_plot_png
    label: "Query UMAP split by condition grouped by cell type"
    doc: |
      Query UMAP split by condition grouped by cell type
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Cell Types'
        Caption: 'Query UMAP split by condition grouped by cell type'

  query_expr_avg_per_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/query_expr_avg_per_ctype_plot_png
    label: "Scaled average log normalized gene expression per predicted cell type of query dataset"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of query dataset
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Gene expression'
        Caption: 'Scaled average log normalized gene expression per predicted cell type of query dataset'

  query_expr_per_ctype_cell_plot_png:
    type: File?
    outputSource: assign_cell_types/query_expr_per_ctype_cell_plot_png
    label: "Log normalized gene expression per cell of query dataset with predicted cell types"
    doc: |
      Log normalized gene expression per cell of query dataset with predicted cell types
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Gene expression'
        Caption: 'Log normalized gene expression per cell of query dataset with predicted cell types'

  query_expr_dnst_per_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/query_expr_dnst_per_ctype_plot_png
    label: "Log normalized gene expression densities per predicted cell type of query dataset"
    doc: |
      Log normalized gene expression densities per predicted cell type of query dataset
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Gene expression'
        Caption: 'Log normalized gene expression densities per predicted cell type of query dataset'

  seurat_ctype_data_rds:
    type: File
    outputSource: assign_cell_types/seurat_ctype_data_rds
    label: "Query Seurat data with assigned cell types"
    doc: |
      Query Seurat data with assigned cell types.
      RDS format

  compressed_cellbrowser_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data

  cellbrowser_html_data:
    type: Directory
    outputSource: assign_cell_types/cellbrowser_html_data
    label: "Directory with UCSC Cellbrowser formatted html data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputSource: assign_cell_types/cellbrowser_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  assign_cell_types_stdout_log:
    type: File
    outputSource: assign_cell_types/stdout_log
    label: stdout log generated by 'assign_cell_types' step
    doc: |
      stdout log generated by 'assign_cell_types' step

  assign_cell_types_stderr_log:
    type: File
    outputSource: assign_cell_types/stderr_log
    label: stderr log generated by 'assign_cell_types' step
    doc: |
      stderr log generated by 'assign_cell_types' step


steps:

  assign_cell_types:
    run: ../tools/sc-assign-cell-types.cwl
    in:
      query_data_rds: query_data_rds
      ref_data_rds: ref_data_rds
      ref_source_column: ref_source_column
      ref_reduction: ref_reduction
      query_target_column:
        default: "transferred_cell_type"
      high_var_features_count: high_var_features_count
      query_source_column:
        default: "integrated_snn_res.0.4"
      downsample:
        default: 100
      regress_mito_perc: regress_mito_perc
      dimensionality: dimensionality
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      export_pdf_plots:
        default: false
      verbose:
        default: true
      memory: memory
      threads: threads
    out:
    - ref_elbow_plot_png
    - ref_depth_corr_plot_png
    - query_umap_ctype_plot_png
    - query_umap_ctype_spl_by_cond_plot_png
    - query_expr_avg_per_ctype_plot_png
    - query_expr_per_ctype_cell_plot_png
    - query_expr_dnst_per_ctype_plot_png
    - seurat_ctype_data_rds
    - cellbrowser_config_data
    - cellbrowser_html_data
    - cellbrowser_html_file
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: assign_cell_types/cellbrowser_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Reference Cell Types Assignment"
s:name: "Single-cell Reference Cell Types Assignment"
s:alternateName: "Assigns cell types to Seurat clusters based on provided reference Seureat RDS file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-reference-mapping.cwl
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
  Single-cell Reference Cell Types Assignment
  ===========================================

  Assigns cell types to Seurat clusters based on provided reference Seureat RDS file.