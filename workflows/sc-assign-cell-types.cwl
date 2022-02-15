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
    - var get_query_source_column = function(resolution, query_data_type) {
          if (query_data_type=="Seurat WNN Analysis") {
            return "wsnn_res."+resolution;
          } else if (query_data_type=="Seurat Cluster (single)") {
            return "RNA_snn_res."+resolution;
          } else if (query_data_type=="Seurat Cluster (integrated)") {
            return "integrated_snn_res."+resolution;
          }
      };
    - var get_query_target_column = function(resolution) {
          return "cluster_ext_type_res."+resolution;
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

  query_data_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "Seurat WNN Analysis"
      - "Seurat Cluster (single)"
      - "Seurat Cluster (integrated)"
    default: "Seurat Cluster (single)"
    label: "Seurat Cluster Experiment Type"
    doc: |
      If set to 'Seurat WNN Analysis', then 'query_source_column' will have prefix 'wsnn_res'.
      If set to 'Seurat Cluster (single)', then 'query_source_column' will have prefix 'RNA_snn_res'.
      If set to 'Seurat Cluster (integrated)', then 'query_source_column' will have prefix 'integrated_snn_res'.

  resolution:
    type: string
    label: "Clustering resolution to assign cell types to"
    doc: |
      Clustering resolution defines 'query_source_column' and 'query_target_column'
      inputs for 'assign_cell_types' step

  cell_type_data:
    type: File
    label: "TSV/CSV cell types metadata file with 'cluster' and 'type' columns"
    doc: |
      Path to the cell types metadata TSV/CSV file with
      "cluster" and "type" columns

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

  query_umap_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/query_umap_ctype_plot_png
    label: "Query UMAP grouped by cell type"
    doc: |
      Query UMAP grouped by cell type.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Cell Types'
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
        tab: 'Step 1. Cell Types'
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
        tab: 'Step 2. Gene expression'
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
        tab: 'Step 2. Gene expression'
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
        tab: 'Step 2. Gene expression'
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
      cell_type_data: cell_type_data
      query_source_column:
        source: [resolution, query_data_type]
        valueFrom: $(get_query_source_column(self[0], self[1]))
      query_target_column:
        source: resolution
        valueFrom: $(get_query_target_column(self))
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

label: "Single-cell Manual Cell Types Assignment"
s:name: "Single-cell Manual Cell Types Assignment"
s:alternateName: "Assigns cell types to Seurat clusters based on provided metadata file"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/sc-assign-cell-types.cwl
s:codeRepository: https://github.com/datirium/workflows
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
  Single-cell Manual Cell Types Assignment
  ========================================

  Assigns cell types to Seurat clusters based on provided metadata file.