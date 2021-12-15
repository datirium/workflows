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
    - var get_source_column = function(resolution, from_aggregated) {
          if (from_aggregated) {
            return "integrated_snn_res."+resolution;
          } else {
            return "RNA_snn_res."+resolution;
          }
      };
    - var get_target_column = function(resolution) {
          return "cluster_ext_type_res."+resolution;
      };


'sd:upstream':
  seurat_cluster_sample:
  - "https://github.com/datirium/workflows/workflows/seurat-cluster.cwl"
  - "seurat-cluster.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  seurat_data_rds:
    type: File
    label: "Seurat Cluster Experiment"
    doc: |
      Path to the RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R script.
    'sd:upstreamSource': "seurat_cluster_sample/seurat_clst_data_rds"
    'sd:localLabel': true

  from_aggregated:
    type: boolean?
    default: true
    label: "Treat Seurat Cluster Experiment as aggregated"
    doc: |
      If set to true the 'source_column' and 'target_column' inputs will have
      prefix 'integrated_snn_res.{resolution}', otherwise 'RNA_res.{resolution}'

  resolution:
    type: string
    label: "Clustering resolution to assign cell types to"
    doc: |
      Clustering resolution define 'source_column' and 'target_column'
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

  threads:
    type: int?
    default: 2
    label: "Threads number to use"
    doc: |
      Threads number
    'sd:layout':
      advanced: true


outputs:

  umap_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/umap_ctype_plot_png
    label: "Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Cell Types'
        Caption: 'Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets'

  umap_ctype_plot_pdf:
    type: File?
    outputSource: assign_cell_types/umap_ctype_plot_pdf
    label: "Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  umap_ctype_spl_by_cond_plot_png:
    type: File?
    outputSource: assign_cell_types/umap_ctype_spl_by_cond_plot_png
    label: "Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Cell Types'
        Caption: 'Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets'

  umap_ctype_spl_by_cond_plot_pdf:
    type: File?
    outputSource: assign_cell_types/umap_ctype_spl_by_cond_plot_pdf
    label: "Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets"
    doc: |
      Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets
      PDF format

  expr_avg_per_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/expr_avg_per_ctype_plot_png
    label: "Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Cell Types'
        Caption: 'Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets'

  expr_avg_per_ctype_plot_pdf:
    type: File?
    outputSource: assign_cell_types/expr_avg_per_ctype_plot_pdf
    label: "Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets
      PDF format

  expr_per_ctype_cell_plot_png:
    type: File?
    outputSource: assign_cell_types/expr_per_ctype_cell_plot_png
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Cell Types'
        Caption: 'Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types'

  expr_per_ctype_cell_plot_pdf:
    type: File?
    outputSource: assign_cell_types/expr_per_ctype_cell_plot_pdf
    label: "Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types
      PDF format

  expr_dnst_per_ctype_plot_png:
    type: File?
    outputSource: assign_cell_types/expr_dnst_per_ctype_plot_png
    label: "Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Cell Types'
        Caption: 'Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets'

  expr_dnst_per_ctype_plot_pdf:
    type: File?
    outputSource: assign_cell_types/expr_dnst_per_ctype_plot_pdf
    label: "Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets
      PDF format

  seurat_ctype_data_rds:
    type: File
    outputSource: assign_cell_types/seurat_ctype_data_rds
    label: "Clustered filtered integrated/scaled Seurat data with assigned cell types"
    doc: |
      Clustered filtered integrated/scaled Seurat data with assigned cell types.
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
      seurat_data_rds: seurat_data_rds
      cell_type_data: cell_type_data
      source_column:
        source: [resolution, from_aggregated]
        valueFrom: $(get_source_column(self[0], self[1]))
      target_column:
        source: resolution
        valueFrom: $(get_target_column(self))
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      export_pdf_plots:
        default: true
      threads: threads
    out:
    - umap_ctype_plot_png
    - umap_ctype_plot_pdf
    - umap_ctype_spl_by_cond_plot_png
    - umap_ctype_spl_by_cond_plot_pdf
    - expr_avg_per_ctype_plot_png
    - expr_avg_per_ctype_plot_pdf
    - expr_per_ctype_cell_plot_png
    - expr_per_ctype_cell_plot_pdf
    - expr_dnst_per_ctype_plot_png
    - expr_dnst_per_ctype_plot_pdf
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

label: "Single-cell Assign Cell Types"
s:name: "Single-cell Assign Cell Types"
s:alternateName: "Assigns cell types to Seurat clusters"

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
  Single-cell Assign Cell Types
  =============================

  Assigns cell types to Seurat clusters.