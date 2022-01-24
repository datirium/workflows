cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.17


inputs:

  seurat_data_rds:
    type: File
    inputBinding:
      prefix: "--rds"
    doc: |
      Path to the RDS file to load Seurat object from.
      RDS file produced by run_seurat.R script.

  cell_type_data:
    type: File
    inputBinding:
      prefix: "--ctype"
    doc: |
      Path to the cell types metadata TSV/CSV file with
      "cluster" and "type" columns

  source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column name to select clusters for cell type assignment

  target_column:
    type: string
    inputBinding:
      prefix: "--target"
    doc: |
      Column name to store assigned cell types

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--features"
    doc: |
      Features of interest to evaluate expression.
      Default: None
  
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

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
    doc: |
      Threads number
      Default: 1


outputs:

  umap_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_ctype.png"
    doc: |
      Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets.
      PNG format

  umap_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_ctype.pdf"
    doc: |
      Grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets.
      PDF format

  umap_ctype_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_ctype_spl_by_cond.png"
    doc: |
      Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets
      PNG format

  umap_ctype_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_ctype_spl_by_cond.pdf"
    doc: |
      Split by condition grouped by cell type UMAP projected PCA of filtered integrated/scaled datasets
      PDF format

  expr_avg_per_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_avg_per_ctype.png"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets
      PNG format

  expr_avg_per_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_avg_per_ctype.pdf"
    doc: |
      Scaled average log normalized gene expression per predicted cell type of filtered integrated/scaled datasets
      PDF format

  expr_per_ctype_cell_plot_png:
    type: File?
    outputBinding:
      glob: "*_per_ctype_cell.png"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types
      PNG format

  expr_per_ctype_cell_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_per_ctype_cell.pdf"
    doc: |
      Log normalized gene expression per cell of clustered filtered integrated/scaled datasets with predicted cell types
      PDF format

  expr_dnst_per_ctype_plot_png:
    type: File?
    outputBinding:
      glob: "*_dnst_per_ctype.png"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets
      PNG format

  expr_dnst_per_ctype_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dnst_per_ctype.pdf"
    doc: |
      Log normalized gene expression densities per predicted cell type of filtered integrated/scaled datasets
      PDF format

  seurat_ctype_data_rds:
    type: File?
    outputBinding:
      glob: "*_ctype_data.rds"
    doc: |
      Clustered filtered integrated/scaled Seurat data with assigned cell types.
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
  usage: assign_cell_types.R
        [-h] --rds RDS --ctype CTYPE --source SOURCE --target TARGET
        [--features [FEATURES [FEATURES ...]]] [--output OUTPUT] [--pdf]
        [--threads THREADS]

  Assigns cell types to clusters

  optional arguments:
    -h, --help            show this help message and exit
    --rds RDS             Path to the RDS file to load Seurat object from. RDS
                          file produced by run_seurat.R script
    --ctype CTYPE         Path to the cell types metadata TSV/CSV file with
                          cluster and type columns
    --source SOURCE       Column name to select clusters for cell type
                          assignment
    --target TARGET       Column name to store assigned cell types
    --features [FEATURES [FEATURES ...]]
                          Features of interest to highlight. Default: None
    --output OUTPUT       Output prefix. Default: ./seurat
    --pdf                 Export plots in PDF. Default: false
    --threads THREADS     Threads. Default: 1
