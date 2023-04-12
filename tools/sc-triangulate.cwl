cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.18


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and/or chromatin accessibility information stored in the RNA
      and/or ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata be selected barcodes. First column should be named as 'barcode'.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added

  query_source_column:
    type: string[]
    inputBinding:
      prefix: "--source"
    doc: |
      Columns from the metadata of the loaded Seurat object to select
      conflicting cells annotations.

  query_target_column:
    type: string?
    inputBinding:
      prefix: "--target"
    doc: |
      Suffix to be used as part of the columns names to save label
      integration result.
      Default: sctri

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
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: "--h5seurat"
    doc: |
      Save Seurat data to h5seurat file.
      Default: false

  export_h5ad_data:
    type: boolean?
    inputBinding:
      prefix: "--h5ad"
    doc: |
      Save Seurat data to h5ad file.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
    doc: |
      Export results to UCSC Cell Browser. Default: false

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
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
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


outputs:

  umap_tril_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_rnaumap.png"
    doc: |
      Cells UMAP with integrated labels (rnaumap dim. reduction).
      PNG format

  umap_tril_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_rnaumap.pdf"
    doc: |
      Cells UMAP with integrated labels (rnaumap dim. reduction).
      PDF format

  umap_tril_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_atacumap.png"
    doc: |
      Cells UMAP with integrated labels (atacumap dim. reduction).
      PNG format

  umap_tril_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_atacumap.pdf"
    doc: |
      Cells UMAP with integrated labels (atacumap dim. reduction).
      PDF format

  umap_tril_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_wnnumap.png"
    doc: |
      Cells UMAP with integrated labels (wnnumap dim. reduction).
      PNG format

  umap_tril_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tril_rd_wnnumap.pdf"
    doc: |
      Cells UMAP with integrated labels (wnnumap dim. reduction).
      PDF format

  umap_tria_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_rnaumap.png"
    doc: |
      Cells UMAP with winning annotations (rnaumap dim. reduction).
      PNG format

  umap_tria_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_rnaumap.pdf"
    doc: |
      Cells UMAP with winning annotations (rnaumap dim. reduction).
      PDF format

  umap_tria_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_atacumap.png"
    doc: |
      Cells UMAP with winning annotations (atacumap dim. reduction).
      PNG format

  umap_tria_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_atacumap.pdf"
    doc: |
      Cells UMAP with winning annotations (atacumap dim. reduction).
      PDF format

  umap_tria_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_wnnumap.png"
    doc: |
      Cells UMAP with winning annotations (wnnumap dim. reduction).
      PNG format

  umap_tria_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tria_rd_wnnumap.pdf"
    doc: |
      Cells UMAP with winning annotations (wnnumap dim. reduction).
      PDF format

  umap_tric_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_rnaumap.png"
    doc: |
      Cells UMAP with integration confidence scores (rnaumap dim. reduction).
      PNG format

  umap_tric_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_rnaumap.pdf"
    doc: |
      Cells UMAP with integration confidence scores (rnaumap dim. reduction).
      PDF format

  umap_tric_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_atacumap.png"
    doc: |
      Cells UMAP with integration confidence scores (atacumap dim. reduction).
      PNG format

  umap_tric_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_atacumap.pdf"
    doc: |
      Cells UMAP with integration confidence scores (atacumap dim. reduction).
      PDF format

  umap_tric_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_wnnumap.png"
    doc: |
      Cells UMAP with integration confidence scores (wnnumap dim. reduction).
      PNG format

  umap_tric_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_tric_rd_wnnumap.pdf"
    doc: |
      Cells UMAP with integration confidence scores (wnnumap dim. reduction).
      PDF format

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat data in h5seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_data.h5ad"
    doc: |
      Seurat data in h5ad format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_triangulate.R"]

stdout: sc_triangulate_stdout.log
stderr: sc_triangulate_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Label Integration Analysis"
s:name: "Single-cell Label Integration Analysis"
s:alternateName: "Harmonizes conflicting annotations in single-cell genomics studies"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-triangulate.cwl
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
  Single-cell Label Integration Analysis

  Harmonizes conflicting annotations in single-cell genomics studies.


s:about: |
  usage: sc_triangulate.R
        [-h] --query QUERY [--barcodes BARCODES] --source SOURCE [SOURCE ...]
        [--target TARGET] [--pdf] [--verbose] [--h5seurat] [--h5ad] [--cbbuild]
        [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell Label Integration Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression and/or chromatin
                          accessibility information stored in the RNA and/or
                          ATAC assays correspondingly. Additionally, 'rnaumap',
                          and/or 'atacumap', and/or 'wnnumap' dimensionality
                          reductions should be present.
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --source SOURCE [SOURCE ...]
                          Columns from the metadata of the loaded Seurat object
                          to select conflicting cells annotations.
    --target TARGET       Suffix to be used as part of the columns names to save
                          label integration result. Default: sctri
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save Seurat data to h5ad file. Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32