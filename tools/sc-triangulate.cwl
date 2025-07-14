cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.41
inputs:
  query_data_rds:
    type: File
    inputBinding:
      prefix: --query
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and/or chromatin accessibility information stored in the RNA
      and/or ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.
  barcodes_data:
    type: File?
    inputBinding:
      prefix: --barcodes
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata be selected barcodes. First column should be named as 'barcode'.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added
  query_source_column:
    type: string[]
    inputBinding:
      prefix: --source
    doc: |
      Columns from the metadata of the loaded Seurat object to select
      conflicting cells annotations.
  query_target_column:
    type: string?
    inputBinding:
      prefix: --target
    doc: |
      Suffix to be used as part of the columns names to save label
      integration result.
      Default: sctri
  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: --pdf
    doc: |
      Export plots in PDF.
      Default: false
  color_theme:
    type:
    - 'null'
    - type: enum
      symbols:
      - gray
      - bw
      - linedraw
      - light
      - dark
      - minimal
      - classic
      - void
    inputBinding:
      prefix: --theme
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic
  verbose:
    type: boolean?
    inputBinding:
      prefix: --verbose
    doc: |
      Print debug information.
      Default: false
  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: --h5seurat
    doc: |
      Save Seurat data to h5seurat file.
      Default: false
  export_h5ad_data:
    type: boolean?
    inputBinding:
      prefix: --h5ad
    doc: |
      Save raw counts from the RNA and/or ATAC assay(s) to h5ad file(s).
      Default: false
  export_loupe_data:
    type: boolean?
    inputBinding:
      prefix: --loupe
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: --cbbuild
    doc: |
      Export results to UCSC Cell Browser. Default: false
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
      prefix: --output
    doc: |
      Output prefix.
      Default: ./sc
  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: --memory
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
      prefix: --cpus
    doc: |
      Number of cores/cpus to use.
      Default: 1
  seed:
    type: int?
    inputBinding:
      prefix: --seed
    doc: |
      Seed number for random values.
      Default: 42
outputs:
  umap_tril_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_rnaumap.png'
    doc: |
      Cells UMAP with integrated labels (rnaumap dim. reduction).
      PNG format
  umap_tril_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_rnaumap.pdf'
    doc: |
      Cells UMAP with integrated labels (rnaumap dim. reduction).
      PDF format
  umap_tril_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_atacumap.png'
    doc: |
      Cells UMAP with integrated labels (atacumap dim. reduction).
      PNG format
  umap_tril_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_atacumap.pdf'
    doc: |
      Cells UMAP with integrated labels (atacumap dim. reduction).
      PDF format
  umap_tril_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_wnnumap.png'
    doc: |
      Cells UMAP with integrated labels (wnnumap dim. reduction).
      PNG format
  umap_tril_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tril_rd_wnnumap.pdf'
    doc: |
      Cells UMAP with integrated labels (wnnumap dim. reduction).
      PDF format
  umap_tria_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_rnaumap.png'
    doc: |
      Cells UMAP with winning annotations (rnaumap dim. reduction).
      PNG format
  umap_tria_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_rnaumap.pdf'
    doc: |
      Cells UMAP with winning annotations (rnaumap dim. reduction).
      PDF format
  umap_tria_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_atacumap.png'
    doc: |
      Cells UMAP with winning annotations (atacumap dim. reduction).
      PNG format
  umap_tria_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_atacumap.pdf'
    doc: |
      Cells UMAP with winning annotations (atacumap dim. reduction).
      PDF format
  umap_tria_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_wnnumap.png'
    doc: |
      Cells UMAP with winning annotations (wnnumap dim. reduction).
      PNG format
  umap_tria_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tria_rd_wnnumap.pdf'
    doc: |
      Cells UMAP with winning annotations (wnnumap dim. reduction).
      PDF format
  umap_tric_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_rnaumap.png'
    doc: |
      Cells UMAP with integration confidence scores (rnaumap dim. reduction).
      PNG format
  umap_tric_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_rnaumap.pdf'
    doc: |
      Cells UMAP with integration confidence scores (rnaumap dim. reduction).
      PDF format
  umap_tric_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_atacumap.png'
    doc: |
      Cells UMAP with integration confidence scores (atacumap dim. reduction).
      PNG format
  umap_tric_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_atacumap.pdf'
    doc: |
      Cells UMAP with integration confidence scores (atacumap dim. reduction).
      PDF format
  umap_tric_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_wnnumap.png'
    doc: |
      Cells UMAP with integration confidence scores (wnnumap dim. reduction).
      PNG format
  umap_tric_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_umap_tric_rd_wnnumap.pdf'
    doc: |
      Cells UMAP with integration confidence scores (wnnumap dim. reduction).
      PDF format
  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: '*_cellbrowser'
    doc: |
      UCSC Cell Browser configuration data.
  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: '*_cellbrowser/html_data'
    doc: |
      UCSC Cell Browser html data.
  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: '*_cellbrowser/html_data/index.html'
    doc: |
      UCSC Cell Browser html index.
  seurat_data_rds:
    type: File
    outputBinding:
      glob: '*_data.rds'
    doc: |
      Seurat object.
      RDS format
  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: '*_data.h5seurat'
    doc: |
      Seurat object.
      h5Seurat format
  seurat_rna_data_h5ad:
    type: File?
    outputBinding:
      glob: '*_rna_counts.h5ad'
    doc: |
      Seurat object.
      RNA counts.
      H5AD format.
  seurat_atac_data_h5ad:
    type: File?
    outputBinding:
      glob: '*_atac_counts.h5ad'
    doc: |
      Seurat object.
      ATAC counts.
      H5AD format.
  seurat_rna_data_cloupe:
    type: File?
    outputBinding:
      glob: '*_rna_counts.cloupe'
    doc: |
      Seurat object.
      RNA counts.
      Loupe format
  sc_report_html_file:
    type: File?
    outputBinding:
      glob: sc_report.html
    doc: |
      Tehcnical report.
      HTML format.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- Rscript
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_triangulate.R"]:"/usr/local/bin/sc_triangulate.R")
stdout: sc_triangulate_stdout.log
stderr: sc_triangulate_stderr.log
label: Single-Cell Label Integration Analysis
doc: |
  Single-Cell Label Integration Analysis

  Harmonizes conflicting annotations in single-cell genomics studies.
