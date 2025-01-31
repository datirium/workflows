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
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and/or chromatin accessibility information stored in the RNA
      and ATAC assays correspondingly. Both dimensionality reductions selected in
      the --reduction and --embeddings parameters should be present in the loaded
      Seurat object.

  reduction:
    type: string
    inputBinding:
      prefix: "--reduction"
    doc: |
      Dimensionality reduction to be used for generating UMAP plots.

  embeddings:
    type: string?
    inputBinding:
      prefix: "--embeddings"
    doc: |
      Dimensionality reduction to extract embeddings
      for differential abundance analysis run with DAseq.
      Default: automatically selected based on the
      --reduction parameter.

  dimensions:
    type: int?
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to be used when running differential
      abundance analysis with DAseq (from 1 to 50).
      Default: 10

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten.
      Default: no extra metadata is added.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and extend Seurat object
      metadata by selected barcodes. First column should be named as barcode.
      If file includes any other columns they will be added to the Seurat object
      metadata ovewriting the existing ones if those are present.
      Default: all cells used, no extra metadata is added.

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells into two groups
      to run --second vs --first differential abundance analysis. May include
      columns from the extra metadata added with --metadata or --barcodes
      parameters.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the first group of cells for differential abundance analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the second group of cells for differential abundance analysis.

  groupby:
    type: string
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to group cells
      by categories, such as clusters, cell types, etc., when
      generating UMAP and composition plots.

  ranges:
    type:
    - "null"
    - float[]
    inputBinding:
      prefix: "--ranges"
    doc: |
      Minimum and maximum thresholds to filter out cells with
      the low (by absolute values) differential abundance scores.
      Default: calculated based on the permutation test.

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
      Save raw counts from the RNA and/or
      ATAC assay(s) to h5ad file(s).
      Default: false

  export_loupe_data:
    type: boolean?
    inputBinding:
      prefix: "--loupe"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false

  export_scope_data:
    type: boolean?
    inputBinding:
      prefix: "--scope"
    doc: |
      Save Seurat data to SCope compatible loom file. Only
      not normalized raw counts from the RNA assay will be
      saved. If loaded Seurat object doesn't have RNA assay
      this parameter will be ignored.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
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

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  cmp_bp_gr_tst_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_bp_gr_tst_spl_clst.png"
    doc: |
      Composition box plot colored by tested condition.
      Split by cluster; downsampled to the smallest
      dataset.
      PNG format.

  umap_gr_tst_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_tst.png"
    doc: |
      UMAP colored by tested condition. First downsampled
      to the smallest dataset, then downsampled to the
      smallest tested condition group.
      PNG format.

  umap_da_scr_ctg_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_da_scr_ctg.png"
    doc: |
      UMAP colored by differential abundance score,
      categorical scale. All cells; categories are
      defined based on the selected ranges for
      the differential abundance score.
      PNG format.

  umap_da_scr_cnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_da_scr_cnt.png"
    doc: |
      UMAP colored by differential abundance score,
      continuous scale. All cells.
      PNG format.

  umap_gr_clst_spl_tst_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_clst_spl_tst.png"
    doc: |
      UMAP colored by cluster. Split by tested condition;
      first downsampled to the smallest dataset, then
      downsampled to the smallest tested condition group.
      PNG format.

  cmp_gr_clst_spl_tst_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_clst_spl_tst.png"
    doc: |
      Composition plot colored by cluster. Split by tested
      condition; first downsampled to the smallest dataset,
      then downsampled to the smallest tested condition group.
      PNG format.

  cmp_gr_tst_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_tst_spl_clst.png"
    doc: |
      Composition plot colored by tested condition. Split by
      cluster; first downsampled to the smallest dataset,
      then downsampled to the smallest tested condition group.
      PNG format.

  rank_da_scr_plot_png:
    type: File?
    outputBinding:
      glob: "*_rank_da_scr.png"
    doc: |
      Estimated thresholds for
      differential abundance score.
      All cells.
      PNG format.

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

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      UCSC Cell Browser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      UCSC Cell Browser html index.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat object.
      h5Seurat format.

  seurat_rna_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_rna_counts.h5ad"
    doc: |
      Seurat object.
      RNA counts.
      H5AD format.

  seurat_atac_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_atac_counts.h5ad"
    doc: |
      Seurat object.
      ATAC counts.
      H5AD format.

  seurat_rna_data_cloupe:
    type: File?
    outputBinding:
      glob: "*_rna_counts.cloupe"
    doc: |
      Seurat object.
      RNA counts.
      Loupe format.

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

  human_log:
    type: File
    outputBinding:
      glob: "*_hlog.txt"
    doc: |
      Human readable error log.
      TXT format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["Rscript"]
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_da_cells.R"]:"/usr/local/bin/sc_rna_da_cells.R")

stdout: sc_rna_da_cells_stdout.log
stderr: sc_rna_da_cells_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell Differential Abundance Analysis"
s:name: "Single-Cell Differential Abundance Analysis"
s:alternateName: "Single-Cell Differential Abundance Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-da-cells.cwl
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
  Single-Cell Differential Abundance Analysis

  Compares the composition of cell types between
  two tested conditions


s:about: |
  usage: /usr/local/bin/sc_rna_da_cells.R [-h] --query QUERY --reduction
                                          REDUCTION [--embeddings EMBEDDINGS]
                                          [--dimensions DIMENSIONS]
                                          [--metadata METADATA]
                                          [--barcodes BARCODES] --splitby
                                          SPLITBY --first FIRST --second SECOND
                                          --groupby GROUPBY
                                          [--ranges RANGES RANGES] [--pdf]
                                          [--verbose] [--h5seurat] [--h5ad]
                                          [--loupe] [--cbbuild] [--scope]
                                          [--output OUTPUT]
                                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                          [--cpus CPUS] [--memory MEMORY]
                                          [--seed SEED]

  Single-Cell Differential Abundance Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression and/or chromatin
                          accessibility information stored in the RNA and ATAC
                          assays correspondingly. Both dimensionality reductions
                          selected in the --reduction and --embeddings
                          parameters should be present in the loaded Seurat
                          object.
    --reduction REDUCTION
                          Dimensionality reduction to be used for generating
                          UMAP plots.
    --embeddings EMBEDDINGS
                          Dimensionality reduction to extract embeddings for
                          differential abundance analysis run with DAseq.
                          Default: automatically selected based on the
                          --reduction parameter.
    --dimensions DIMENSIONS
                          Dimensionality to be used when running differential
                          abundance analysis with DAseq (from 1 to 50). Default:
                          10
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          Default: no extra metadata is added.
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata by selected barcodes.
                          First column should be named as barcode. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added.
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first
                          differential abundance analysis. May include columns
                          from the extra metadata added with --metadata or
                          --barcodes parameters.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby to define the first group of cells for
                          differential abundance analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby to define the second group of cells for
                          differential abundance analysis.
    --groupby GROUPBY     Column from the Seurat object metadata to group cells
                          by categories, such as clusters, cell types, etc.,
                          when generating UMAP and composition plots.
    --ranges RANGES RANGES
                          Minimum and maximum thresholds to filter out cells
                          with the low (by absolute values) differential
                          abundance scores. Default: calculated based on the
                          permutation test.
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA and/or ATAC assay(s) to
                          h5ad file(s). Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --scope               Save Seurat data to SCope compatible loom file. Only
                          not normalized raw counts from the RNA assay will be
                          saved. If loaded Seurat object doesn't have RNA assay
                          this parameter will be ignored. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42