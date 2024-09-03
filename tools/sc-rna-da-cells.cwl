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
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay and selected with the --reduction
      parameter dimensionality reduction. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.

  reduction:
    type: string?
    inputBinding:
      prefix: "--reduction"
    doc: |
      Dimensionality reduction to be used for DA analysis.
      Default: pca

  dimensions:
    type: int?
    inputBinding:
      prefix: "--dimensions"
    doc: |
      Dimensionality to use when running DA analysis (from 1 to 50).
      Default: 10

  score_vector_knn:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--knn"
    doc: |
      Array of k values for kNN graph construction when calculating the
      score vector for each cell to represent the DA behavior in the
      neighborhood.
      Default: calculated based on the cells number

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
      Default: no extra metadata is added

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells into two groups
      to run --second vs --first DA analysis. May include columns from the
      extra metadata added with --metadata parameter.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the first group of cells for DA analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the second group of cells for DA analysis.

  resolution:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--resolution"
    doc: |
      Clustering resolution applied to DA cells to identify DA cells populations.
      Can be set as an array.
      Default: 0.01, 0.03, 0.05

  ranges:
    type:
    - "null"
    - float[]
    inputBinding:
      prefix: "--ranges"
    doc: |
      DA scores ranges for to filter out not significant cells.
      Default: calculated based on the permutation test

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
      Save raw counts from the RNA assay to h5ad file.
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

  da_perm_plot_png:
    type: File?
    outputBinding:
      glob: "*_da_perm.png"
    doc: |
      DA scores random permutations plot for second
      vs first biological conditions comparison.
      PNG format

  da_perm_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_da_perm.pdf"
    doc: |
      DA scores random permutations plot for second
      vs first biological conditions comparison.
      PDF format

  umap_rd_rnaumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_rnaumap_res_*.png"
    doc: |
      Clustered DA cells subpopulations UMAP (rnaumap dim. reduction).
      PNG format

  umap_rd_rnaumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_rnaumap_res_*.pdf"
    doc: |
      Clustered DA cells subpopulations UMAP (rnaumap dim. reduction).
      PDF format

  umap_rd_atacumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_atacumap_res_*.png"
    doc: |
      Clustered DA cells subpopulations UMAP (atacumap dim. reduction).
      PNG format

  umap_rd_atacumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_atacumap_res_*.pdf"
    doc: |
      Clustered DA cells subpopulations UMAP (atacumap dim. reduction).
      PDF format

  umap_rd_wnnumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_wnnumap_res_*.png"
    doc: |
      Clustered DA cells subpopulations UMAP (wnnumap dim. reduction).
      PNG format

  umap_rd_wnnumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_rd_wnnumap_res_*.pdf"
    doc: |
      Clustered DA cells subpopulations UMAP (wnnumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_rnaumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_rnaumap_res_*.png"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (rnaumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_rnaumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_rnaumap_res_*.pdf"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (rnaumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_atacumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_atacumap_res_*.png"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (atacumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_atacumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_atacumap_res_*.pdf"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (atacumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_wnnumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_wnnumap_res_*.png"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (wnnumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_wnnumap_res_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_umap_spl_cnd_rd_wnnumap_res_*.pdf"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (wnnumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_rnaumap_da_scr_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_rnaumap_da_scr.png"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (rnaumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_rnaumap_da_scr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_rnaumap_da_scr.pdf"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (rnaumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_atacumap_da_scr_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_atacumap_da_scr.png"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (atacumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_atacumap_da_scr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_atacumap_da_scr.pdf"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (atacumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_wnnumap_da_scr_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_wnnumap_da_scr.png"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (wnnumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_wnnumap_da_scr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_wnnumap_da_scr.pdf"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (wnnumap dim. reduction).
      PDF format

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
      RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat object.
      h5Seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_counts.h5ad"
    doc: |
      Seurat object.
      H5AD format

  seurat_data_cloupe:
    type: File?
    outputBinding:
      glob: "*_counts.cloupe"
    doc: |
      Seurat object.
      Loupe format

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

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
s:alternateName: "Detects cell subpopulations with differential abundance between datasets split by biological condition"

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

  Detects cell subpopulations with differential abundance
  between datasets split by biological condition.


s:about: |
  usage: /usr/local/bin/sc_rna_da_cells.R [-h] --query QUERY
                                          [--reduction REDUCTION]
                                          [--dimensions DIMENSIONS]
                                          [--knn [KNN [KNN ...]]]
                                          [--metadata METADATA] --splitby
                                          SPLITBY --first FIRST --second SECOND
                                          [--resolution [RESOLUTION [RESOLUTION ...]]]
                                          [--ranges RANGES RANGES] [--pdf]
                                          [--verbose] [--h5seurat] [--h5ad]
                                          [--cbbuild] [--output OUTPUT]
                                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                          [--cpus CPUS] [--memory MEMORY]
                                          [--seed SEED]

  Single-Cell Differential Abundance Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression information
                          stored in the RNA assay and selected with the
                          --reduction parameter dimensionality reduction.
                          Additionally, 'rnaumap', and/or 'atacumap', and/or
                          'wnnumap' dimensionality reductions should be present.
    --reduction REDUCTION
                          Dimensionality reduction to be used for DA analysis.
                          Default: pca
    --dimensions DIMENSIONS
                          Dimensionality to use when running DA analysis (from 1
                          to 50). Default: 10
    --knn [KNN [KNN ...]]
                          Array of k values for kNN graph construction when
                          calculating the score vector for each cell to
                          represent the DA behavior in the neighborhood.
                          Default: calculated based on the cells number
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          Default: no extra metadata is added
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first DA
                          analysis. May include columns from the extra metadata
                          added with --metadata parameter.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby to define the first group of cells for DA
                          analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby to define the second group of cells for DA
                          analysis.
    --resolution [RESOLUTION [RESOLUTION ...]]
                          Clustering resolution applied to DA cells to identify
                          DA cells populations. Can be set as an array. Default:
                          0.01, 0.03, 0.05
    --ranges RANGES RANGES
                          DA scores ranges for to filter out not significant
                          cells. Default: calculated based on the permutation
                          test
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA assay to h5ad file.
                          Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42