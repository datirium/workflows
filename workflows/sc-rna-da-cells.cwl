cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_numbers = function(line) {
          let splitted_line = line?line.split(/[\s,]+/).map(parseFloat):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  sc_tools_sample:
  - "sc-rna-reduce.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-ctype-assign.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through Single-cell RNA-Seq Dimensionality Reduction Analysis"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay and selected with the --reduction
      parameter dimensionality reduction. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.
    'sd:upstreamSource': "sc_tools_sample/seurat_data_rds"
    'sd:localLabel': true

  splitby:
    type: string
    label: "Column from the Seurat object metadata to split cells into two groups"
    doc: |
      Column from the Seurat object metadata to split cells into two groups
      to run --second vs --first DA analysis. May include columns from the
      extra metadata added with --metadata parameter.

  first_cond:
    type: string
    label: "Value from the Seurat object metadata column to define the first group of cells"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the first group of cells for DA analysis.

  second_cond:
    type: string
    label: "Value from the Seurat object metadata column to define the second group of cells"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define
      the second group of cells for DA analysis.

  dimensions:
    type: int?
    default: 20
    label: "Dimensionality to use when running DA analysis (from 1 to 50)"
    doc: |
      Dimensionality to use when running DA analysis (from 1 to 50).
      If single value N is provided, use from 1 to N PCs. If multiple
      values are provided, subset to only selected PCs.
      Default: from 1 to 10

  resolution:
    type: string?
    default: "0.05 0.1 0.15"
    label: "Clustering resolution applied to DA cells to identify DA cells populations"
    doc: |
      Clustering resolution applied to DA cells to identify DA cells populations.
      Can be set as an array.
      Default: 0.01, 0.03, 0.05

  ranges:
    type: string?
    default: "-0.5 0.5"
    label: " DA scores ranges for to filter out not significant cells"
    doc: |
      DA scores ranges for to filter out not significant cells.
      Default: calculated based on the permutation test

  datasets_metadata:
    type: File?
    label: "Path to the TSV/CSV file to optionally extend Seurat object metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten.
      Default: no extra metadata is added

  parallel_memory_limit:
    type:
    - "null"
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
    - "null"
    - type: enum
      symbols:
      - "64"
    default: "64"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Forced to 64 GB
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
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

  da_perm_plot_png:
    type: File?
    outputSource: da_cells/da_perm_plot_png
    label: "DA scores random permutations plot"
    doc: |
      DA scores random permutations plot for second
      vs first biological conditions comparison.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'DA scores random permutations plot'

  umap_rd_rnaumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_rd_rnaumap_res_plot_png
    label: "Clustered DA cells subpopulations RNA UMAP"
    doc: |
      Clustered DA cells subpopulations UMAP (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered DA cells subpopulations RNA UMAP'

  umap_rd_atacumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_rd_atacumap_res_plot_png
    label: "Clustered DA cells subpopulations ATAC UMAP"
    doc: |
      Clustered DA cells subpopulations UMAP (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered DA cells subpopulations ATAC UMAP'

  umap_rd_wnnumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_rd_wnnumap_res_plot_png
    label: "Clustered DA cells subpopulations WNN UMAP"
    doc: |
      Clustered DA cells subpopulations UMAP (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Overall'
        Caption: 'Clustered DA cells subpopulations WNN UMAP'

  umap_spl_cnd_rd_rnaumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_spl_cnd_rd_rnaumap_res_plot_png
    label: "Split by grouping condition clustered DA cells subpopulations RNA UMAP"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered DA cells subpopulations RNA UMAP'

  umap_spl_cnd_rd_atacumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_spl_cnd_rd_atacumap_res_plot_png
    label: "Split by grouping condition clustered DA cells subpopulations ATAC UMAP"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered DA cells subpopulations ATAC UMAP'

  umap_spl_cnd_rd_wnnumap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: da_cells/umap_spl_cnd_rd_wnnumap_res_plot_png
    label: "Split by grouping condition clustered DA cells subpopulations WNN UMAP"
    doc: |
      Split by grouping condition clustered DA cells subpopulations UMAP
      (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per group'
        Caption: 'Split by grouping condition clustered DA cells subpopulations WNN UMAP'

  umap_spl_idnt_rd_rnaumap_da_scr_plot_png:
    type: File?
    outputSource: da_cells/umap_spl_idnt_rd_rnaumap_da_scr_plot_png
    label: "Split by dataset cells RNA UMAP with DA scores"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (rnaumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset cells RNA UMAP with DA scores'

  umap_spl_idnt_rd_atacumap_da_scr_plot_png:
    type: File?
    outputSource: da_cells/umap_spl_idnt_rd_atacumap_da_scr_plot_png
    label: "Split by dataset cells ATAC UMAP with DA scores"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (atacumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset cells ATAC UMAP with DA scores'

  umap_spl_idnt_rd_wnnumap_da_scr_plot_png:
    type: File?
    outputSource: da_cells/umap_spl_idnt_rd_wnnumap_da_scr_plot_png
    label: "Split by dataset cells WNN UMAP with DA scores"
    doc: |
      Split by dataset cells UMAP with DA scores for second vs first
      biological conditions comparison (wnnumap dim. reduction).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Per dataset'
        Caption: 'Split by dataset cells WNN UMAP with DA scores'

  ucsc_cb_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory
    outputSource: da_cells/ucsc_cb_html_data
    label: "Directory with UCSC Cellbrowser html data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File
    outputSource: da_cells/ucsc_cb_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: da_cells/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  da_cells_stdout_log:
    type: File
    outputSource: da_cells/stdout_log
    label: "stdout log generated by da_cells step"
    doc: |
      stdout log generated by da_cells step

  da_cells_stderr_log:
    type: File
    outputSource: da_cells/stderr_log
    label: "stderr log generated by da_cells step"
    doc: |
      stderr log generated by da_cells step


steps:

  da_cells:
    run: ../tools/sc-rna-da-cells.cwl
    in:
      query_data_rds: query_data_rds
      datasets_metadata: datasets_metadata
      dimensions: dimensions
      splitby: splitby
      first_cond: first_cond
      second_cond: second_cond
      resolution:
        source: resolution
        valueFrom: $(split_numbers(self))
      ranges:
        source: ranges
        valueFrom: $(split_numbers(self))
      verbose:
        default: true
      export_ucsc_cb:
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
    - da_perm_plot_png
    - umap_rd_rnaumap_res_plot_png
    - umap_rd_atacumap_res_plot_png
    - umap_rd_wnnumap_res_plot_png
    - umap_spl_cnd_rd_rnaumap_res_plot_png
    - umap_spl_cnd_rd_atacumap_res_plot_png
    - umap_spl_cnd_rd_wnnumap_res_plot_png
    - umap_spl_idnt_rd_rnaumap_da_scr_plot_png
    - umap_spl_idnt_rd_atacumap_da_scr_plot_png
    - umap_spl_idnt_rd_wnnumap_da_scr_plot_png
    - ucsc_cb_config_data
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: da_cells/ucsc_cb_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Differential Abundance Analysis"
s:name: "Single-cell Differential Abundance Analysis"
s:alternateName: "Detects cell subpopulations with differential abundance between datasets split by biological condition"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-da-cells.cwl
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
  Single-cell Differential Abundance Analysis

  Detects cell subpopulations with differential abundance
  between datasets split by biological condition.