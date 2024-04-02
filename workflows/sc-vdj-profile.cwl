cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement


"sd:upstream":
  sc_tools_sample:
  - "sc-rna-cluster.cwl"
  - "sc-ctype-assign.cwl"
  sc_vdj_sample:
  - "cellranger-multi.cwl"
  - "cellranger-aggr.cwl"


inputs:

  alias_:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with Clustered RNA-Seq Datasets"
    doc: |
      Analysis that includes single-cell
      RNA-Seq datasets run through either
      "Single-Cell Manual Cell Type
      Assignment" or "Single-Cell RNA-Seq
      Cluster Analysis" at any of the
      processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  contigs_data:
    type: File
    label: "Cell Ranger RNA+VDJ Sample"
    doc: |
      Any "Cell Ranger RNA+VDJ Sample" to
      load high level annotations of each
      high-confidence contig from the
      cell-associated barcodes. This sample
      can be obtained from either "Cell
      Ranger Count (RNA+VDJ)" or "Cell Ranger
      Aggregate (RNA, RNA+VDJ)" pipeline.
    "sd:upstreamSource": "sc_vdj_sample/filtered_contig_annotations_csv"
    "sd:localLabel": true

  query_source_column:
    type: string
    label: "Cells grouping"
    doc: |
      Single cell metadata column to group
      cells into clusters. Usually, in a form
      of "[rna|atac|wsnn]_res.X", where X is
      the clustering resolution. If cell types
      are available, add "custom_" prefix to
      the column name.

  cloneby:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene"
      - "nt"
      - "aa"
      - "strict"
    default: "gene"
    label: "Clonotype calling"
    doc: |
      Defines how to call the clonotype.
      gene: based on VDJC gene sequence.
      nt: based on the nucleotide sequence.
      aa: based on the amino acid sequence.
      strict: based on the combination of
      the nucleotide and gene sequences.
      Default: gene

  strictness:
    type:
    - "null"
    - type: enum
      symbols:
      - "removemulti"
      - "filtermulti"
      - "none"
    default: "none"
    label: "Stringency filter"
    doc: |
      Apply stringency filters. removemulti:
      remove any cell with more than 2 immune
      receptor chains. filtermulti: isolate
      the top 2 expressed chains in cell with
      multiple chains. none: do not apply any
      filters. Default: none

  datasets_metadata:
    type: File?
    label: "Datasets metadata (optional)"
    doc: |
      If the selected single-cell analysis
      includes multiple aggregated datasets,
      each of them can be assigned to a
      separate group by one or multiple
      categories. This can be achieved by
      providing a TSV/CSV file with
      "library_id" as the first column and
      any number of additional columns with
      unique names, representing the desired
      grouping categories.
# To obtain a proper                               this is not available yet, because we didn't refactor sc-rna-filter pipeline
# template of this file, download
# "datasets_metadata.tsv" output from the
# "Files" tab of the selected "Single-cell
# Analysis with Filtered RNA-Seq Datasets"
# and add extra columns as needed.

  barcodes_data:
    type: File?
    label: "Selected cell barcodes (optional)"
    doc: |
      A TSV/CSV file to optionally prefilter
      the single cell data by including only
      the cells with the selected barcodes.
      The provided file should include at
      least one column named "barcode", with
      one cell barcode per line. All other
      columns, except for "barcode", will be
      added to the single cell metadata loaded
      from "Single-cell Analysis with Clustered
      RNA-Seq Datasets" and can be utilized in
      the current or future steps of analysis.

  export_loupe_data:
    type: boolean?
    default: false
    label: "Save raw counts to Loupe file by accepting the EULA available at https://10xgen.com/EULA"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
    "sd:layout":
      advanced: true

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
    default: "classic"
    label: "Plots color theme"
    doc: |
      Color theme for all plots saved
      as PNG files.
      Default: classic
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "6"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 6
    "sd:layout":
      advanced: true


outputs:

  count_spl_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/count_spl_idnt_plot_png
    label: "Unique clonotypes, split by dataset"
    doc: |
      Unique clonotypes,
      split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Unique clonotypes, split by dataset"

  hmst_spl_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/hmst_spl_idnt_plot_png
    label: "Clonal space homeostasis, split by dataset"
    doc: |
      Clonal space homeostasis,
      split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Clonal space homeostasis, split by dataset"

  vrlp_spl_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/vrlp_spl_idnt_plot_png
    label: "Clonotypes similarity, split by dataset"
    doc: |
      Clonotypes similarity,
      split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Clonotypes similarity, split by dataset"

  ntwr_gr_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/ntwr_gr_idnt_plot_png
    label: "Clonotypes network, colored by dataset"
    doc: |
      Clonotypes network,
      colored by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Clonotypes network, colored by dataset"

  dvrs_gr_clst_spl_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/dvrs_gr_clst_spl_idnt_plot_png
    label: "Clonotypes diversity, colored by cluster, split by dataset"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Clonotypes diversity, colored by cluster, split by dataset"

  chrd_gr_idnt_plot_png:
    type: File?
    outputSource: vdj_profile/chrd_gr_idnt_plot_png
    label: "Shared clonotype, colored by dataset"
    doc: |
      Shared clonotype,
      colored by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Shared clonotype, colored by dataset"

  gene_spl_idnt_vdjc_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: vdj_profile/gene_spl_idnt_vdjc_plot_png
    label: "Relative usage of V, D, J, C genes, split by dataset"
    doc: |
      Relative usage of V, D, J, C
      genes, split by dataset
    "sd:visualPlugins":
    - image:
        tab: "Per dataset"
        Caption: "Relative usage of V, D, J, C genes, split by dataset"

  count_spl_clst_plot_png:
    type: File?
    outputSource: vdj_profile/count_spl_clst_plot_png
    label: "Unique clonotypes, split by cluster"
    doc: |
      Unique clonotypes,
      split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Unique clonotypes, split by cluster"

  hmst_spl_clst_plot_png:
    type: File?
    outputSource: vdj_profile/hmst_spl_clst_plot_png
    label: "Clonal space homeostasis, split by cluster"
    doc: |
      Clonal space homeostasis,
      split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Clonal space homeostasis, split by cluster"

  vrlp_spl_clst_plot_png:
    type: File?
    outputSource: vdj_profile/vrlp_spl_clst_plot_png
    label: "Clonotypes similarity, split by cluster"
    doc: |
      Clonotypes similarity,
      split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Clonotypes similarity, split by cluster"

  ntwr_gr_clst_plot_png:
    type: File?
    outputSource: vdj_profile/ntwr_gr_clst_plot_png
    label: "Clonotypes network, colored by cluster"
    doc: |
      Clonotypes network,
      colored by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Clonotypes network, colored by cluster"

  dvrs_gr_idnt_spl_clst_plot_png:
    type: File?
    outputSource: vdj_profile/dvrs_gr_idnt_spl_clst_plot_png
    label: "Clonotypes diversity, colored by dataset, split by cluster"
    doc: |
      Clonotypes diversity,
      colored by dataset,
      split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Clonotypes diversity, colored by dataset, split by cluster"

  chrd_gr_clst_plot_png:
    type: File?
    outputSource: vdj_profile/chrd_gr_clst_plot_png
    label: "Shared clonotype, colored by cluster"
    doc: |
      Shared clonotype,
      colored by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Shared clonotype, colored by cluster"

  gene_spl_clst_vdjc_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: vdj_profile/gene_spl_clst_vdjc_plot_png
    label: "Relative usage of V, D, J, C genes, split by cluster"
    doc: |
      Relative usage of V, D, J, C
      genes, split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per cluster"
        Caption: "Relative usage of V, D, J, C genes, split by cluster"

  count_spl_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/count_spl_cnd_plot_png
    label: "Unique clonotypes, split by grouping condition"
    doc: |
      Unique clonotypes,
      split by grouping
      condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Unique clonotypes, split by grouping condition"

  hmst_spl_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/hmst_spl_cnd_plot_png
    label: "Clonal space homeostasis, split by grouping condition"
    doc: |
      Clonal space homeostasis,
      split by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Clonal space homeostasis, split by grouping condition"

  vrlp_spl_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/vrlp_spl_cnd_plot_png
    label: "Clonotypes similarity, split by grouping condition"
    doc: |
      Clonotypes similarity,
      split by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Clonotypes similarity, split by grouping condition"

  ntwr_gr_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/ntwr_gr_cnd_plot_png
    label: "Clonotypes network, colored by grouping condition"
    doc: |
      Clonotypes network,
      colored by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Clonotypes network, colored by grouping condition"

  dvrs_gr_clst_spl_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/dvrs_gr_clst_spl_cnd_plot_png
    label: "Clonotypes diversity, colored by cluster, split by grouping condition"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by grouping condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Clonotypes diversity, colored by cluster, split by grouping condition"

  dvrs_gr_cnd_spl_clst_plot_png:
    type: File?
    outputSource: vdj_profile/dvrs_gr_cnd_spl_clst_plot_png
    label: "Clonotypes diversity, colored by grouping condition, split by cluster"
    doc: |
      Clonotypes diversity,
      colored by grouping condition,
      split by cluster
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Clonotypes diversity, colored by grouping condition, split by cluster"

  chrd_gr_cnd_plot_png:
    type: File?
    outputSource: vdj_profile/chrd_gr_cnd_plot_png
    label: "Shared clonotype, colored by grouping condition"
    doc: |
      Shared clonotype,
      colored by grouping
      condition
    "sd:visualPlugins":
    - image:
        tab: "Per group"
        Caption: "Shared clonotype, colored by grouping condition"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: vdj_profile/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: vdj_profile/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: vdj_profile/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_scope:
    type: File?
    outputSource: vdj_profile/seurat_data_scope
    label: "Seurat object in SCope compatible loom format"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  seurat_data_cloupe:
    type: File?
    outputSource: vdj_profile/seurat_data_cloupe
    label: "Seurat object in Loupe format"
    doc: |
      Seurat object.
      Loupe format.

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  vdj_profile_stdout_log:
    type: File
    outputSource: vdj_profile/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the vdj_profile step.

  vdj_profile_stderr_log:
    type: File
    outputSource: vdj_profile/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the vdj_profile step.


steps:

  vdj_profile:
    run: ../tools/sc-vdj-profile.cwl
    in:
      query_data_rds: query_data_rds
      contigs_data: contigs_data
      datasets_metadata: datasets_metadata
      barcodes_data: barcodes_data
      query_source_column: query_source_column
      cloneby: cloneby
      groupby:
        default: "new.ident"
      strictness:
        source: strictness
        valueFrom: $(self=="none"?null:self)
      color_theme: color_theme
      export_loupe_data: export_loupe_data
      export_pdf_plots:
        default: true
      verbose:
        default: true
      export_ucsc_cb:
        default: true
      export_scope_data:
        default: true
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 96
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - count_spl_idnt_plot_png
    - count_spl_idnt_plot_pdf
    - count_spl_clst_plot_png
    - count_spl_clst_plot_pdf
    - hmst_spl_idnt_plot_png
    - hmst_spl_idnt_plot_pdf
    - hmst_spl_clst_plot_png
    - hmst_spl_clst_plot_pdf
    - vrlp_spl_clst_plot_png
    - vrlp_spl_clst_plot_pdf
    - vrlp_spl_idnt_plot_png
    - vrlp_spl_idnt_plot_pdf
    - ntwr_gr_clst_plot_png
    - ntwr_gr_clst_plot_pdf
    - ntwr_gr_idnt_plot_png
    - ntwr_gr_idnt_plot_pdf
    - dvrs_gr_clst_spl_idnt_plot_png
    - dvrs_gr_clst_spl_idnt_plot_pdf
    - dvrs_gr_idnt_spl_clst_plot_png
    - dvrs_gr_idnt_spl_clst_plot_pdf
    - gene_spl_clst_vdjc_plot_png
    - gene_spl_clst_vdjc_plot_pdf
    - gene_spl_idnt_vdjc_plot_png
    - gene_spl_idnt_vdjc_plot_pdf
    - chrd_gr_clst_plot_png
    - chrd_gr_clst_plot_pdf
    - chrd_gr_idnt_plot_png
    - chrd_gr_idnt_plot_pdf
    - chrd_gr_cnd_plot_png
    - chrd_gr_cnd_plot_pdf
    - count_spl_cnd_plot_png
    - count_spl_cnd_plot_pdf
    - hmst_spl_cnd_plot_png
    - hmst_spl_cnd_plot_pdf
    - vrlp_spl_cnd_plot_png
    - vrlp_spl_cnd_plot_pdf
    - ntwr_gr_cnd_plot_png
    - ntwr_gr_cnd_plot_pdf
    - dvrs_gr_clst_spl_cnd_plot_png
    - dvrs_gr_clst_spl_cnd_plot_pdf
    - dvrs_gr_cnd_spl_clst_plot_png
    - dvrs_gr_cnd_spl_clst_plot_pdf
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_cloupe
    - seurat_data_scope
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - vdj_profile/count_spl_idnt_plot_pdf
        - vdj_profile/count_spl_clst_plot_pdf
        - vdj_profile/hmst_spl_idnt_plot_pdf
        - vdj_profile/hmst_spl_clst_plot_pdf
        - vdj_profile/vrlp_spl_clst_plot_pdf
        - vdj_profile/vrlp_spl_idnt_plot_pdf
        - vdj_profile/ntwr_gr_clst_plot_pdf
        - vdj_profile/ntwr_gr_idnt_plot_pdf
        - vdj_profile/dvrs_gr_clst_spl_idnt_plot_pdf
        - vdj_profile/dvrs_gr_idnt_spl_clst_plot_pdf
        - vdj_profile/gene_spl_clst_vdjc_plot_pdf
        - vdj_profile/gene_spl_idnt_vdjc_plot_pdf
        - vdj_profile/chrd_gr_clst_plot_pdf
        - vdj_profile/chrd_gr_idnt_plot_pdf
        - vdj_profile/chrd_gr_cnd_plot_pdf
        - vdj_profile/count_spl_cnd_plot_pdf
        - vdj_profile/hmst_spl_cnd_plot_pdf
        - vdj_profile/vrlp_spl_cnd_plot_pdf
        - vdj_profile/ntwr_gr_cnd_plot_pdf
        - vdj_profile/dvrs_gr_clst_spl_cnd_plot_pdf
        - vdj_profile/dvrs_gr_cnd_spl_clst_plot_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell Immune Profiling Analysis"
s:name: "Single-Cell Immune Profiling Analysis"
s:alternateName: "TCR/BCR clonotype dynamics analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-vdj-profile.cwl
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
  Single-Cell Immune Profiling Analysis

  Estimates clonotype diversity and dynamics from V(D)J
  sequencing data assembled into contigs