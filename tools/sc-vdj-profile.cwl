cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.37


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This
      file should include gene expression information stored
      in the RNA assay, as well as pca and rnaumap dimensionality
      reductions applied to that assay. If loaded Seurat object
      includes multiple datasets, it should have a donor column
      to define grouping for clonotype calling.

  contigs_data:
    type: File
    inputBinding:
      prefix: "--contigs"
    doc: |
      Path to the file with high-level annotations of each
      high-confidence contig from cell-associated barcodes
      from the Cell Ranger Multi or Cell Ranger Aggregate
      experiments in TSV/CSV format.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - 'library_id' should correspond
      to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file
      columns are already present in the Seurat object metadata,
      they will be overwritten. When combined with --barcodes
      parameter, first the metadata will be extended, then barcode
      filtering will be applied. Default: no extra metadata is added

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones
      if those are present. Default: all cells used, no
      extra metadata is added

  cloneby:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene"
      - "nt"
      - "aa"
      - "strict"
    inputBinding:
      prefix: "--cloneby"
    doc: |
      Defines how to call the clonotype. gene: based on VDJC gene
      sequence. nt: based on the nucleotide sequence. aa: based on
      the amino acid sequence. strict: based on the combination of
      the nucleotide and gene sequences. Default: gene

  minimum_frequency:
    type: int?
    inputBinding:
      prefix: "--minfrequency"
    doc: |
      Minimum frequency (number of cells) per
      clonotype to be reported.
      Default: 3

  filterby:
    type:
    - "null"
    - type: enum
      symbols:
      - "cells"
      - "chains"
    inputBinding:
      prefix: "--filter"
    doc: |
      Stringency filters to be applied. 1) cells: remove
      cells with more than 2 chains. 2) chains: remove
      chains exceeding 2 (select the most expressed ones).
      Default: do not apply any filters.

  remove_partial:
    type: boolean?
    inputBinding:
      prefix: "--removepartial"
    doc: |
      Remove cells with only one chain detected.
      Default: keep all cells if at least one chain detected

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

  export_scope_data:
    type: boolean?
    inputBinding:
      prefix: "--scope"
    doc: |
      Save Seurat data to SCope compatible
      loom file. Only not normalized raw
      counts from the RNA assay will be
      saved. If loaded Seurat object doesn't
      have RNA assay this parameter will be
      ignored. Default: false

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

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  cl_qnt_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_qnt_gr_idnt_spl_ch.png"
    doc: |
      Percentage of unique clonotypes per dataset.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  cl_dnst_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_dnst_gr_idnt_spl_ch.png"
    doc: |
      Distribution of clonotype frequencies per dataset.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  allu_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_allu_gr_idnt_spl_ch.png"
    doc: |
      Proportion of top shared clonotypes between datasets.
      Split by chain; filtered by minimum clonotype
      frequency per donor; top clonotypes selected from
      each dataset.
      PNG format.

  hmst_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_gr_idnt_spl_ch.png"
    doc: |
      Proportion of clonotype frequencies per dataset.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  vrlp_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_gr_idnt_spl_ch.png"
    doc: |
      Overlap of clonotypes between datasets.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  dvrs_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_idnt_spl_ch.png"
    doc: |
      Diversity of clonotypes per dataset.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  gene_gr_idnt_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_gr_idnt_spl_ch.png"
    doc: |
      Distribution of gene usage per dataset.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  umap_cl_freq_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_cl_freq_spl_ch.png"
    doc: |
      UMAP colored by clonotype frequency.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  cl_qnt_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_qnt_gr_dnr_spl_ch.png"
    doc: |
      Percentage of unique clonotypes per donor.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  cl_dnst_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_dnst_gr_dnr_spl_ch.png"
    doc: |
      Distribution of clonotype frequencies per donor.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  allu_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_allu_gr_dnr_spl_ch.png"
    doc: |
      Proportion of top shared clonotypes between donors.
      Split by chain; filtered by minimum clonotype
      frequency per donor; top clonotypes selected from
      each donor.
      PNG format.

  hmst_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_gr_dnr_spl_ch.png"
    doc: |
      Proportion of clonotype frequencies per donor.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  vrlp_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_gr_dnr_spl_ch.png"
    doc: |
      Overlap of clonotypes between donors.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  dvrs_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_dnr_spl_ch.png"
    doc: |
      Diversity of clonotypes per donor.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  gene_gr_dnr_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_gr_dnr_spl_ch.png"
    doc: |
      Distribution of gene usage per donor.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  cl_qnt_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_qnt_gr_cnd_spl_ch.png"
    doc: |
      Percentage of unique clonotypes per
      grouping condition.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  cl_dnst_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_cl_dnst_gr_cnd_spl_ch.png"
    doc: |
      Distribution of clonotype frequencies
      per grouping condition.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  allu_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_allu_gr_cnd_spl_ch.png"
    doc: |
      Proportion of top shared clonotypes between
      grouping conditions.
      Split by chain; filtered by minimum clonotype
      frequency per donor; top clonotypes selected from
      each grouping condition.
      PNG format.

  hmst_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_gr_cnd_spl_ch.png"
    doc: |
      Proportion of clonotype frequencies per
      grouping condition.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  vrlp_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_gr_cnd_spl_ch.png"
    doc: |
      Overlap of clonotypes between grouping conditions.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
      PNG format.

  dvrs_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_cnd_spl_ch.png"
    doc: |
      Diversity of clonotypes per grouping condition.
      Split by chain; not filtered by clonotype frequency.
      PNG format.

  gene_gr_cnd_spl_ch_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_gr_cnd_spl_ch.png"
    doc: |
      Distribution of gene usage per grouping condition.
      Split by chain; filtered by minimum clonotype
      frequency per donor.
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

  clonotypes_data_tsv:
    type: File?
    outputBinding:
      glob: "*_clonotypes.tsv"
    doc: |
      Clonotypes data.
      Filtered by minimum clonotype
      frequency per donor.
      TSV format.

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

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_counts.h5ad"
    doc: |
      Seurat object.
      H5AD format.

  seurat_data_cloupe:
    type: File?
    outputBinding:
      glob: "*_counts.cloupe"
    doc: |
      Seurat object.
      Loupe format.

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_vdj_profile.R"]

stdout: sc_vdj_profile_stdout.log
stderr: sc_vdj_profile_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell Immune Profiling Analysis"
s:name: "Single-Cell Immune Profiling Analysis"
s:alternateName: "Single-Cell Immune Profiling Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-vdj-profile.cwl
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
  Single-Cell Immune Profiling Analysis

  Estimates clonotype diversity and dynamics from V(D)J
  sequencing data assembled into contigs.


s:about: |
  sc_vdj_profile.R [-h] --query QUERY --contigs CONTIGS
                                        [--metadata METADATA]
                                        [--barcodes BARCODES]
                                        [--cloneby {gene,nt,aa,strict}]
                                        [--minfrequency MINFREQUENCY]
                                        [--filter {cells,chains}]
                                        [--removepartial] [--pdf] [--verbose]
                                        [--h5seurat] [--h5ad] [--loupe]
                                        [--cbbuild] [--scope] [--output OUTPUT]
                                        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                        [--cpus CPUS] [--memory MEMORY]
                                        [--seed SEED]

  Single-Cell Immune Profiling Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include gene expression information stored
                          in the RNA assay, as well as pca and rnaumap
                          dimensionality reductions applied to that assay. If
                          loaded Seurat object includes multiple datasets, it
                          should have a donor column to define grouping for
                          clonotype calling.
    --contigs CONTIGS     Path to the file with high-level annotations of each
                          high-confidence contig from cell-associated barcodes
                          from the Cell Ranger Multi or Cell Ranger Aggregate
                          experiments in TSV/CSV format.
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          When combined with --barcodes parameter, first the
                          metadata will be extended, then barcode filtering will
                          be applied. Default: no extra metadata is added
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --cloneby {gene,nt,aa,strict}
                          Defines how to call the clonotype. gene: based on VDJC
                          gene sequence. nt: based on the nucleotide sequence.
                          aa: based on the amino acid sequence. strict: based on
                          the combination of the nucleotide and gene sequences.
                          Default: gene
    --minfrequency MINFREQUENCY
                          Minimum frequency (number of cells) per clonotype to
                          be reported. Default: 3
    --filter {cells,chains}
                          Stringency filters to be applied. cells: remove cells
                          with more than 2 chains. chains: remove chains
                          exceeding 2 (select the most expressed ones). Default:
                          do not apply any filters.
    --removepartial       Remove cells with only one chain detected. Default:
                          keep all cells if at least one chain detected
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
    --scope               Save Seurat data to SCope compatible loom file.
                          Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42