cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.33


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This
      file should include gene expression information stored
      in the RNA assay, as well as 'pca' and 'rnaumap'
      dimensionality reductions applied to that assay.

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

  query_source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column from the metadata of the loaded Seurat
      object to select clusters from.

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

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the metadata of the loaded Seurat object
      to group cells for clonotype frequency calculation.
      Default: group by dataset

  strictness:
    type:
    - "null"
    - type: enum
      symbols:
      - "removemulti"
      - "filtermulti"
    inputBinding:
      prefix: "--strictness"
    doc: |
      Apply stringency filters. Removemulti: remove any cell
      with more than 2 immune receptor chains. Filtermulti:
      isolate the top 2 expressed chains in cell with multiple
      chains. Default: do not apply any filters

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


outputs:

  count_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_count_spl_idnt.png"
    doc: |
      Unique clonotypes,
      split by dataset
      PNG format

  count_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_count_spl_idnt.pdf"
    doc: |
      Unique clonotypes,
      split by dataset
      PDF format

  count_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_count_spl_clst.png"
    doc: |
      Unique clonotypes,
      split by cluster
      PNG format

  count_spl_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_count_spl_clst.pdf"
    doc: |
      Unique clonotypes,
      split by cluster
      PDF format

  hmst_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_idnt.png"
    doc: |
      Clonal space homeostasis,
      split by dataset
      PNG format

  hmst_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_idnt.pdf"
    doc: |
      Clonal space homeostasis,
      split by dataset
      PDF format

  hmst_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_clst.png"
    doc: |
      Clonal space homeostasis,
      split by cluster
      PNG format

  hmst_spl_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_clst.pdf"
    doc: |
      Clonal space homeostasis,
      split by cluster
      PDF format

  vrlp_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_clst.png"
    doc: |
      Clonotypes similarity,
      split by cluster
      PNG format

  vrlp_spl_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_clst.pdf"
    doc: |
      Clonotypes similarity,
      split by cluster
      PDF format

  vrlp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_idnt.png"
    doc: |
      Clonotypes similarity,
      split by dataset
      PNG format

  vrlp_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_idnt.pdf"
    doc: |
      Clonotypes similarity,
      split by dataset
      PDF format

  ntwr_gr_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_clst.png"
    doc: |
      Clonotypes network,
      colored by cluster
      PNG format

  ntwr_gr_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_clst.pdf"
    doc: |
      Clonotypes network,
      colored by cluster
      PDF format

  ntwr_gr_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_idnt.png"
    doc: |
      Clonotypes network,
      colored by dataset
      PNG format

  ntwr_gr_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_idnt.pdf"
    doc: |
      Clonotypes network,
      colored by dataset
      PDF format

  dvrs_gr_clst_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_clst_spl_idnt.png"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by dataset
      PNG format

  dvrs_gr_clst_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_clst_spl_idnt.pdf"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by dataset
      PDF format

  dvrs_gr_idnt_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_idnt_spl_clst.png"
    doc: |
      Clonotypes diversity,
      colored by dataset,
      split by cluster
      PNG format

  dvrs_gr_idnt_spl_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_idnt_spl_clst.pdf"
    doc: |
      Clonotypes diversity,
      colored by dataset,
      split by cluster
      PDF format

  gene_spl_clst_vdjc_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_gene_spl_clst_*.png"
    doc: |
      Relative usage of V, D, J, C
      genes, split by cluster
      PNG format

  gene_spl_clst_vdjc_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_gene_spl_clst_*.pdf"
    doc: |
      Relative usage of V, D, J, C
      genes, split by cluster
      PDF format

  gene_spl_idnt_vdjc_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_gene_spl_idnt_*.png"
    doc: |
      Relative usage of V, D, J, C
      genes, split by dataset
      PNG format

  gene_spl_idnt_vdjc_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_gene_spl_idnt_*.pdf"
    doc: |
      Relative usage of V, D, J, C
      genes, split by dataset
      PDF format

  chrd_gr_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_clst.png"
    doc: |
      Shared clonotype,
      colored by cluster
      PNG format

  chrd_gr_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_clst.pdf"
    doc: |
      Shared clonotype,
      colored by cluster
      PDF format

  chrd_gr_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_idnt.png"
    doc: |
      Shared clonotype,
      colored by dataset
      PNG format

  chrd_gr_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_idnt.pdf"
    doc: |
      Shared clonotype,
      colored by dataset
      PDF format

  chrd_gr_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_cnd.png"
    doc: |
      Shared clonotype,
      colored by grouping condition
      PNG format

  chrd_gr_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_chrd_gr_cnd.pdf"
    doc: |
      Shared clonotype,
      colored by grouping condition
      PDF format

  count_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_count_spl_cnd.png"
    doc: |
      Unique clonotypes,
      split by grouping condition
      PNG format

  count_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_count_spl_cnd.pdf"
    doc: |
      Unique clonotypes,
      split by grouping condition
      PDF format

  hmst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_cnd.png"
    doc: |
      Clonal space homeostasis,
      split by grouping condition
      PNG format

  hmst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_hmst_spl_cnd.pdf"
    doc: |
      Clonal space homeostasis,
      split by grouping condition
      PDF format

  vrlp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_cnd.png"
    doc: |
      Clonotypes similarity,
      split by grouping condition
      PNG format

  vrlp_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_vrlp_spl_cnd.pdf"
    doc: |
      Clonotypes similarity,
      split by grouping condition
      PDF format

  ntwr_gr_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_cnd.png"
    doc: |
      Clonotypes network,
      colored by grouping condition
      PNG format

  ntwr_gr_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_ntwr_gr_cnd.pdf"
    doc: |
      Clonotypes network,
      colored by grouping condition
      PDF format

  dvrs_gr_clst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_clst_spl_cnd.png"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by grouping condition
      PNG format

  dvrs_gr_clst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_clst_spl_cnd.pdf"
    doc: |
      Clonotypes diversity,
      colored by cluster,
      split by grouping condition
      PDF format

  dvrs_gr_cnd_spl_clst_plot_png:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_cnd_spl_clst.png"
    doc: |
      Clonotypes diversity,
      colored by grouping condition,
      split by cluster
      PNG format

  dvrs_gr_cnd_spl_clst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dvrs_gr_cnd_spl_clst.pdf"
    doc: |
      Clonotypes diversity,
      colored by grouping condition,
      split by cluster
      PDF format

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      Directory with UCSC Cellbrowser
      configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      Directory with UCSC Cellbrowser
      html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      HTML index file from the directory
      with UCSC Cellbrowser html data.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Reduced Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Reduced Seurat data in h5seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_data.h5ad"
    doc: |
      Reduced Seurat data in h5ad format

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Reduced Seurat data in SCope
      compatible loom format

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


label: "Single-cell Immune Profiling Analysis"
s:name: "Single-cell Immune Profiling Analysis"
s:alternateName: "TCR/BCR clonotype dynamics analysis"

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
  Single-cell Immune Profiling Analysis

  TCR/BCR clonotype dynamics analysis


s:about: |
  usage: sc_vdj_profile.R [-h] --query QUERY --contigs CONTIGS
                          [--metadata METADATA] [--barcodes BARCODES]
                          --source SOURCE
                          [--cloneby {gene,nt,aa,strict}] [--groupby GROUPBY]
                          [--strictness {removemulti,filtermulti}] [--pdf]
                          [--verbose] [--h5seurat] [--h5ad] [--cbbuild]
                          [--scope] [--output OUTPUT]
                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                          [--cpus CPUS] [--memory MEMORY]

  Single-cell Immune Profiling Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include gene expression information stored
                          in the RNA assay, as well as 'pca' and 'rnaumap'
                          dimensionality reductions applied to that assay.
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
    --source SOURCE       Column from the metadata of the loaded Seurat object
                          to select clusters from.
    --cloneby {gene,nt,aa,strict}
                          Defines how to call the clonotype. gene: based on VDJC
                          gene sequence. nt: based on the nucleotide sequence.
                          aa: based on the amino acid sequence. strict: based on
                          the combination of the nucleotide and gene sequences.
                          Default: gene
    --groupby GROUPBY     Column from the metadata of the loaded Seurat object
                          to group cells for clonotype frequency calculation.
                          Default: group by dataset
    --strictness {removemulti,filtermulti}
                          Apply stringency filters. Removemulti: remove any cell
                          with more than 2 immune receptor chains. Filtermulti:
                          isolate the top 2 expressed chains in cell with
                          multiple chains. Default: do not apply any filters.
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save Seurat data to h5ad file. Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --scope               Save Seurat data to SCope compatible loom file.
                          Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32