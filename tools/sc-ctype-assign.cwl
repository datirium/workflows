cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.21


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include
      genes expression and/or chromatin accessibility information stored in the RNA
      and ATAC assays correspondingly. Additionally, 'rnaumap', and/or 'atacumap',
      and/or 'wnnumap' dimensionality reductions should be present.

  cell_type_data:
    type: File
    inputBinding:
      prefix: "--celltypes"
    doc: |
      Path to the TSV/CSV file for manual cell type assignment for each of the clusters.
      First column - 'cluster', second column may have arbitrary name.

  query_source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column from the metadata of the loaded Seurat object to select clusters from.

  query_target_column:
    type: string
    inputBinding:
      prefix: "--target"
    doc: |
      Column from the metadata of the loaded Seurat object to save manually
      assigned cell types. Should start with 'custom_', otherwise, it won't
      be shown in UCSC Cell Browser.

  identify_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--diffgenes"
    doc: |
      Identify differentially expressed genes (putative gene markers) for
      assigned cell types. Ignored if loaded Seurat object doesn't include
      genes expression information stored in the RNA assay.
      Default: false

  identify_diff_peaks:
    type: boolean?
    inputBinding:
      prefix: "--diffpeaks"
    doc: |
      Identify differentially accessible peaks for assigned cell types. Ignored
      if loaded Seurat object doesn't include chromatin accessibility information
      stored in the ATAC assay.
      Default: false

  rna_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--rnalogfc"
    doc: |
      For putative gene markers identification include only those genes that
      on average have log fold change difference in expression between every
      tested pair of cell types not lower than this value. Ignored if '--diffgenes'
      is not set or RNA assay is not present.
      Default: 0.25

  rna_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--rnaminpct"
    doc: |
      For putative gene markers identification include only those genes that
      are detected in not lower than this fraction of cells in either of the
      two tested cell types. Ignored if '--diffgenes' is not set or RNA assay
      is not present.
      Default: 0.1

  only_positive_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--rnaonlypos"
    doc: |
      For putative gene markers identification return only positive markers.
      Ignored if '--diffgenes' is not set or RNA assay is not present.
      Default: false

  rna_test_to_use:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcox"
      - "bimod"
      - "roc"
      - "t"
      - "negbinom"
      - "poisson"
      - "LR"
      - "MAST"
      - "DESeq2"
    inputBinding:
      prefix: "--rnatestuse"
    doc: |
      Statistical test to use for putative gene markers identification.
      Ignored if '--diffgenes' is not set or RNA assay is not present.
      Default: wilcox

  atac_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--ataclogfc"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      on average have log fold change difference in the chromatin accessibility between
      every tested pair of cell types not lower than this value. Ignored if '--diffpeaks'
      is not set or ATAC assay is not present.
      Default: 0.25

  atac_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--atacminpct"
    doc: |
      For differentially accessible peaks identification include only those peaks that
      are detected in not lower than this fraction of cells in either of the two tested
      cell types. Ignored if '--diffpeaks' is not set or ATAC assay is not present.
      Default: 0.05

  atac_test_to_use:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcox"
      - "bimod"
      - "roc"
      - "t"
      - "negbinom"
      - "poisson"
      - "LR"
      - "MAST"
      - "DESeq2"
    inputBinding:
      prefix: "--atactestuse"
    doc: |
      Statistical test to use for differentially accessible peaks identification.
      Ignored if '--diffpeaks' is not set or ATAC assay is not present.
      Default: LR

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment used in the loaded Seurat
      object. File should be saved in TSV format with tbi-index file. Ignored if the
      loaded Seurat object doesn't include ATAC assay.

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to build gene expression and/or Tn5 insertion frequency plots
      for the nearest peaks. To build gene expression plots the loaded Seurat object
      should include RNA assay. To build Tn5 insertion frequency plots for the nearest
      peaks the loaded Seurat object should include ATAC assay as well as the --fragments
      file should be provided.
      Default: None

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

  umap_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.png"
    doc: |
      Cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  umap_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.pdf"
    doc: |
      Cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  umap_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.png"
    doc: |
      Cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  umap_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.pdf"
    doc: |
      Cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  umap_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.png"
    doc: |
      Cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  umap_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.pdf"
    doc: |
      Cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_rnaumap.png"
    doc: |
      Split by dataset cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_rnaumap.pdf"
    doc: |
      Split by dataset cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_atacumap.png"
    doc: |
      Split by dataset cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_atacumap.pdf"
    doc: |
      Split by dataset cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  umap_spl_idnt_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_wnnumap.png"
    doc: |
      Split by dataset cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  umap_spl_idnt_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_idnt_rd_wnnumap.pdf"
    doc: |
      Split by dataset cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_rnaumap.png"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_rnaumap.pdf"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_atacumap.png"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_atacumap.pdf"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  umap_spl_cnd_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_wnnumap.png"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  umap_spl_cnd_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_cnd_rd_wnnumap.pdf"
    doc: |
      Split by grouping condition cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  umap_spl_ph_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_rnaumap.png"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  umap_spl_ph_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_rnaumap.pdf"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  umap_spl_ph_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_atacumap.png"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  umap_spl_ph_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_atacumap.pdf"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  umap_spl_ph_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_wnnumap.png"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  umap_spl_ph_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_spl_ph_rd_wnnumap.pdf"
    doc: |
      Split by cell cycle phase cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_idnt.png"
    doc: |
      Grouped by cell type split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_ctyp_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_idnt.pdf"
    doc: |
      Grouped by cell type split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_idnt_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_ctyp.png"
    doc: |
      Grouped by dataset split by cell type cells composition plot. Downsampled.
      PNG format

  cmp_gr_idnt_spl_ctyp_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_ctyp.pdf"
    doc: |
      Grouped by dataset split by cell type cells composition plot. Downsampled.
      PDF format

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_idnt_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.pdf"
    doc: |
      Grouped by cell cycle phase split by dataset cells composition plot. Downsampled.
      PDF format

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_cnd.png"
    doc: |
      Grouped by cell type split by condition cells composition plot. Downsampled.
      PNG format

  cmp_gr_ctyp_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_cnd.pdf"
    doc: |
      Grouped by cell type split by condition cells composition plot. Downsampled.
      PDF format

  cmp_gr_cnd_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_ctyp.png"
    doc: |
      Grouped by condition split by cell type cells composition plot. Downsampled.
      PNG format

  cmp_gr_cnd_spl_ctyp_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_ctyp.pdf"
    doc: |
      Grouped by condition split by cell type cells composition plot. Downsampled.
      PDF format

  cmp_gr_ph_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_ctyp.png"
    doc: |
      Grouped by cell cycle phase split by cell type cells composition plot. Downsampled.
      PNG format

  cmp_gr_ph_spl_ctyp_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_ctyp.pdf"
    doc: |
      Grouped by cell cycle phase split by cell type cells composition plot. Downsampled.
      PDF format

  xpr_avg_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_avg.png"
    doc: |
      Log normalized scaled average gene expression per cell type.
      PNG format

  xpr_avg_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_xpr_avg.pdf"
    doc: |
      Log normalized scaled average gene expression per cell type.
      PDF format

  xpr_dnst_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_*.png"
    doc: |
      Log normalized gene expression density per cell type.
      PNG format

  xpr_dnst_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_*.pdf"
    doc: |
      Log normalized gene expression density per cell type.
      PDF format

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_rnaumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  xpr_per_cell_rd_rnaumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_rnaumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_atacumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  xpr_per_cell_rd_atacumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_atacumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_wnnumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  xpr_per_cell_rd_wnnumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_wnnumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  xpr_per_cell_sgnl_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_rnaumap_*.png"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PNG format

  xpr_per_cell_sgnl_rd_rnaumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_rnaumap_*.pdf"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (rnaumap dim. reduction).
      PDF format

  xpr_per_cell_sgnl_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_atacumap_*.png"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (atacumap dim. reduction).
      PNG format

  xpr_per_cell_sgnl_rd_atacumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_atacumap_*.pdf"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (atacumap dim. reduction).
      PDF format

  xpr_per_cell_sgnl_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_wnnumap_*.png"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PNG format

  xpr_per_cell_sgnl_rd_wnnumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_rd_wnnumap_*.pdf"
    doc: |
      Log normalized gene expression density on cells UMAP with assigned cell types (wnnumap dim. reduction).
      PDF format

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_*.png"
    doc: |
      Tn5 insertion frequency plot around gene.
      PNG format

  cvrg_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_*.pdf"
    doc: |
      Tn5 insertion frequency plot around gene.
      PDF format

  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.png"
    doc: |
      Normalized gene expression heatmap grouped by cell type.
      PNG format

  xpr_htmp_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.pdf"
    doc: |
      Normalized gene expression heatmap grouped by cell type.
      PDF format

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Differentially expressed genes between each pair of cell types.
      TSV format

  peak_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_peak_markers.tsv"
    doc: |
      Differentially accessible peaks between each pair of cell types.
      TSV format

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

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_ctype_assign.R"]

stdout: sc_ctype_assign_stdout.log
stderr: sc_ctype_assign_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Manual Cell Type Assignment"
s:name: "Single-cell Manual Cell Type Assignment"
s:alternateName: "Assigns cell types for clusters based on the provided metadata file"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-ctype-assign.cwl
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
  Single-cell Manual Cell Type Assignment

  Assigns cell types for clusters based on the provided metadata file.


s:about: |
  usage: sc_ctype_assign.R
        [-h] --query QUERY --celltypes CELLTYPES --source SOURCE --target
        TARGET [--diffgenes] [--diffpeaks] [--rnalogfc RNALOGFC]
        [--rnaminpct RNAMINPCT] [--rnaonlypos]
        [--rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--ataclogfc ATACLOGFC] [--atacminpct ATACMINPCT]
        [--atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--fragments FRAGMENTS] [--genes [GENES [GENES ...]]] [--pdf]
        [--verbose] [--h5seurat] [--h5ad] [--cbbuild] [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell Manual Cell Type Assignment

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression and/or chromatin
                          accessibility information stored in the RNA and ATAC
                          assays correspondingly. Additionally, 'rnaumap',
                          and/or 'atacumap', and/or 'wnnumap' dimensionality
                          reductions should be present.
    --celltypes CELLTYPES
                          Path to the TSV/CSV file for manual cell type
                          assignment for each of the clusters. First column -
                          'cluster', second column may have arbitrary name.
    --source SOURCE       Column from the metadata of the loaded Seurat object
                          to select clusters from.
    --target TARGET       Column from the metadata of the loaded Seurat object
                          to save manually assigned cell types. Should start
                          with 'custom_', otherwise, it won't be shown in UCSC
                          Cell Browser.
    --diffgenes           Identify differentially expressed genes (putative gene
                          markers) for assigned cell types. Ignored if loaded
                          Seurat object doesn't include genes expression
                          information stored in the RNA assay. Default: false
    --diffpeaks           Identify differentially accessible peaks for assigned
                          cell types. Ignored if loaded Seurat object doesn't
                          include chromatin accessibility information stored in
                          the ATAC assay. Default: false
    --rnalogfc RNALOGFC   For putative gene markers identification include only
                          those genes that on average have log fold change
                          difference in expression between every tested pair of
                          cell types not lower than this value. Ignored if '--
                          diffgenes' is not set or RNA assay is not present.
                          Default: 0.25
    --rnaminpct RNAMINPCT
                          For putative gene markers identification include only
                          those genes that are detected in not lower than this
                          fraction of cells in either of the two tested cell
                          types. Ignored if '--diffgenes' is not set or RNA
                          assay is not present. Default: 0.1
    --rnaonlypos          For putative gene markers identification return only
                          positive markers. Ignored if '--diffgenes' is not set
                          or RNA assay is not present. Default: false
    --rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for putative gene markers
                          identification. Ignored if '--diffgenes' is not set or
                          RNA assay is not present. Default: wilcox
    --ataclogfc ATACLOGFC
                          For differentially accessible peaks identification
                          include only those peaks that on average have log fold
                          change difference in the chromatin accessibility
                          between every tested pair of cell types not lower than
                          this value. Ignored if '--diffpeaks' is not set or
                          ATAC assay is not present. Default: 0.25
    --atacminpct ATACMINPCT
                          For differentially accessible peaks identification
                          include only those peaks that are detected in not
                          lower than this fraction of cells in either of the two
                          tested cell types. Ignored if '--diffpeaks' is not set
                          or ATAC assay is not present. Default: 0.05
    --atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differentially accessible
                          peaks identification. Ignored if '--diffpeaks' is not
                          set or ATAC assay is not present. Default: LR
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the loaded Seurat object. File should be saved
                          in TSV format with tbi-index file. Ignored if the
                          loaded Seurat object doesn't include ATAC assay.
    --genes [GENES [GENES ...]]
                          Genes of interest to build gene expression and/or Tn5
                          insertion frequency plots for the nearest peaks. To
                          build gene expression plots the loaded Seurat object
                          should include RNA assay. To build Tn5 insertion
                          frequency plots for the nearest peaks the loaded
                          Seurat object should include ATAC assay as well as the
                          --fragments file should be provided. Default: None
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