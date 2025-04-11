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
      and ATAC assays correspondingly.

  cell_type_data:
    type: File
    inputBinding:
      prefix: "--celltypes"
    doc: |
      Path to the TSV/CSV file for manual cell type assignment for each of the clusters.
      First column - 'cluster', second column may have arbitrary name.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata
      by the selected barcodes. First column should be named as 'barcode'.
      Other columns will be added to the Seurat object metadata ovewriting
      the existing ones if those are present.
      Default: no extra metadata is added

  query_source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column from the metadata of the loaded Seurat object to
      select clusters from. May be one of the columns added
      with the --barcodes parameter.

  query_target_column:
    type: string
    inputBinding:
      prefix: "--target"
    doc: |
      Column from the metadata of the loaded Seurat object to save manually
      assigned cell types. Should start with 'custom_', otherwise, it won't
      be shown in UCSC Cell Browser.

  reduction:
    type: string?
    inputBinding:
      prefix: "--reduction"
    doc: |
      Dimensionality reduction to be used in the generated plots. If not
      provided it will be automatically defined on the basis of the --source
      parameter as follows: rna_res.* - rnaumap, atac_res.* - atacumap,
      wsnn_res.* - wnnumap.
      Default: defined automatically

  query_splitby_column:
    type: string?
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to additionally split
      every cluster selected with --source into smaller groups.
      May be one of the columns added with the --barcodes parameter.
      Default: do not split

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

  genesets_data:
    type: File?
    inputBinding:
      prefix: "--genesets"
    doc: |
      Path to the GMT file for calculating average expression levels
      (module scores) per gene set. This file can be downloaded from
      the Molecular Signatures Database (MSigDB) following the link
      https://www.gsea-msigdb.org/gsea/msigdb. To calculate module
      scores the loaded Seurat object should include RNA assay.
      Default: do not calculate gene set expression scores.

  cvrg_upstream_bp:
    type: int?
    inputBinding:
      prefix: "--upstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene upstream. Ignored if --genes or --fragments
      parameters are not provided. Default: 2500

  cvrg_downstream_bp:
    type: int?
    inputBinding:
      prefix: "--downstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene downstream. Ignored if --genes or --fragments
      parameters are not provided. Default: 2500

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
      Save raw counts from the RNA and/or ATAC assay(s) to h5ad file(s).
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

  export_azimuth_ref:
    type: boolean?
    inputBinding:
      prefix: "--azimuth"
    doc: |
      Save Seurat object with the assigned cell
      types as model for the reference mapping
      in Azimuth. Both RDS and annoy index files
      will be created.
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

  cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cell_cnts_gr_ctyp.png"
    doc: |
      Number of cells per cell type.
      All cells.
      PNG format.

  gene_umi_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_umi_spl_ctyp.png"
    doc: |
      Genes vs RNA reads per cell.
      Split by cell type; all cells.
      PNG format.

  umi_mito_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_umi_mito_spl_ctyp.png"
    doc: |
      RNA reads vs mitochondrial % per cell.
      Split by cell type; all cells.
      PNG format.

  rnadbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_rnadbl_gr_ctyp.png"
    doc: |
      Percentage of RNA doublets per cell type.
      All cells.
      PNG format.

  tss_frgm_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_tss_frgm_spl_ctyp.png"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Split by cell type; all cells.
      PNG format.

  atacdbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_atacdbl_gr_ctyp.png"
    doc: |
      Percentage of ATAC doublets per cell type.
      All cells.
      PNG format.

  rna_atac_cnts_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_rna_atac_cnts_spl_ctyp.png"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      Split by cell type; all cells.
      PNG format.

  vrlpdbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlpdbl_gr_ctyp.png"
    doc: |
      Percentage of RNA and ATAC doublets
      per cell type.
      All cells.
      PNG format.

  qc_mtrcs_dnst_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_qc_mtrcs_dnst_gr_ctyp.png"
    doc: |
      Distribution of QC metrics per cell
      colored by cell type.
      All cells.
      PNG format.

  umap_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp.png"
    doc: |
      UMAP colored by cell type.
      All cells.
      PNG format.

  umap_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_idnt.png"
    doc: |
      UMAP colored by cell type.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_idnt.png"
    doc: |
      Composition plot colored by cell type.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.

  cmp_gr_idnt_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_ctyp.png"
    doc: |
      Composition plot colored by dataset.
      Split by cell type; downsampled to
      the smallest dataset.
      PNG format.

  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ph_spl_idnt.png"
    doc: |
      UMAP colored by cell cycle phase.
      Split by dataset; downsampled to the
      smallest dataset.
      PNG format.

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by dataset; downsampled to the smallest
      dataset.
      PNG format.

  umap_gr_ctyp_spl_ph_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_ph.png"
    doc: |
      UMAP colored by cell type.
      Split by cell cycle phase; downsampled
      to the smallest dataset (if multiple
      datasets are analyzed jointly).
      PNG format.

  cmp_gr_ph_spl_ctyp_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_ctyp.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by cell type; downsampled to the
      smallest dataset (if multiple datasets are
      analyzed jointly).
      PNG format.

  umap_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_cnd.png"
    doc: |
      UMAP colored by cell type.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_cnd.png"
    doc: |
      Composition plot colored by cell type.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  cmp_gr_cnd_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_ctyp.png"
    doc: |
      Composition plot colored by grouping condition.
      Split by cell type; first downsampled to the
      smallest dataset, then downsampled to the
      smallest group.
      PNG format.

  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ph_spl_cnd.png"
    doc: |
      UMAP colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  cmp_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_cnd.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Split by grouping condition; first downsampled
      to the smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  gse_per_cell_plot_png:
    type: File?
    outputBinding:
      glob: "*_gse_per_cell.png"
    doc: |
      UMAP colored by gene set expression score.
      PNG format.

  gse_avg_plot_png:
    type: File?
    outputBinding:
      glob: "*_gse_avg.png"
    doc: |
      Average gene set expression score.
      PNG format.

  gse_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_gse_dnst.png"
    doc: |
      Gene set expression score density.
      PNG format.

  xpr_avg_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_avg.png"
    doc: |
      Average gene expression.
      PNG format.

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.png"
    doc: |
      UMAP colored by gene expression.
      All genes of interest.
      PNG format.

  xpr_per_cell_sgnl_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.png"
    doc: |
      UMAP colored by gene expression density.
      All genes of interest.
      PNG format.

  xpr_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_dnst.png"
    doc: |
      Gene expression density.
      PNG format.

  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.png"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      PNG format.

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_cvrg_*.png"
    doc: |
      ATAC fragment coverage.
      All genes of interest.
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

  xpr_htmp_tsv:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.tsv"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      TSV format.

  xpr_htmp_gct:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.gct"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      GCT format.

  xpr_htmp_html:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.html"
    doc: |
      Gene expression heatmap.
      Top gene markers.
      HTML format.

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Gene markers.
      TSV format.

  peak_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_peak_markers.tsv"
    doc: |
      Peak markers.
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
      glob: "*[!_ref]_data.rds"
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
      Loupe format

  seurat_data_scope:
    type: File?
    outputBinding:
      glob: "*_data.loom"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  reference_data_rds:
    type: File?
    outputBinding:
      glob: "*_ref_data.rds"
    doc: |
      Seurat object with assigned cell
      types formatted as an Azimuth
      reference model.
      RDS format.

  reference_data_index:
    type: File?
    outputBinding:
      glob: "*_ref_data.annoy"
    doc: |
      Annoy index generated for the
      Azimuth reference model.
      Annoy format.

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

  human_log:
    type: File?
    outputBinding:
      glob: "error_report.txt"
    doc: |
      Human readable error log.
      TXT format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["Rscript"]
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_ctype_assign.R"]:"/usr/local/bin/sc_ctype_assign.R")


stdout: sc_ctype_assign_stdout.log
stderr: sc_ctype_assign_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell Manual Cell Type Assignment"
s:name: "Single-Cell Manual Cell Type Assignment"
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
  Single-Cell Manual Cell Type Assignment

  Assigns cell types for clusters based on
  the provided metadata file.


s:about: |
  usage: /usr/local/bin/sc_ctype_assign.R [-h] --query QUERY --celltypes
                                          CELLTYPES [--barcodes BARCODES]
                                          --source SOURCE --target TARGET
                                          [--splitby SPLITBY]
                                          [--reduction REDUCTION] [--diffgenes]
                                          [--diffpeaks] [--rnalogfc RNALOGFC]
                                          [--rnaminpct RNAMINPCT] [--rnaonlypos]
                                          [--rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                          [--ataclogfc ATACLOGFC]
                                          [--atacminpct ATACMINPCT]
                                          [--atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                          [--fragments FRAGMENTS]
                                          [--genes [GENES [GENES ...]]]
                                          [--genesets GENESETS]
                                          [--upstream UPSTREAM]
                                          [--downstream DOWNSTREAM] [--pdf]
                                          [--verbose] [--h5seurat] [--h5ad]
                                          [--loupe] [--azimuth] [--cbbuild]
                                          [--scope] [--output OUTPUT]
                                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                          [--cpus CPUS] [--memory MEMORY]
                                          [--seed SEED]

  Single-Cell Manual Cell Type Assignment

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression and/or chromatin
                          accessibility information stored in the RNA and ATAC
                          assays correspondingly.
    --celltypes CELLTYPES
                          Path to the TSV/CSV file for manual cell type
                          assignment for each of the clusters. First column -
                          'cluster', second column may have arbitrary name.
    --barcodes BARCODES   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata by the selected barcodes. First column
                          should be named as 'barcode'. Other columns will be
                          added to the Seurat object metadata ovewriting the
                          existing ones if those are present. Default: no extra
                          metadata is added
    --source SOURCE       Column from the metadata of the loaded Seurat object
                          to select clusters from. May be one of the columns
                          added with the --barcodes parameter.
    --target TARGET       Column from the metadata of the loaded Seurat object
                          to save manually assigned cell types. Should start
                          with 'custom_', otherwise, it won't be shown in UCSC
                          Cell Browser.
    --splitby SPLITBY     Column from the Seurat object metadata to additionally
                          split every cluster selected with --source into
                          smaller groups. May be one of the columns added with
                          the --barcodes parameter. Default: do not split
    --reduction REDUCTION
                          Dimensionality reduction to be used in the generated
                          plots. If not provided it will be automatically
                          defined on the basis of the --source parameter as
                          follows: rna_res.* - rnaumap, atac_res.* - atacumap,
                          wsnn_res.* - wnnumap. Default: defined automatically
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
    --genesets GENESETS   Path to the GMT file for calculating average
                          expression levels (module scores) per gene set. This
                          file can be downloaded from the Molecular Signatures
                          Database (MSigDB) following the link https://www.gsea-
                          msigdb.org/gsea/msigdb. To calculate module scores the
                          loaded Seurat object should include RNA assay.
                          Default: do not calculate gene set expression scores.
    --upstream UPSTREAM   Number of bases to extend the genome coverage region
                          for a specific gene upstream. Ignored if --genes or
                          --fragments parameters are not provided. Default: 2500
    --downstream DOWNSTREAM
                          Number of bases to extend the genome coverage region
                          for a specific gene downstream. Ignored if --genes or
                          --fragments parameters are not provided. Default: 2500
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA and/or ATAC assay(s) to
                          h5ad file(s). Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --azimuth             Save Seurat object with the assigned cell types as
                          model for the reference mapping in Azimuth. Both RDS
                          and annoy index files will be created. Default: false
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