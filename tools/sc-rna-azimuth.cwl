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
      Path to the RDS file to load the query Seurat
      object from. This file should include genes
      expression information stored in the RNA assay
      and, optionally, chromatin accessibility
      information stored in the ATAC assay. The later
      is used only for plots.

  reference_data_rds:
    type: File
    inputBinding:
      prefix: "--reference"
    doc: |
      Path to the RDS file to load the reference Seurat
      object from. This file can be downloaded as
      ref.Rds from the
      https://azimuth.hubmapconsortium.org/references/

  reference_data_index:
    type: File
    inputBinding:
      prefix: "--annoyidx"
    doc: |
      Path to the annoy index file generated for the
      provided reference Seurat object. This file can
      be downloaded as idx.annoy from the
      https://azimuth.hubmapconsortium.org/references/

  reference_source_column:
    type: string
    inputBinding:
      prefix: "--source"
    doc: |
      Column from the metadata of the reference Seurat
      object to select the reference annotations.

  minimum_confidence_score:
    type: float?
    inputBinding:
      prefix: "--minconfscore"
    doc: |
      The minimum threshold for a prediction
      confidence score is calculated at the cell
      level. This metric ranges from 0 to 1 and
      reflects the confidence associated with each
      annotation. Only cells that meet both the
      minimum prediction confidence score and the
      minimum prediction mapping score thresholds
      will be included in the analysis.
      Default: 0.75

  minimum_mapping_score:
    type: float?
    inputBinding:
      prefix: "--minmapscore"
    doc: |
      The minimum threshold for a prediction
      mapping score is calculated at the cell.
      This metric ranges from 0 to 1 and reflects
      how well the unique structure of a cell's
      local neighborhood is preserved during
      reference mapping. Only cells that meet both
      the minimum prediction mapping score and the
      minimum prediction confidence score thresholds
      will be included in the analysis.
      Default: 0.75

  identify_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--diffgenes"
    doc: |
      Identify differentially expressed genes (putative
      gene markers) for the predicted cell types.
      Default: false

  identify_diff_peaks:
    type: boolean?
    inputBinding:
      prefix: "--diffpeaks"
    doc: |
      Identify differentially accessible peaks for the
      predicted cell types. Ignored if the query Seurat
      object doesn't include chromatin accessibility
      information stored in the ATAC assay. 
      Default: false

  rna_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--rnalogfc"
    doc: |
      For putative gene markers identification include only
      those genes that on average have a log fold change
      difference in the expression between every tested
      pair of the predicted cell types not lower than this
      value. Ignored if '--diffgenes' is not set.
      Default: 0.25

  rna_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--rnaminpct"
    doc: |
      For putative gene markers identification include only
      those genes that are detected in not lower than this
      fraction of cells in either of the two tested predicted
      cell types. Ignored if '--diffgenes' is not set.
      Default: 0.1

  only_positive_diff_genes:
    type: boolean?
    inputBinding:
      prefix: "--rnaonlypos"
    doc: |
      For putative gene markers identification return only
      upregulated markers. Ignored if '--diffgenes' is not
      set. Default: false

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
      Statistical test to use for putative gene markers
      identification. Ignored if '--diffgenes' is not set.
      Default: wilcox

  atac_minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--ataclogfc"
    doc: |
      For differentially accessible peaks identification include
      only those peaks that on average have a log fold change
      difference in the chromatin accessibility between every
      tested pair of the predicted cell types not lower than this
      value. Ignored if '--diffpeaks' is not set or if the query
      Seurat object doesn't include ATAC assay.
      Default: 0.25

  atac_minimum_pct:
    type: float?
    inputBinding:
      prefix: "--atacminpct"
    doc: |
      For differentially accessible peaks identification include
      only those peaks that are detected in not lower than this
      fraction of cells in either of the two tested predicted
      cell types. Ignored if '--diffpeaks' is not set or if the
      query Seurat object doesn't include ATAC assay.
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
      Statistical test to use for differentially accessible peaks
      identification. Ignored if '--diffpeaks' is not set or if
      the query Seurat object doesn't include ATAC assay.
      Default: LR

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment
      used in the query Seurat object. File should be saved
      in TSV format with a tbi-index file. Ignored if the
      query Seurat object doesn't include ATAC assay.

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to build gene expression and Tn5
      insertion frequency plots for the nearest peaks. To
      build Tn5 insertion frequency plots for the nearest
      peaks, the query Seurat object should include ATAC assay
      and the --fragments file should be provided.
      Default: None

  cvrg_upstream_bp:
    type: int?
    inputBinding:
      prefix: "--upstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene upstream. Ignored if --genes or --fragments
      parameters are not provided or when the query Seurat object
      doesn't include ATAC assay.
      Default: 2500

  cvrg_downstream_bp:
    type: int?
    inputBinding:
      prefix: "--downstream"
    doc: |
      Number of bases to extend the genome coverage region for
      a specific gene downstream. Ignored if --genes or --fragments
      parameters are not provided or when the query Seurat object
      doesn't include ATAC assay.
      Default: 2500

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
      Color theme for all generated plots.
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
      Save raw counts from the RNA and ATAC (if present)
      assays to h5ad file(s).
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
      Save Seurat data to SCope compatible loom file. Only not
      normalized raw counts from the RNA assay will be saved.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
    doc: |
      Export results to UCSC Cell Browser.
      Default: false

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
      Maximum memory in GB allowed to be shared between
      the workers when using multiple --cpus.
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

  ref_cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_ref_cell_cnts_gr_ctyp.png"
    doc: |
      Number of cells per cell type.
      All reference cells.
      PNG format.

  ref_umap_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_ref_umap_gr_ctyp.png"
    doc: |
      Reference UMAP colored by cell type.
      All reference cells.
      PNG format.

  cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_ref]_cell_cnts_gr_ctyp.png"
    doc: |
      Number of cells per cell type.
      All query cells.
      PNG format.

  umap_cnf_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_cnf.png"
    doc: |
      Projected UMAP colored by
      prediction confidence score.
      All query cells.
      PNG format.

  umap_map_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_map.png"
    doc: |
      Projected UMAP colored by
      prediction mapping score.
      All query cells.
      PNG format.

  qc_mtrcs_dnst_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_qc_mtrcs_dnst_gr_ctyp.png"
    doc: |
      Distribution of QC metrics per
      cell colored by cell type.
      All query cells.
      PNG format.

  gene_umi_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_umi_gr_cnf_spl_ctyp.png"
    doc: |
      Genes vs RNA reads per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.

  gene_umi_gr_map_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_gene_umi_gr_map_spl_ctyp.png"
    doc: |
      Genes vs RNA reads per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.

  umi_mito_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_umi_mito_gr_cnf_spl_ctyp.png"
    doc: |
      RNA reads vs mitochondrial % per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.

  umi_mito_gr_map_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_umi_mito_gr_map_spl_ctyp.png"
    doc: |
      RNA reads vs mitochondrial % per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.

  tss_frgm_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_tss_frgm_gr_cnf_spl_ctyp.png"
    doc: |
      TSS enrichment score vs ATAC fragments
      in peaks per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.

  tss_frgm_gr_map_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_tss_frgm_gr_map_spl_ctyp.png"
    doc: |
      TSS enrichment score vs ATAC fragments
      in peaks per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.

  rna_atac_cnts_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_rna_atac_cnts_gr_cnf_spl_ctyp.png"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      All query cells; split by cell type; colored
      by prediction confidence score.
      PNG format.

  rna_atac_cnts_gr_map_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_rna_atac_cnts_gr_map_spl_ctyp.png"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      All query cells; split by cell type; colored
      by prediction mapping score.
      PNG format.

  rnadbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_rnadbl_gr_ctyp.png"
    doc: |
      Percentage of RNA doublets
      per cell type.
      All query cells.
      PNG format.

  atacdbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_atacdbl_gr_ctyp.png"
    doc: |
      Percentage of ATAC doublets
      per cell type.
      All query cells.
      PNG format.

  vrlpdbl_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_vrlpdbl_gr_ctyp.png"
    doc: |
      Percentage of RNA and ATAC
      doublets per cell type.
      All query cells.
      PNG format.

  umap_gr_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_ref]_umap_gr_ctyp.png"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells.
      PNG format.

  umap_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_idnt.png"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_idnt.png"
    doc: |
      Composition plot colored by cell type.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.

  cmp_gr_idnt_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_idnt_spl_ctyp.png"
    doc: |
      Composition plot colored by dataset.
      Filtered query cells; split by cell type;
      downsampled to the smallest dataset.
      PNG format.

  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ph_spl_idnt.png"
    doc: |
      Projected UMAP colored by cell cycle phase.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_idnt.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.

  umap_gr_ctyp_spl_ph_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_ph.png"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by cell
      cycle phase; downsampled to the
      smallest dataset (if multiple datasets
      are analyzed jointly).
      PNG format.

  cmp_gr_ph_spl_ctyp_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_ctyp.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by cell type;
      downsampled to the smallest dataset (if
      multiple datasets are analyzed jointly).
      PNG format.

  umap_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ctyp_spl_cnd.png"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by grouping
      condition; first downsampled to the
      smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ctyp_spl_cnd.png"
    doc: |
      Composition plot colored by cell type.
      Filtered query cells; split by grouping
      condition; first downsampled to the
      smallest dataset, then downsampled to
      the smallest group.
      PNG format.

  cmp_gr_cnd_spl_ctyp_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_cnd_spl_ctyp.png"
    doc: |
      Composition plot colored by grouping condition.
      Filtered query cells; split by cell type;
      first downsampled to the smallest dataset,
      then downsampled to the smallest group.
      PNG format.

  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_gr_ph_spl_cnd.png"
    doc: |
      Projected UMAP colored by cell cycle phase.
      Filtered query cells; split by grouping
      condition; first downsampled to the smallest
      dataset, then downsampled to the smallest group.
      PNG format.

  cmp_gr_ph_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_cmp_gr_ph_spl_cnd.png"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by grouping condition;
      first downsampled to the smallest dataset, then
      downsampled to the smallest group.
      PNG format.

  xpr_avg_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_avg.png"
    doc: |
      Average gene expression.
      Filtered query cells.
      PNG format.

  xpr_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_dnst.png"
    doc: |
      Gene expression density.
      Filtered query cells.
      PNG format.

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_[!sgnl_]*.png"
    doc: |
      Projected UMAP colored by gene expression.
      Filtered query cells; all genes of interest.
      PNG format.

  xpr_per_cell_sgnl_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_sgnl_*.png"
    doc: |
      Projected UMAP colored by gene
      expression density.
      Filtered query cells; all genes
      of interest.
      PNG format.

  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.png"
    doc: |
      Gene expression heatmap from
      the filtered query cells.
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
      Filtered query cells;
      all genes of interest.
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
      Gene expression heatmap from
      the filtered query cells.
      Top gene markers.
      TSV format.

  gene_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_gene_markers.tsv"
    doc: |
      Gene markers from the filtered
      query cells.
      TSV format.

  peak_markers_tsv:
    type: File?
    outputBinding:
      glob: "*_peak_markers.tsv"
    doc: |
      Peak markers from the filtered
      query cells.
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
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_rna_azimuth.R"]:"/usr/local/bin/sc_rna_azimuth.R")


stdout: sc_rna_azimuth_stdout.log
stderr: sc_rna_azimuth_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell RNA-Seq Reference Mapping"
s:name: "Single-Cell RNA-Seq Reference Mapping"
s:alternateName: "Single-Cell RNA-Seq Reference Mapping"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-azimuth.cwl
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
  Single-Cell RNA-Seq Reference Mapping

  Predicts cell types on the cell level based on
  the reference annotation using Azimuth R package.
  Reference models can be downloaded from the
  https://azimuth.hubmapconsortium.org/


s:about: |
  usage: /usr/local/bin/sc_rna_azimuth.R [-h] --query QUERY --reference
                                        REFERENCE --annoyidx ANNOYIDX --source
                                        SOURCE [--minconfscore MINCONFSCORE]
                                        [--minmapscore MINMAPSCORE]
                                        [--diffgenes] [--diffpeaks]
                                        [--rnalogfc RNALOGFC]
                                        [--rnaminpct RNAMINPCT] [--rnaonlypos]
                                        [--rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                        [--ataclogfc ATACLOGFC]
                                        [--atacminpct ATACMINPCT]
                                        [--atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                                        [--fragments FRAGMENTS]
                                        [--genes [GENES [GENES ...]]]
                                        [--upstream UPSTREAM]
                                        [--downstream DOWNSTREAM] [--pdf]
                                        [--verbose] [--h5seurat] [--h5ad]
                                        [--loupe] [--cbbuild] [--scope]
                                        [--tmpdir TMPDIR] [--output OUTPUT]
                                        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                        [--cpus CPUS] [--memory MEMORY]
                                        [--seed SEED]

  Single-Cell RNA-Seq Reference Mapping

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load the query Seurat object
                          from. This file should include genes expression
                          information stored in the RNA assay and, optionally,
                          chromatin accessibility information stored in the ATAC
                          assay. The later is used only for plots.
    --reference REFERENCE
                          Path to the RDS file to load the reference Seurat
                          object from. This file can be downloaded as ref.Rds
                          from the
                          https://azimuth.hubmapconsortium.org/references/
    --annoyidx ANNOYIDX   Path to the annoy index file generated for the
                          provided reference Seurat object. This file can be
                          downloaded as idx.annoy from the
                          https://azimuth.hubmapconsortium.org/references/
    --source SOURCE       Column from the metadata of the reference Seurat
                          object to select the reference annotations.
    --minconfscore MINCONFSCORE
                          The minimum threshold for a prediction confidence
                          score is calculated at the cell level. This metric
                          ranges from 0 to 1 and reflects the confidence
                          associated with each annotation. Only cells that meet
                          both the minimum prediction confidence score and the
                          minimum prediction mapping score thresholds will be
                          included in the analysis. Default: 0.75
    --minmapscore MINMAPSCORE
                          The minimum threshold for a prediction mapping score
                          is calculated at the cell. This metric ranges from 0
                          to 1 and reflects how well the unique structure of a
                          cell’s local neighborhood is preserved during
                          reference mapping. Only cells that meet both the
                          minimum prediction mapping score and the minimum
                          prediction confidence score thresholds will be
                          included in the analysis. Default: 0.75
    --diffgenes           Identify differentially expressed genes (putative gene
                          markers) for the predicted cell types. Default: false
    --diffpeaks           Identify differentially accessible peaks for the
                          predicted cell types. Ignored if the query Seurat
                          object doesn't include chromatin accessibility
                          information stored in the ATAC assay. Default: false
    --rnalogfc RNALOGFC   For putative gene markers identification include only
                          those genes that on average have a log fold change
                          difference in the expression between every tested pair
                          of the predicted cell types not lower than this value.
                          Ignored if '--diffgenes is not set. Default: 0.25
    --rnaminpct RNAMINPCT
                          For putative gene markers identification include only
                          those genes that are detected in not lower than this
                          fraction of cells in either of the two tested
                          predicted cell types. Ignored if '--diffgenes' is not
                          set. Default: 0.1
    --rnaonlypos          For putative gene markers identification return only
                          upregulated markers. Ignored if '--diffgenes' is not
                          set. Default: false
    --rnatestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for putative gene markers
                          identification. Ignored if '--diffgenes' is not set.
                          Default: wilcox
    --ataclogfc ATACLOGFC
                          For differentially accessible peaks identification
                          include only those peaks that on average have a log
                          fold change difference in the chromatin accessibility
                          between every tested pair of the predicted cell types
                          not lower than this value. Ignored if '--diffpeaks is
                          not set or if the query Seurat object doesn't include
                          ATAC assay. Default: 0.25
    --atacminpct ATACMINPCT
                          For differentially accessible peaks identification
                          include only those peaks that are detected in not
                          lower than this fraction of cells in either of the two
                          tested predicted cell types. Ignored if '--diffpeaks'
                          is not set or if the query Seurat object doesn't
                          include ATAC assay. Default: 0.05
    --atactestuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differentially accessible
                          peaks identification. Ignored if '--diffpeaks' is not
                          set or if the query Seurat object doesn't include ATAC
                          assay. Default: LR
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the query Seurat object. File should be saved
                          in TSV format with tbi-index file. Ignored if the
                          query Seurat object doesn't include ATAC assay.
    --genes [GENES [GENES ...]]
                          Genes of interest to build gene expression and Tn5
                          insertion frequency plots for the nearest peaks. To
                          build Tn5 insertion frequency plots for the nearest
                          peaks the query Seurat object should include ATAC
                          assay as well as the --fragments file should be
                          provided. Default: None
    --upstream UPSTREAM   Number of bases to extend the genome coverage region
                          for a specific gene upstream. Ignored if --genes or
                          --fragments parameters are not provided or when the
                          query Seurat object doesn't include ATAC assay.
                          Default: 2500
    --downstream DOWNSTREAM
                          Number of bases to extend the genome coverage region
                          for a specific gene downstream. Ignored if --genes or
                          --fragments parameters are not provided or when the
                          query Seurat object doesn't include ATAC assay.
                          Default: 2500
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA and ATAC (if present)
                          assays to h5ad file(s). Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --scope               Save Seurat data to SCope compatible loom file. Only
                          not normalized raw counts from the RNA assay will be
                          saved. Default: false
    --tmpdir TMPDIR       Directory to keep temporary files. Default: either
                          /tmp or defined by the environment variables TMPDIR,
                          TMP, TEMP.
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --seed SEED           Seed number for random values. Default: 42