cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: dummy_metadata.csv
    entry: |
      library_id
      Experiment
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.14


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate
      experiment in MEX format. The rows consist of all the genes and peaks concatenated
      together and the columns are restricted to those barcodes that are identified as cells.

  aggregation_metadata:
    type: File?
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger ARC Aggregate outputs, the aggr.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment observed in the experiment in TSV
      format. Tbi-index file is required.

  annotation_gtf_file:
    type: File
    inputBinding:
      prefix: "--annotations"
    doc: |
      Path to the genome annotation file in GTF format.

  grouping_data:
    type: File?
    inputBinding:
      prefix: "--grouping"
    doc: |
      Path to the TSV/CSV file to define datasets grouping.
      First column - 'library_id' with the values and order
      that correspond to the 'library_id' column from the '
      --identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  blacklist_regions_file:
    type: File?
    inputBinding:
      prefix: "--blacklist"
    doc: |
      Path to the optional BED file with the genomic blacklist regions.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

  rna_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--rnamincells"
    doc: |
      Include only genes detected in at least this many cells.
      Default: 5 (applied to all datasets)

  minimum_genes:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--mingenes"
    doc: |
      Include cells where at least this many genes are detected. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file.
      Default: 250 (applied to all datasets)

  maximum_genes:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--maxgenes"
    doc: |
      Include cells with the number of genes not bigger than this value. If multiple
      values provided, each of them will be applied to the correspondent dataset from
      the '--mex' input based on the '--identity' file.
      Default: 5000 (applied to all datasets)

  rna_minimum_umi:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--rnaminumi"
    doc: |
      Include cells where at least this many UMI (RNA transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 500 (applied to all datasets)

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^mt-|^MT-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5 (applied to all datasets)

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      as log10(genes)/log10(UMI) for RNA assay. If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)

  atac_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--atacmincells"
    doc: |
      Include only peaks detected in at least this many cells.
      Default: 5 (applied to all datasets)

  atac_minimum_umi:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--atacminumi"
    doc: |
      Include cells where at least this many UMI (ATAC transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 1000 (applied to all datasets)

  maximum_nucl_signal:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--maxnuclsignal"
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments. If multiple values provided, each of
      them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 4 (applied to all datasets)

  minimum_tss_enrich:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--mintssenrich"
    doc: |
      Include cells with the TSS enrichment score not lower than this value.
      Score is calculated based on the ratio of fragments centered at the TSS
      to fragments in TSS-flanking regions. If multiple values provided, each
      of them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 2 (applied to all datasets)

  minimum_frip:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minfrip"
    doc: |
      Include cells with the FRiP not lower than this value. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file. FRiP is calculated for fragments.
      Default: 0.15 (applied to all datasets)

  maximum_blacklist_fraction:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--maxblacklist"
    doc: |
      Include cells with the fraction of fragments in
      genomic blacklist regions not bigger than this value.
      If multiple values provided, each of them will be
      applied to the correspondent dataset from the '--mex'
      input based on the '--identity' file.
      Default: 0.05 (applied to all datasets)

  call_by:
    type: string?
    inputBinding:
      prefix: "--callby"
    doc: |
      Replace Cell Ranger ARC peaks with MACS2 peaks called
      for cells grouped by the column from the optionally
      provided --barcodes file. If --barcodes file was not
      provided MACS2 peaks can be still called per dataset
      by setting --callby to new.ident. Peaks are called
      only after applying all RNA related thresholds,
      maximum nucleosome signal, and minimum TSS enrichment
      scores filters.
      Default: do not call peaks

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

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.png"
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PNG format

  raw_1_2_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.pdf"
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PDF format

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.png"
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PNG format

  raw_2_3_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.pdf"
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PDF format

  raw_cells_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cells_count.png"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_cells_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cells_count.pdf"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format

  raw_rna_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rna_umi_dnst.png"
    doc: |
      UMI per cell density for RNA assay (not filtered).
      PNG format

  raw_rna_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_rna_umi_dnst.pdf"
    doc: |
      UMI per cell density for RNA assay (not filtered).
      PDF format

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Genes per cell density (not filtered).
      PNG format

  raw_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.pdf"
    doc: |
      Genes per cell density (not filtered).
      PDF format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.png"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (not filtered).
      PNG format

  raw_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.pdf"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (not filtered).
      PDF format

  raw_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.png"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.pdf"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PDF format

  raw_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.png"
    doc: |
      Novelty score per cell density for RNA assay (not filtered).
      PNG format

  raw_nvlt_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.pdf"
    doc: |
      Novelty score per cell density for RNA assay (not filtered).
      PDF format

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.png"
    doc: |
      UMI per cell density for ATAC assay (not filtered).
      PNG format

  raw_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst.pdf"
    doc: |
      UMI per cell density for ATAC assay (not filtered).
      PDF format

  raw_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.png"
    doc: |
      Peaks per cell density (not filtered).
      PNG format

  raw_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.pdf"
    doc: |
      Peaks per cell density (not filtered).
      PDF format

  raw_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst.png"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered).
      PNG format

  raw_blck_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst.pdf"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered).
      PDF format

  raw_rna_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rna_atac_umi_corr.png"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (not filtered).
      PNG format

  raw_rna_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_rna_atac_umi_corr.pdf"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (not filtered).
      PDF format

  raw_tss_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_tss_atac_umi_corr.png"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (not filtered).
      PNG format

  raw_tss_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_tss_atac_umi_corr.pdf"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (not filtered).
      PDF format

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.png"
    doc: |
      QC metrics per cell density (not filtered).
      PNG format

  raw_qc_mtrcs_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.pdf"
    doc: |
      QC metrics per cell density (not filtered).
      PDF format

  raw_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_tss_nrch.png"
    doc: |
      TSS enrichment score (not filtered).
      PNG format

  raw_tss_nrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_tss_nrch.pdf"
    doc: |
      TSS enrichment score (not filtered).
      PDF format

  raw_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*_raw_frgm_hist.png"
    doc: |
      Fragments length histogram (not filtered).
      PNG format

  raw_frgm_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_frgm_hist.pdf"
    doc: |
      Fragments length histogram (not filtered).
      PDF format

  raw_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rna_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (not filtered).
      PNG format

  raw_rna_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_rna_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (not filtered).
      PDF format

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PNG format

  raw_gene_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PDF format

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PDF format

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (not filtered).
      PNG format

  raw_nvlt_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (not filtered).
      PDF format

  raw_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (not filtered).
      PNG format

  raw_atac_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_atac_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (not filtered).
      PDF format

  raw_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition peaks per cell density (not filtered).
      PNG format

  raw_peak_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition peaks per cell density (not filtered).
      PDF format

  raw_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (not filtered).
      PNG format

  raw_blck_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (not filtered).
      PDF format

  mid_fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      PC1 and PC2 from the QC metrics PCA (intermediate filtered).
      PNG format

  mid_fltr_1_2_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_1_2_qc_mtrcs_pca.pdf"
    doc: |
      PC1 and PC2 from the QC metrics PCA (intermediate filtered).
      PDF format

  mid_fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      PC2 and PC3 from the QC metrics PCA (intermediate filtered).
      PNG format

  mid_fltr_2_3_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_2_3_qc_mtrcs_pca.pdf"
    doc: |
      PC2 and PC3 from the QC metrics PCA (intermediate filtered).
      PDF format

  mid_fltr_cells_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_cells_count.png"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PNG format

  mid_fltr_cells_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_cells_count.pdf"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PDF format

  mid_fltr_rna_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_umi_dnst.png"
    doc: |
      UMI per cell density for RNA assay (intermediate filtered).
      PNG format

  mid_fltr_rna_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_umi_dnst.pdf"
    doc: |
      UMI per cell density for RNA assay (intermediate filtered).
      PDF format

  mid_fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst.png"
    doc: |
      Genes per cell density (intermediate filtered).
      PNG format

  mid_fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst.pdf"
    doc: |
      Genes per cell density (intermediate filtered).
      PDF format

  mid_fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_umi_corr.png"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (intermediate filtered).
      PNG format

  mid_fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (intermediate filtered).
      PDF format

  mid_fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst.png"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (intermediate filtered).
      PNG format

  mid_fltr_mito_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst.pdf"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (intermediate filtered).
      PDF format

  mid_fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst.png"
    doc: |
      Novelty score per cell density for RNA assay (intermediate filtered).
      PNG format

  mid_fltr_nvlt_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst.pdf"
    doc: |
      Novelty score per cell density for RNA assay (intermediate filtered).
      PDF format

  mid_fltr_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst.png"
    doc: |
      UMI per cell density for ATAC assay (intermediate filtered).
      PNG format

  mid_fltr_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst.pdf"
    doc: |
      UMI per cell density for ATAC assay (intermediate filtered).
      PDF format

  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst.png"
    doc: |
      Peaks per cell density (intermediate filtered).
      PNG format

  mid_fltr_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst.pdf"
    doc: |
      Peaks per cell density (intermediate filtered).
      PDF format

  mid_fltr_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst.png"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (intermediate filtered).
      PNG format

  mid_fltr_blck_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst.pdf"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (intermediate filtered).
      PDF format

  mid_fltr_rna_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_atac_umi_corr.png"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (intermediate filtered).
      PNG format

  mid_fltr_rna_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_atac_umi_corr.pdf"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (intermediate filtered).
      PDF format

  mid_fltr_tss_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_atac_umi_corr.png"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (intermediate filtered).
      PNG format

  mid_fltr_tss_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_atac_umi_corr.pdf"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (intermediate filtered).
      PDF format

  mid_fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_qc_mtrcs_dnst.png"
    doc: |
      QC metrics per cell density (intermediate filtered).
      PNG format

  mid_fltr_qc_mtrcs_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_qc_mtrcs_dnst.pdf"
    doc: |
      QC metrics per cell density (intermediate filtered).
      PDF format

  mid_fltr_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_nrch.png"
    doc: |
      TSS enrichment score (intermediate filtered).
      PNG format

  mid_fltr_tss_nrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_nrch.pdf"
    doc: |
      TSS enrichment score (intermediate filtered).
      PDF format

  mid_fltr_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frgm_hist.png"
    doc: |
      Fragments length histogram (intermediate filtered).
      PNG format

  mid_fltr_frgm_hist_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frgm_hist.pdf"
    doc: |
      Fragments length histogram (intermediate filtered).
      PDF format

  mid_fltr_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (intermediate filtered).
      PNG format

  mid_fltr_rna_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (intermediate filtered).
      PDF format

  mid_fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition genes per cell density (intermediate filtered).
      PNG format

  mid_fltr_gene_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition genes per cell density (intermediate filtered).
      PDF format

  mid_fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (intermediate filtered).
      PNG format

  mid_fltr_mito_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (intermediate filtered).
      PDF format

  mid_fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (intermediate filtered).
      PNG format

  mid_fltr_nvlt_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (intermediate filtered).
      PDF format

  mid_fltr_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (intermediate filtered).
      PNG format

  mid_fltr_atac_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atac_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (intermediate filtered).
      PDF format

  mid_fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition peaks per cell density (intermediate filtered).
      PNG format

  mid_fltr_peak_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition peaks per cell density (intermediate filtered).
      PDF format

  mid_fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (intermediate filtered).
      PNG format

  mid_fltr_blck_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (intermediate filtered).
      PDF format

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PNG format

  fltr_1_2_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_1_2_qc_mtrcs_pca.pdf"
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PDF format

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PNG format

  fltr_2_3_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_2_3_qc_mtrcs_pca.pdf"
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PDF format

  fltr_cells_count_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_cells_count.png"
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_cells_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_cells_count.pdf"
    doc: |
      Number of cells per dataset (filtered).
      PDF format

  fltr_rna_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_umi_dnst.png"
    doc: |
      UMI per cell density for RNA assay (filtered).
      PNG format

  fltr_rna_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_umi_dnst.pdf"
    doc: |
      UMI per cell density for RNA assay (filtered).
      PDF format

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst.png"
    doc: |
      Genes per cell density (filtered).
      PNG format

  fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst.pdf"
    doc: |
      Genes per cell density (filtered).
      PDF format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_umi_corr.png"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (filtered).
      PNG format

  fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (filtered).
      PDF format

  fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst.png"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst.pdf"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PDF format

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst.png"
    doc: |
      Novelty score per cell density for RNA assay (filtered).
      PNG format

  fltr_nvlt_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst.pdf"
    doc: |
      Novelty score per cell density for RNA assay (filtered).
      PDF format

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst.png"
    doc: |
      UMI per cell density for ATAC assay (filtered).
      PNG format

  fltr_atac_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst.pdf"
    doc: |
      UMI per cell density for ATAC assay (filtered).
      PDF format

  fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst.png"
    doc: |
      Peaks per cell density (filtered).
      PNG format

  fltr_peak_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst.pdf"
    doc: |
      Peaks per cell density (filtered).
      PDF format

  fltr_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst.png"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (filtered).
      PNG format

  fltr_blck_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst.pdf"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (filtered).
      PDF format

  fltr_rna_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_atac_umi_corr.png"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (filtered).
      PNG format

  fltr_rna_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_atac_umi_corr.pdf"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (filtered).
      PDF format

  fltr_tss_atac_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_atac_umi_corr.png"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (filtered).
      PNG format

  fltr_tss_atac_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_atac_umi_corr.pdf"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (filtered).
      PDF format

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_qc_mtrcs_dnst.png"
    doc: |
      QC metrics per cell density (filtered).
      PNG format

  fltr_qc_mtrcs_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_qc_mtrcs_dnst.pdf"
    doc: |
      QC metrics per cell density (filtered).
      PDF format

  fltr_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_nrch.png"
    doc: |
      TSS enrichment score (filtered).
      PNG format

  fltr_tss_nrch_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_nrch.pdf"
    doc: |
      TSS enrichment score (filtered).
      PDF format

  fltr_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frgm_hist.png"
    doc: |
      Fragments length histogram (filtered).
      PNG format

  fltr_frgm_hist_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frgm_hist.pdf"
    doc: |
      Fragments length histogram (filtered).
      PDF format

  fltr_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (filtered).
      PNG format

  fltr_rna_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (filtered).
      PDF format

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PNG format

  fltr_gene_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PDF format

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PDF format

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (filtered).
      PNG format

  fltr_nvlt_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (filtered).
      PDF format

  fltr_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (filtered).
      PNG format

  fltr_atac_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atac_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (filtered).
      PDF format

  fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition peaks per cell density (filtered).
      PNG format

  fltr_peak_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition peaks per cell density (filtered).
      PDF format

  fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (filtered).
      PNG format

  fltr_blck_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (filtered).
      PDF format

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
      Filtered Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Filtered Seurat data in h5seurat format

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


baseCommand: ["sc_multiome_filter.R"]
arguments:
- valueFrom: |
    ${
      if (inputs.aggregation_metadata) {
        return inputs.aggregation_metadata;
      } else {
        return runtime.outdir + "/dummy_metadata.csv"
      }
    }
  prefix: "--identity"


stdout: sc_multiome_filter_stdout.log
stderr: sc_multiome_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Multiome ATAC and RNA-Seq Filtering Analysis"
s:name: "Single-cell Multiome ATAC and RNA-Seq Filtering Analysis"
s:alternateName: "Filters single-cell multiome ATAC and RNA-Seq datasets based on the common QC metrics"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-multiome-filter.cwl
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
  Single-cell Multiome ATAC and RNA-Seq Filtering Analysis

  Filters single-cell multiome ATAC and RNA-Seq datasets based on the common QC metrics.


s:about: |
  usage: sc_multiome_filter.R
        [-h] --mex MEX --identity IDENTITY --fragments FRAGMENTS --annotations
        ANNOTATIONS [--grouping GROUPING] [--blacklist BLACKLIST]
        [--barcodes BARCODES] [--rnamincells RNAMINCELLS]
        [--mingenes [MINGENES ...]] [--maxgenes [MAXGENES ...]]
        [--rnaminumi [RNAMINUMI ...]] [--mitopattern MITOPATTERN]
        [--maxmt MAXMT] [--minnovelty [MINNOVELTY ...]]
        [--atacmincells ATACMINCELLS] [--atacminumi [ATACMINUMI ...]]
        [--maxnuclsignal [MAXNUCLSIGNAL ...]]
        [--mintssenrich [MINTSSENRICH ...]] [--minfrip [MINFRIP ...]]
        [--maxblacklist [MAXBLACKLIST ...]] [--callby CALLBY] [--pdf]
        [--verbose] [--h5seurat] [--h5ad] [--cbbuild] [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell Multiome ATAC and RNA-Seq Filtering Analysis

  options:
    -h, --help            show this help message and exit
    --mex MEX             Path to the folder with feature-barcode matrix from
                          Cell Ranger ARC Count/Aggregate experiment in MEX
                          format. The rows consist of all the genes and peaks
                          concatenated together and the columns are restricted
                          to those barcodes that are identified as cells.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities. If '--mex' points to the Cell Ranger ARC
                          Aggregate outputs, the aggr.csv file can be used. If
                          Cell Ranger ARC Count outputs have been used in the '
                          --mex' input, the file should include at least one
                          column - 'library_id' and one row with the alias for
                          Cell Ranger ARC Count experiment.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          observed in the experiment in TSV format. Tbi-index
                          file is required.
    --annotations ANNOTATIONS
                          Path to the genome annotation file in GTF format
    --grouping GROUPING   Path to the TSV/CSV file to define datasets grouping.
                          First column - 'library_id' with the values and order
                          that correspond to the 'library_id' column from the '
                          --identity' file, second column 'condition'. Default:
                          each dataset is assigned to its own group.
    --blacklist BLACKLIST
                          Path to the optional BED file with the genomic
                          blacklist regions.
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --rnamincells RNAMINCELLS
                          Include only genes detected in at least this many
                          cells. Default: 5 (applied to all datasets)
    --mingenes [MINGENES ...]
                          Include cells where at least this many genes are
                          detected. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES ...]
                          Include cells with the number of genes not bigger than
                          this value. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Default:
                          5000 (applied to all datasets)
    --rnaminumi [RNAMINUMI ...]
                          Include cells where at least this many UMI (RNA
                          transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the '--mex' input based on
                          the '--identity' file. Default: 500 (applied to all
                          datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial genes.
                          Default: '^mt-|^MT-'
    --maxmt MAXMT         Include cells with the percentage of transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5 (applied to all datasets)
    --minnovelty [MINNOVELTY ...]
                          Include cells with the novelty score not lower than
                          this value, calculated for as log10(genes)/log10(UMI)
                          for RNA assay. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the '--mex' input based on the '--identity' file.
                          Default: 0.8 (applied to all datasets)
    --atacmincells ATACMINCELLS
                          Include only peaks detected in at least this many
                          cells. Default: 5 (applied to all datasets)
    --atacminumi [ATACMINUMI ...]
                          Include cells where at least this many UMI (ATAC
                          transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the '--mex' input based on
                          the '--identity' file. Default: 1000 (applied to all
                          datasets)
    --maxnuclsignal [MAXNUCLSIGNAL ...]
                          Include cells with the nucleosome signal not bigger
                          than this value. Nucleosome signal quantifies the
                          approximate ratio of mononucleosomal to nucleosome-
                          free fragments. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the '--mex' input based on the '--identity' file.
                          Default: 4 (applied to all datasets)
    --mintssenrich [MINTSSENRICH ...]
                          Include cells with the TSS enrichment score not lower
                          than this value. Score is calculated based on the
                          ratio of fragments centered at the TSS to fragments in
                          TSS-flanking regions. If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the '--mex' input based on the '--
                          identity' file. Default: 2 (applied to all datasets)
    --minfrip [MINFRIP ...]
                          Include cells with the FRiP not lower than this value.
                          If multiple values provided, each of them will be
                          applied to the correspondent dataset from the '--mex'
                          input based on the '--identity' file. FRiP is
                          calculated for fragments. Default: 0.15 (applied to
                          all datasets)
    --maxblacklist [MAXBLACKLIST ...]
                          Include cells with the fraction of fragments in
                          genomic blacklist regions not bigger than this value.
                          If multiple values provided, each of them will be
                          applied to the correspondent dataset from the '--mex'
                          input based on the '--identity' file. Default: 0.05
                          (applied to all datasets)
    --callby CALLBY       Replace Cell Ranger ARC peaks with MACS2 peaks called
                          for cells grouped by the column from the optionally
                          provided --barcodes file. If --barcodes file was not
                          provided MACS2 peaks can be still called per dataset
                          by setting --callby to new.ident. Peaks are called
                          only after applying all RNA related thresholds,
                          maximum nucleosome signal, and minimum TSS enrichment
                          scores filters. Default: do not call peaks
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
                          workers when using multiple '--cpus'. Default: 32