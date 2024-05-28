cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.38


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
    inputBinding:
      prefix: "--identity"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities, if '--mex' points to
      the Cell Ranger ARC Aggregate outputs. The aggr.csv file can be used.

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

  chrom_length_file:
    type: File
    inputBinding:
      prefix: "--seqinfo"
    doc: |
      Path to the headerless chromosome length file in TSV format

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
    type:
    - "null"
    - File
    - type: enum
      symbols:
      - "hg19"
      - "hg38"
      - "mm10"
    inputBinding:
      prefix: "--blacklist"
      valueFrom: |
        ${
          if (self.class && self.class == "File"){
            return self;
          } else if (self == "hg19") {
            return "/opt/sc_tools/hg19-blacklist.v2.bed";
          } else if (self == "hg38") {
            return "/opt/sc_tools/hg38-blacklist.v2.bed";
          } else if (self == "mm10") {
            return "/opt/sc_tools/mm10-blacklist.v2.bed";
          } else {
            return null;
          }
        }
    doc: |
      Path to the optional BED file with the genomic blacklist regions.
      If a string value provided, it should be one of the hg19, hg38,
      or mm10 as we replace it with the file location from docker image

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
      '--mex' input based on the '--identity' file. Any 0 will be replaced with the
      auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.
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
      the '--mex' input based on the '--identity' file. Any 0 will be replaced with the
      auto-estimated threshold (median + 5 * MAD) calculated per dataset.
      Default: 5000 (applied to all datasets)

  minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--minumis"
    doc: |
      Include cells where at least this many RNA reads are detected.
      If multiple values provided, each of them will be applied to the
      correspondent dataset from the '--mex' input based on the '--identity' file.
      Any 0 will be replaced with the auto-estimated threshold (median - 2.5 * MAD)
      calculated per dataset.
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
      Include cells with the percentage of RNA reads mapped to mitochondrial
      genes not bigger than this value. Set to 0 for using an auto-estimated
      threshold equal to the maximum among (median + 2 * MAD) values calculated
      per dataset.
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

  minimum_fragments:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--minfragments"
    doc: |
      Include cells where at least this many ATAC fragments in peaks are
      detected. If multiple values provided, each of them will be
      applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file. Any 0 will be replaced with the
      auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.
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
      to nucleosome-free ATAC fragments. If multiple values provided, each of
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
      Score is calculated based on the ratio of ATAC fragments centered at the TSS
      to ATAC fragments in TSS-flanking regions. If multiple values provided, each
      of them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 2 (applied to all datasets)

  minimum_frip:
    type: float?
    inputBinding:
      prefix: "--minfrip"
    doc: |
      Include cells with the FRiP not lower than this
      value. FRiP is calculated for ATAC fragments.
      Default: 0.15 (applied to all datasets)

  maximum_blacklist_fraction:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--maxblacklist"
    doc: |
      Include cells with the fraction of ATAC fragments in
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

  minimum_qvalue:
    type: float?
    inputBinding:
      prefix: "--qvalue"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2 peak detection.
      Ignored if --callby is not provided. Default: 0.05

  remove_doublets:
    type:
    - "null"
    - type: enum
      symbols:
      - "union"
      - "onlyrna"
      - "onlyatac"
      - "intersect"
    inputBinding:
      prefix: "--removedoublets"
    doc: |
      Remove cells that were identified as doublets. For
      RNA assay cells with UMI < 200 will not be evaluated.
      Default: do not remove doublets

  rna_doublet_rate:
    type: float?
    inputBinding:
      prefix: "--rnadbr"
    doc: |
      Expected RNA doublet rate. Default: 1 percent per
      thousand cells captured with 10x genomics

  rna_doublet_rate_sd:
    type: float?
    inputBinding:
      prefix: "--rnadbrsd"
    doc: |
      Uncertainty range in the RNA doublet rate, interpreted as
      a +/- around the value provided in --rnadbr. Set to 0 to
      disable. Set to 1 to make the threshold depend entirely
      on the misclassification rate. Default: 40 percents of the
      value provided in --rnadbr

  atac_doublet_rate:
    type: float?
    inputBinding:
      prefix: "--atacdbr"
    doc: |
      Expected ATAC doublet rate. Default: 1 percent per thousand
      cells captured with 10x genomics

  atac_doublet_rate_sd:
    type: float?
    inputBinding:
      prefix: "--atacdbrsd"
    doc: |
      Uncertainty range in the ATAC doublet rate, interpreted as
      a +/- around the value provided in --atacdbr. Set to 0 to
      disable. Set to 1 to make the threshold depend entirely
      on the misclassification rate. Default: 40 percents of the
      value provided in --atacdbr

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
      Save raw counts from the RNA and ATAC assays to h5ad files.
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

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered; PC1/PC2.
      PNG format.

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered; PC2/PC3.
      PNG format.

  raw_cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cell_cnts.png"
    doc: |
      Number of cells per dataset.
      Unfiltered.
      PNG format.

  raw_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered.
      PNG format.

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered.
      PNG format.

  raw_gene_umi_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi.png"
    doc: |
      Genes vs RNA reads per cell.
      Unfiltered.
      PNG format.

  raw_umi_mito_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_mito.png"
    doc: |
      RNA reads vs mitochondrial percentage
      per cell.
      Unfiltered.
      PNG format.

  raw_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered.
      PNG format.

  raw_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered.
      PNG format.

  raw_frgm_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_frgm_dnst.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Unfiltered.
      PNG format.

  raw_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst.png"
    doc: |
      Distribution of peaks per cell.
      Unfiltered.
      PNG format.

  raw_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Unfiltered.
      PNG format.

  raw_rna_atac_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rna_atac_cnts.png"
    doc: |
      RNA reads vs ATAC fragments
      in peaks per cell.
      Unfiltered.
      PNG format.

  raw_tss_frgm_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_tss_frgm.png"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Unfiltered.
      PNG format.

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.png"
    doc: |
      Distribution of QC metrics per cell.
      Unfiltered.
      PNG format.

  raw_rnadbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rnadbl.png"
    doc: |
      Percentage of RNA doublets.
      Unfiltered.
      PNG format.

  raw_atacdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_atacdbl.png"
    doc: |
      Percentage of ATAC doublets.
      Unfiltered.
      PNG format.

  raw_vrlpdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_vrlpdbl.png"
    doc: |
      Percentage of RNA and ATAC doublets.
      Unfiltered.
      PNG format.

  raw_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_tss_nrch.png"
    doc: |
      Signal enrichment around TSS.
      Unfiltered; split by the minimum
      TSS enrichment score threshold.
      PNG format.

  raw_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*_raw_frgm_hist.png"
    doc: |
      Histogram of ATAC fragment length.
      Unfiltered; split by the maximum
      nucleosome signal threshold.
      PNG format.

  raw_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_frgm_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_frgm_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_peak_dnst_spl_cnd.png"
    doc: |
      Distribution of peaks per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_blck_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  mid_fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered, after MACS2 peak calling;
      PC1/PC2.
      PNG format.

  mid_fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered, after MACS2 peak calling;
      PC2/PC3.
      PNG format.

  mid_fltr_cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_cell_cnts.png"
    doc: |
      Number of cells per dataset.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_umi_dnst.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_gene_umi_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_umi.png"
    doc: |
      Genes vs RNA reads per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_umi_mito_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_umi_mito.png"
    doc: |
      RNA reads vs mitochondrial percentage
      per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_frgm_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frgm_dnst.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst.png"
    doc: |
      Distribution of peaks per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_rna_atac_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rna_atac_cnts.png"
    doc: |
      RNA reads vs ATAC fragments
      in peaks per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_tss_frgm_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_frgm.png"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_qc_mtrcs_dnst.png"
    doc: |
      Distribution of QC metrics per cell.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_rnadbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_rnadbl.png"
    doc: |
      Percentage of RNA doublets.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_atacdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_atacdbl.png"
    doc: |
      Percentage of ATAC doublets.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_vrlpdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_vrlpdbl.png"
    doc: |
      Percentage of RNA and ATAC doublets.
      Unfiltered, after MACS2 peak calling.
      PNG format.

  mid_fltr_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_tss_nrch.png"
    doc: |
      Signal enrichment around TSS.
      Unfiltered, after MACS2 peak calling;
      split by the minimum TSS enrichment
      score threshold.
      PNG format.

  mid_fltr_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frgm_hist.png"
    doc: |
      Histogram of ATAC fragment length.
      Unfiltered, after MACS2 peak calling;
      split by the maximum nucleosome signal
      threshold.
      PNG format.

  mid_fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_umi_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_frgm_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_frgm_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_peak_dnst_spl_cnd.png"
    doc: |
      Distribution of peaks per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  mid_fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_mid_fltr_blck_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Unfiltered, after MACS2 peak calling;
      split by grouping condition.
      PNG format.

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Filtered; PC1/PC2.
      PNG format.

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Filtered; PC2/PC3.
      PNG format.

  fltr_cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_cell_cnts.png"
    doc: |
      Number of cells per dataset.
      Filtered.
      PNG format.

  fltr_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_umi_dnst.png"
    doc: |
      Distribution of RNA reads per cell.
      Filtered.
      PNG format.

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst.png"
    doc: |
      Distribution of genes per cell.
      Filtered.
      PNG format.

  fltr_gene_umi_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_umi.png"
    doc: |
      Genes vs RNA reads per cell.
      Filtered.
      PNG format.

  fltr_umi_mito_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_umi_mito.png"
    doc: |
      RNA reads vs mitochondrial percentage
      per cell.
      Filtered.
      PNG format.

  fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Filtered.
      PNG format.

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst.png"
    doc: |
      Distribution of novelty score per cell.
      Filtered.
      PNG format.

  fltr_frgm_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frgm_dnst.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Filtered.
      PNG format.

  fltr_peak_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst.png"
    doc: |
      Distribution of peaks per cell.
      Filtered.
      PNG format.

  fltr_blck_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Filtered.
      PNG format.

  fltr_rna_atac_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rna_atac_cnts.png"
    doc: |
      RNA reads vs ATAC fragments
      in peaks per cell.
      Filtered.
      PNG format.

  fltr_tss_frgm_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_frgm.png"
    doc: |
      TSS enrichment score vs ATAC
      fragments in peaks per cell.
      Filtered.
      PNG format.

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_qc_mtrcs_dnst.png"
    doc: |
      Distribution of QC metrics per cell.
      Filtered.
      PNG format.

  fltr_rnadbl_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_rnadbl.png"
    doc: |
      Percentage of RNA doublets.
      Filtered.
      PNG format.

  fltr_atacdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_atacdbl.png"
    doc: |
      Percentage of ATAC doublets.
      Filtered.
      PNG format.

  fltr_vrlpdbl_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_vrlpdbl.png"
    doc: |
      Percentage of RNA and ATAC doublets.
      Filtered.
      PNG format.

  fltr_tss_nrch_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_tss_nrch.png"
    doc: |
      Signal enrichment around TSS.
      Filtered; split by the minimum
      TSS enrichment score threshold.
      PNG format.

  fltr_frgm_hist_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frgm_hist.png"
    doc: |
      Histogram of ATAC fragment length.
      Filtered; split by the maximum
      nucleosome signal threshold.
      PNG format.

  fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_umi_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Distribution of genes per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Distribution of novelty score per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_frgm_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_frgm_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments in peaks
      per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_peak_dnst_spl_cnd.png"
    doc: |
      Distribution of peaks per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*[!_mid]_fltr_blck_dnst_spl_cnd.png"
    doc: |
      Distribution of ATAC fragments within
      genomic blacklist regions per cell.
      Filtered; split by grouping condition.
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
      RDS format

  datasets_metadata:
    type: File
    outputBinding:
      glob: "*_meta.tsv"
    doc: |
      Example of datasets metadata file
      in TSV format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat object.
      h5Seurat format

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

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_multiome_filter.R"]

stdout: sc_multiome_filter_stdout.log
stderr: sc_multiome_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis"
s:name: "Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis"
s:alternateName: "Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis"

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
  Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis

  Removes low-quality cells from the outputs of the “Cell Ranger Count
  (RNA+ATAC)” and “Cell Ranger Aggregate (RNA+ATAC)” pipelines. The
  results of this workflow are used in the “Single-Cell RNA-Seq
  Dimensionality Reduction Analysis” and “Single-Cell ATAC-Seq
  Dimensionality Reduction Analysis” pipelines.


s:about: |
  usage: sc_multiome_filter.R [-h] --mex MEX
                                            [--identity IDENTITY] --fragments
                                            FRAGMENTS --annotations ANNOTATIONS
                                            --seqinfo SEQINFO
                                            [--grouping GROUPING]
                                            [--blacklist BLACKLIST]
                                            [--barcodes BARCODES]
                                            [--rnamincells RNAMINCELLS]
                                            [--mingenes [MINGENES [MINGENES ...]]]
                                            [--maxgenes [MAXGENES [MAXGENES ...]]]
                                            [--minumis [MINUMIS [MINUMIS ...]]]
                                            [--mitopattern MITOPATTERN]
                                            [--maxmt MAXMT]
                                            [--minnovelty [MINNOVELTY [MINNOVELTY ...]]]
                                            [--atacmincells ATACMINCELLS]
                                            [--minfragments [MINFRAGMENTS [MINFRAGMENTS ...]]]
                                            [--maxnuclsignal [MAXNUCLSIGNAL [MAXNUCLSIGNAL ...]]]
                                            [--mintssenrich [MINTSSENRICH [MINTSSENRICH ...]]]
                                            [--minfrip MINFRIP]
                                            [--maxblacklist [MAXBLACKLIST [MAXBLACKLIST ...]]]
                                            [--callby CALLBY] [--qvalue QVALUE]
                                            [--removedoublets {union,onlyrna,onlyatac,intersect}]
                                            [--rnadbr RNADBR]
                                            [--rnadbrsd RNADBRSD]
                                            [--atacdbr ATACDBR]
                                            [--atacdbrsd ATACDBRSD] [--pdf]
                                            [--verbose] [--h5seurat] [--h5ad]
                                            [--loupe] [--cbbuild]
                                            [--tmpdir TMPDIR] [--output OUTPUT]
                                            [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                            [--cpus CPUS] [--memory MEMORY]
                                            [--seed SEED]

  Single-Cell Multiome ATAC-Seq and RNA-Seq Filtering Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX             Path to the folder with feature-barcode matrix from
                          Cell Ranger ARC Count/Aggregate experiment in MEX
                          format. The rows consist of all the genes and peaks
                          concatenated together and the columns are restricted
                          to those barcodes that are identified as cells.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities, if '--mex' points to the Cell Ranger ARC
                          Aggregate outputs. The aggr.csv file can be used.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          observed in the experiment in TSV format. Tbi-index
                          file is required.
    --annotations ANNOTATIONS
                          Path to the genome annotation file in GTF format
    --seqinfo SEQINFO     Path to the headerless chromosome length file in TSV
                          format
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
    --mingenes [MINGENES [MINGENES ...]]
                          Include cells where at least this many genes are
                          detected. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median - 2.5 * MAD) calculated per dataset. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES [MAXGENES ...]]
                          Include cells with the number of genes not bigger than
                          this value. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median + 5 * MAD) calculated per dataset. Default:
                          5000 (applied to all datasets)
    --minumis [MINUMIS [MINUMIS ...]]
                          Include cells where at least this many RNA reads are
                          detected. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median - 2.5 * MAD) calculated per dataset. Default:
                          500 (applied to all datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial genes.
                          Default: '^mt-|^MT-'
    --maxmt MAXMT         Include cells with the percentage of RNA reads mapped
                          to mitochondrial genes not bigger than this value. Set
                          to 0 for using an auto-estimated threshold equal to
                          the maximum among (median + 2 * MAD) values calculated
                          per dataset. Default: 5 (applied to all datasets)
    --minnovelty [MINNOVELTY [MINNOVELTY ...]]
                          Include cells with the novelty score not lower than
                          this value, calculated for as log10(genes)/log10(UMI)
                          for RNA assay. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the '--mex' input based on the '--identity' file.
                          Default: 0.8 (applied to all datasets)
    --atacmincells ATACMINCELLS
                          Include only peaks detected in at least this many
                          cells. Default: 5 (applied to all datasets)
    --minfragments [MINFRAGMENTS [MINFRAGMENTS ...]]
                          Include cells where at least this many ATAC fragments
                          in peaks are detected. If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the '--mex' input based on the '--
                          identity' file. Any 0 will be replaced with the auto-
                          estimated threshold (median - 2.5 * MAD) calculated
                          per dataset. Default: 1000 (applied to all datasets)
    --maxnuclsignal [MAXNUCLSIGNAL [MAXNUCLSIGNAL ...]]
                          Include cells with the nucleosome signal not bigger
                          than this value. Nucleosome signal quantifies the
                          approximate ratio of mononucleosomal to nucleosome-
                          free ATAC fragments. If multiple values provided, each
                          of them will be applied to the correspondent dataset
                          from the '--mex' input based on the '--identity' file.
                          Default: 4 (applied to all datasets)
    --mintssenrich [MINTSSENRICH [MINTSSENRICH ...]]
                          Include cells with the TSS enrichment score not lower
                          than this value. Score is calculated based on the
                          ratio of ATAC fragments centered at the TSS to ATAC
                          fragments in TSS-flanking regions. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the '--mex' input based on
                          the '--identity' file. Default: 2 (applied to all
                          datasets)
    --minfrip MINFRIP     Include cells with the FRiP not lower than this value.
                          FRiP is calculated for ATAC fragments. Default: 0.15
                          (applied to all datasets)
    --maxblacklist [MAXBLACKLIST [MAXBLACKLIST ...]]
                          Include cells with the fraction of ATAC fragments in
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
    --qvalue QVALUE       Minimum FDR (q-value) cutoff for MACS2 peak detection.
                          Ignored if --callby is not provided. Default: 0.05
    --removedoublets {union,onlyrna,onlyatac,intersect}
                          Remove cells that were identified as doublets. For RNA
                          assay cells with UMI < 200 will not be evaluated.
                          Default: do not remove doublets
    --rnadbr RNADBR       Expected RNA doublet rate. Default: 1 percent per
                          thousand cells captured with 10x genomics
    --rnadbrsd RNADBRSD   Uncertainty range in the RNA doublet rate, interpreted
                          as a +/- around the value provided in --rnadbr. Set to
                          0 to disable. Set to 1 to make the threshold depend
                          entirely on the misclassification rate. Default: 40
                          percents of the value provided in --rnadbr
    --atacdbr ATACDBR     Expected ATAC doublet rate. Default: 1 percent per
                          thousand cells captured with 10x genomics
    --atacdbrsd ATACDBRSD
                          Uncertainty range in the ATAC doublet rate,
                          interpreted as a +/- around the value provided in
                          --atacdbr. Set to 0 to disable. Set to 1 to make the
                          threshold depend entirely on the misclassification
                          rate. Default: 40 percents of the value provided in
                          --atacdbr
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA and ATAC assays to h5ad
                          files. Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --tmpdir TMPDIR       Directory to keep temporary files. Default: either
                          /tmp or defined by environment variables TMPDIR, TMP,
                          TEMP.
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple '--cpus'. Default: 32
    --seed SEED           Seed number for random values. Default: 42