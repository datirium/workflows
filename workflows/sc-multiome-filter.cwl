cwlVersion: v1.1
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
  sc_arc_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-count.cwl"
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-aggr.cwl"
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  filtered_feature_bc_matrix_folder:
    type: File
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Path to the compressed folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate
      experiment in MEX format. The rows consist of all the genes and peaks concatenated
      together and the columns are restricted to those barcodes that are identified as cells.
    'sd:upstreamSource': "sc_arc_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File?
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger ARC Aggregate outputs, the aggr.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.
    'sd:upstreamSource': "sc_arc_sample/aggregation_metadata"

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Count and barcode information for every ATAC fragment observed in the experiment in TSV
      format. Tbi-index file is required.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"

  annotation_gtf_file:
    type: File
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Path to the genome annotation file in GTF format.
    'sd:upstreamSource': "sc_arc_sample/genome_indices/genome_indices/annotation_gtf"
    'sd:localLabel': true

  grouping_data:
    type: File?
    label: "Optional TSV/CSV file to define datasets grouping with 'library_id' and 'condition' columns. Rows order should correspond to the aggregation metadata."
    doc: |
      Path to the TSV/CSV file to define datasets grouping.
      First column - 'library_id' with the values and order
      that correspond to the 'library_id' column from the '
      --identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  blacklist_regions_file:
    type: File?
    label: "Optional BED file with the genomic blacklist regions"
    doc: |
      Path to the optional BED file with the genomic blacklist regions.

  barcodes_data:
    type: File?
    label: "Optional TSV/CSV file to prefilter and extend metadata be barcodes. First column should be named as 'barcode'"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

  minimum_genes:
    type: string?
    default: "250"
    label: "Include cells where at least this many genes are detected"
    doc: |
      Include cells where at least this many genes are detected. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file.
      Default: 250 (applied to all datasets)
    'sd:layout':
      advanced: true

  maximum_genes:
    type: string?
    default: "5000"
    label: "Include cells with the number of genes not bigger than this value"
    doc: |
      Include cells with the number of genes not bigger than this value. If multiple
      values provided, each of them will be applied to the correspondent dataset from
      the '--mex' input based on the '--identity' file.
      Default: 5000 (applied to all datasets)
    'sd:layout':
      advanced: true

  rna_minimum_umi:
    type: string?
    default: "500"
    label: "Include cells where at least this many UMI (RNA transcripts) are detected"
    doc: |
      Include cells where at least this many UMI (RNA transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 500 (applied to all datasets)
    'sd:layout':
      advanced: true

  mito_pattern:
    type: string?
    default: "^mt-|^MT-"
    label: "Regex pattern to identify mitochondrial genes"
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^mt-|^MT-'
    'sd:layout':
      advanced: true

  maximum_mito_perc:
    type: float?
    default: 5
    label: "Include cells with the percentage of transcripts mapped to mitochondrial genes not bigger than this value"
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5 (applied to all datasets)
    'sd:layout':
      advanced: true

  minimum_novelty_score:
    type: string?
    default: "0.8"
    label: "Include cells with the novelty score not lower than this value, calculated as log10(genes)/log10(UMI) for RNA assay"
    doc: |
      Include cells with the novelty score not lower than this value, calculated
      as log10(genes)/log10(UMI) for RNA assay. If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)
    'sd:layout':
      advanced: true

  atac_minimum_umi:
    type: string?
    default: "1000"
    label: "Include cells where at least this many UMI (ATAC transcripts) are detected"
    doc: |
      Include cells where at least this many UMI (ATAC transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 1000 (applied to all datasets)
    'sd:layout':
      advanced: true

  maximum_nucl_signal:
    type: string?
    default: "4"
    label: "Include cells with the nucleosome signal not bigger than this value"
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments. If multiple values provided, each of
      them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 4 (applied to all datasets)
    'sd:layout':
      advanced: true

  minimum_tss_enrich:
    type: string?
    default: "2"
    label: "Include cells with the TSS enrichment score not lower than this value"
    doc: |
      Include cells with the TSS enrichment score not lower than this value.
      Score is calculated based on the ratio of fragments centered at the TSS
      to fragments in TSS-flanking regions. If multiple values provided, each
      of them will be applied to the correspondent dataset from the '--mex' input
      based on the '--identity' file.
      Default: 2 (applied to all datasets)
    'sd:layout':
      advanced: true

  minimum_frip:
    type: string?
    default: "0.15"
    label: "Include cells with the FRiP not lower than this value"
    doc: |
      Include cells with the FRiP not lower than this value. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file. FRiP is calculated for fragments.
      Default: 0.15 (applied to all datasets)
    'sd:layout':
      advanced: true

  maximum_blacklist_fraction:
    type: string?
    default: "0.05"
    label: "Include cells with the fraction of fragments in genomic blacklist regions not bigger than this value"
    doc: |
      Include cells with the fraction of fragments in
      genomic blacklist regions not bigger than this value.
      If multiple values provided, each of them will be
      applied to the correspondent dataset from the '--mex'
      input based on the '--identity' file.
      Default: 0.05 (applied to all datasets)
    'sd:layout':
      advanced: true

  call_by:
    type: string?
    default: null
    label: "Replace Cell Ranger ARC peaks with MACS2 peaks called for cells grouped by selected column"
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
    'sd:layout':
      advanced: true

  remove_doublets:
    type: boolean?
    default: false
    label: "Remove cells that were identified as doublets in either RNA or ATAC assays"
    doc: |
      Remove cells that were identified as doublets in either
      RNA or ATAC assays.
      Default: do not remove doublets
    'sd:layout':
      advanced: true

  rna_doublet_rate:
    type: float?
    default: null
    label: "Expected RNA doublet rate"
    doc: |
      Expected RNA doublet rate. Default: 1 percent per
      thousand cells captured with 10x genomics
    'sd:layout':
      advanced: true

  rna_doublet_rate_sd:
    type: float?
    default: null
    label: "Uncertainty range in the RNA doublet rate"
    doc: |
      Uncertainty range in the RNA doublet rate, interpreted as
      a +/- around the value provided in --rnadbr. Set to 0 to
      disable. Set to 1 to make the threshold depend entirely
      on the misclassification rate. Default: 40 percents of the
      value provided in --rnadbr
    'sd:layout':
      advanced: true

  atac_doublet_rate:
    type: float?
    default: null
    label: "Expected ATAC doublet rate"
    doc: |
      Expected ATAC doublet rate. Default: 1 percent per thousand
      cells captured with 10x genomics
    'sd:layout':
      advanced: true

  atac_doublet_rate_sd:
    type: float?
    default: null
    label: "Uncertainty range in the ATAC doublet rate"
    doc: |
      Uncertainty range in the ATAC doublet rate, interpreted as
      a +/- around the value provided in --atacdbr. Set to 0 to
      disable. Set to 1 to make the threshold depend entirely
      on the misclassification rate. Default: 40 percents of the
      value provided in --atacdbr
    'sd:layout':
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
    label: "Color theme for all generated plots"
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic
    'sd:layout':
      advanced: true

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
      - "96"
    default: "96"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Forced to 96 GB
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
      Number of cores/cpus to use.
      Forced to 1
    'sd:layout':
      advanced: true


outputs:


  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_1_2_qc_mtrcs_pca_plot_png
    label: "PC1 and PC2 from the QC metrics PCA (not filtered)"
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'PC1 and PC2 from the QC metrics PCA'

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_2_3_qc_mtrcs_pca_plot_png
    label: "PC2 and PC3 from the QC metrics PCA (not filtered)"
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'PC2 and PC3 from the QC metrics PCA'

  raw_cells_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_cells_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Number of cells per dataset'

  raw_rna_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_umi_dnst_plot_png
    label: "UMI per cell density for RNA assay (not filtered)"
    doc: |
      UMI per cell density for RNA assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'UMI per cell density for RNA assay'

  raw_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_dnst_plot_png
    label: "Genes per cell density (not filtered)"
    doc: |
      Genes per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Genes per cell density'

  raw_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_umi_corr_plot_png
    label: "Genes vs UMI per cell correlation for RNA assay (not filtered)"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Genes vs UMI per cell correlation for RNA assay'

  raw_mito_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_mito_dnst_plot_png
    label: "Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered)"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Percentage of transcripts mapped to mitochondrial genes per cell density'

  raw_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_nvlt_dnst_plot_png
    label: "Novelty score per cell density for RNA assay (not filtered)"
    doc: |
      Novelty score per cell density for RNA assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Novelty score per cell density for RNA assay'

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atac_umi_dnst_plot_png
    label: "UMI per cell density for ATAC assay (not filtered)"
    doc: |
      UMI per cell density for ATAC assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'UMI per cell density for ATAC assay'

  raw_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_peak_dnst_plot_png
    label: "Peaks per cell density (not filtered)"
    doc: |
      Peaks per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Peaks per cell density'

  raw_blck_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_blck_dnst_plot_png
    label: "Fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered)"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Fraction of ATAC fragments within genomic blacklist regions per cell density'

  raw_rna_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_atac_umi_corr_plot_png
    label: "UMI per cell correlation for RNA vs ATAC assays (not filtered)"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'UMI per cell correlation for RNA vs ATAC assays'

  raw_tss_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_tss_atac_umi_corr_plot_png
    label: "TSS enrichment score vs UMI per cell correlation for ATAC assay (not filtered)"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'TSS enrichment score vs UMI per cell correlation for ATAC assay'

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_qc_mtrcs_dnst_plot_png
    label: "QC metrics per cell density (not filtered)"
    doc: |
      QC metrics per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'QC metrics per cell density'

  raw_rnadbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rnadbl_plot_png
    label: "Percentage of RNA doublets per dataset (not filtered)"
    doc: |
      Percentage of RNA doublets per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Percentage of RNA doublets per dataset'

  raw_atacdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atacdbl_plot_png
    label: "Percentage of ATAC doublets per dataset (not filtered)"
    doc: |
      Percentage of ATAC doublets per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Percentage of ATAC doublets per dataset'

  raw_vrlpdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_vrlpdbl_plot_png
    label: "Doublets overlap for RNA and ATAC assays per dataset (not filtered)"
    doc: |
      Doublets overlap for RNA and ATAC assays per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Doublets overlap for RNA and ATAC assays per dataset'

  raw_tss_nrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_tss_nrch_plot_png
    label: "TSS enrichment score (not filtered)"
    doc: |
      TSS enrichment score (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'TSS enrichment score'

  raw_frgm_hist_png:
    type: File?
    outputSource: sc_multiome_filter/raw_frgm_hist_png
    label: "Fragments length histogram (not filtered)"
    doc: |
      Fragments length histogram (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Fragments length histogram'

  raw_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_rna_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for RNA assay (not filtered)"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for RNA assay'

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_dnst_spl_cnd_plot_png
    label: "Split by grouping condition genes per cell density (not filtered)"
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition genes per cell density'

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_mito_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density (not filtered)"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density'

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_nvlt_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the novelty score per cell density for RNA assay (not filtered)"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition the novelty score per cell density for RNA assay'

  raw_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atac_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for ATAC assay (not filtered)"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for ATAC assay'

  raw_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_peak_dnst_spl_cnd_plot_png
    label: "Split by grouping condition peaks per cell density (not filtered)"
    doc: |
      Split by grouping condition peaks per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition peaks per cell density'

  raw_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_blck_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density (not filtered)"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density'

  mid_fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_1_2_qc_mtrcs_pca_plot_png
    label: "PC1 and PC2 from the QC metrics PCA (intermediate filtered)"
    doc: |
      PC1 and PC2 from the QC metrics PCA (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'PC1 and PC2 from the QC metrics PCA'

  mid_fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_2_3_qc_mtrcs_pca_plot_png
    label: "PC2 and PC3 from the QC metrics PCA (intermediate filtered)"
    doc: |
      PC2 and PC3 from the QC metrics PCA (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'PC2 and PC3 from the QC metrics PCA'

  mid_fltr_cells_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_cells_count_plot_png
    label: "Number of cells per dataset (intermediate filtered)"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Number of cells per dataset'

  mid_fltr_rna_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_rna_umi_dnst_plot_png
    label: "UMI per cell density for RNA assay (intermediate filtered)"
    doc: |
      UMI per cell density for RNA assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'UMI per cell density for RNA assay'

  mid_fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gene_dnst_plot_png
    label: "Genes per cell density (intermediate filtered)"
    doc: |
      Genes per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Genes per cell density'

  mid_fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gene_umi_corr_plot_png
    label: "Genes vs UMI per cell correlation for RNA assay (intermediate filtered)"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Genes vs UMI per cell correlation for RNA assay'

  mid_fltr_mito_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_mito_dnst_plot_png
    label: "Percentage of transcripts mapped to mitochondrial genes per cell density (intermediate filtered)"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Percentage of transcripts mapped to mitochondrial genes per cell density'

  mid_fltr_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_nvlt_dnst_plot_png
    label: "Novelty score per cell density for RNA assay (intermediate filtered)"
    doc: |
      Novelty score per cell density for RNA assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Novelty score per cell density for RNA assay'

  mid_fltr_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_atac_umi_dnst_plot_png
    label: "UMI per cell density for ATAC assay (intermediate filtered)"
    doc: |
      UMI per cell density for ATAC assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'UMI per cell density for ATAC assay'

  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_peak_dnst_plot_png
    label: "Peaks per cell density (intermediate filtered)"
    doc: |
      Peaks per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Peaks per cell density'

  mid_fltr_blck_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_blck_dnst_plot_png
    label: "Fraction of ATAC fragments within genomic blacklist regions per cell density (intermediate filtered)"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Fraction of ATAC fragments within genomic blacklist regions per cell density'

  mid_fltr_rna_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_rna_atac_umi_corr_plot_png
    label: "UMI per cell correlation for RNA vs ATAC assays (intermediate filtered)"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'UMI per cell correlation for RNA vs ATAC assays'

  mid_fltr_tss_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_tss_atac_umi_corr_plot_png
    label: "TSS enrichment score vs UMI per cell correlation for ATAC assay (intermediate filtered)"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'TSS enrichment score vs UMI per cell correlation for ATAC assay'

  mid_fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_qc_mtrcs_dnst_plot_png
    label: "QC metrics per cell density (intermediate filtered)"
    doc: |
      QC metrics per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'QC metrics per cell density'

  mid_fltr_rnadbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_rnadbl_plot_png
    label: "Percentage of RNA doublets per dataset (intermediate filtered)"
    doc: |
      Percentage of RNA doublets per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Percentage of RNA doublets per dataset'

  mid_fltr_atacdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_atacdbl_plot_png
    label: "Percentage of ATAC doublets per dataset (intermediate filtered)"
    doc: |
      Percentage of ATAC doublets per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Percentage of ATAC doublets per dataset'

  mid_fltr_vrlpdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_vrlpdbl_plot_png
    label: "Doublets overlap for RNA and ATAC assays per dataset (intermediate filtered)"
    doc: |
      Doublets overlap for RNA and ATAC assays per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Doublets overlap for RNA and ATAC assays per dataset'

  mid_fltr_tss_nrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_tss_nrch_plot_png
    label: "TSS enrichment score (intermediate filtered)"
    doc: |
      TSS enrichment score (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'TSS enrichment score'

  mid_fltr_frgm_hist_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_frgm_hist_png
    label: "Fragments length histogram (intermediate filtered)"
    doc: |
      Fragments length histogram (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Fragments length histogram'

  mid_fltr_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_rna_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for RNA assay (intermediate filtered)"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for RNA assay'

  mid_fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gene_dnst_spl_cnd_plot_png
    label: "Split by grouping condition genes per cell density (intermediate filtered)"
    doc: |
      Split by grouping condition genes per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition genes per cell density'

  mid_fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_mito_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density (intermediate filtered)"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density'

  mid_fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_nvlt_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the novelty score per cell density for RNA assay (intermediate filtered)"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition the novelty score per cell density for RNA assay'

  mid_fltr_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_atac_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for ATAC assay (intermediate filtered)"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for ATAC assay'

  mid_fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_peak_dnst_spl_cnd_plot_png
    label: "Split by grouping condition peaks per cell density (intermediate filtered)"
    doc: |
      Split by grouping condition peaks per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition peaks per cell density'

  mid_fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_blck_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density (intermediate filtered)"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Mid. filtered QC'
        Caption: 'Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density'

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_1_2_qc_mtrcs_pca_plot_png
    label: "PC1 and PC2 from the QC metrics PCA (filtered)"
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'PC1 and PC2 from the QC metrics PCA'

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_2_3_qc_mtrcs_pca_plot_png
    label: "PC2 and PC3 from the QC metrics PCA (filtered)"
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'PC2 and PC3 from the QC metrics PCA'

  fltr_cells_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_cells_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Number of cells per dataset'

  fltr_rna_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_umi_dnst_plot_png
    label: "UMI per cell density for RNA assay (filtered)"
    doc: |
      UMI per cell density for RNA assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'UMI per cell density for RNA assay'

  fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_dnst_plot_png
    label: "Genes per cell density (filtered)"
    doc: |
      Genes per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Genes per cell density'

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_umi_corr_plot_png
    label: "Genes vs UMI per cell correlation for RNA assay (filtered)"
    doc: |
      Genes vs UMI per cell correlation for RNA assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Genes vs UMI per cell correlation for RNA assay'

  fltr_mito_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_mito_dnst_plot_png
    label: "Percentage of transcripts mapped to mitochondrial genes per cell density (filtered)"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Percentage of transcripts mapped to mitochondrial genes per cell density'

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_nvlt_dnst_plot_png
    label: "Novelty score per cell density for RNA assay (filtered)"
    doc: |
      Novelty score per cell density for RNA assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Novelty score per cell density for RNA assay'

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atac_umi_dnst_plot_png
    label: "UMI per cell density for ATAC assay (filtered)"
    doc: |
      UMI per cell density for ATAC assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'UMI per cell density for ATAC assay'

  fltr_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_peak_dnst_plot_png
    label: "Peaks per cell density (filtered)"
    doc: |
      Peaks per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Peaks per cell density'

  fltr_blck_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_blck_dnst_plot_png
    label: "Fraction of ATAC fragments within genomic blacklist regions per cell density (filtered)"
    doc: |
      Fraction of ATAC fragments within genomic blacklist regions per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Fraction of ATAC fragments within genomic blacklist regions per cell density'

  fltr_rna_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_atac_umi_corr_plot_png
    label: "UMI per cell correlation for RNA vs ATAC assays (filtered)"
    doc: |
      UMI per cell correlation for RNA vs ATAC assays (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'UMI per cell correlation for RNA vs ATAC assays'

  fltr_rnadbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rnadbl_plot_png
    label: "Percentage of RNA doublets per dataset (filtered)"
    doc: |
      Percentage of RNA doublets per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Percentage of RNA doublets per dataset'

  fltr_atacdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atacdbl_plot_png
    label: "Percentage of ATAC doublets per dataset (filtered)"
    doc: |
      Percentage of ATAC doublets per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Percentage of ATAC doublets per dataset'

  fltr_vrlpdbl_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_vrlpdbl_plot_png
    label: "Doublets overlap for RNA and ATAC assays per dataset (filtered)"
    doc: |
      Doublets overlap for RNA and ATAC assays per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Doublets overlap for RNA and ATAC assays per dataset'

  fltr_tss_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_tss_atac_umi_corr_plot_png
    label: "TSS enrichment score vs UMI per cell correlation for ATAC assay (filtered)"
    doc: |
      TSS enrichment score vs UMI per cell correlation for ATAC assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'TSS enrichment score vs UMI per cell correlation for ATAC assay'

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_qc_mtrcs_dnst_plot_png
    label: "QC metrics per cell density (filtered)"
    doc: |
      QC metrics per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'QC metrics per cell density'

  fltr_tss_nrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_tss_nrch_plot_png
    label: "TSS enrichment score (filtered)"
    doc: |
      TSS enrichment score (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'TSS enrichment score'

  fltr_frgm_hist_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_frgm_hist_png
    label: "Fragments length histogram (filtered)"
    doc: |
      Fragments length histogram (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Fragments length histogram'

  fltr_rna_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_rna_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for RNA assay (filtered)"
    doc: |
      Split by grouping condition UMI per cell density for RNA assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for RNA assay'

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_dnst_spl_cnd_plot_png
    label: "Split by grouping condition genes per cell density (filtered)"
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition genes per cell density'

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_mito_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density (filtered)"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density'

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_nvlt_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the novelty score per cell density for RNA assay (filtered)"
    doc: |
      Split by grouping condition the novelty score per cell density for RNA assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition the novelty score per cell density for RNA assay'

  fltr_atac_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atac_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density for ATAC assay (filtered)"
    doc: |
      Split by grouping condition UMI per cell density for ATAC assay (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition UMI per cell density for ATAC assay'

  fltr_peak_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_peak_dnst_spl_cnd_plot_png
    label: "Split by grouping condition peaks per cell density (filtered)"
    doc: |
      Split by grouping condition peaks per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition peaks per cell density'

  fltr_blck_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_blck_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density (filtered)"
    doc: |
      Split by grouping condition the fraction of ATAC fragments within genomic
      blacklist regions per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition the fraction of ATAC fragments within genomic blacklist regions per cell density'

  ucsc_cb_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data.

  ucsc_cb_html_data:
    type: Directory
    outputSource: sc_multiome_filter/ucsc_cb_html_data
    label: "Directory with UCSC Cellbrowser html data"
    doc: |
      Directory with UCSC Cellbrowser html data.

  ucsc_cb_html_file:
    type: File
    outputSource: sc_multiome_filter/ucsc_cb_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser html data.
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: sc_multiome_filter/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  sc_multiome_filter_stdout_log:
    type: File
    outputSource: sc_multiome_filter/stdout_log
    label: "stdout log generated by sc_multiome_filter step"
    doc: |
      stdout log generated by sc_multiome_filter step

  sc_multiome_filter_stderr_log:
    type: File
    outputSource: sc_multiome_filter/stderr_log
    label: "stderr log generated by sc_multiome_filter step"
    doc: |
      stderr log generated by sc_multiome_filter step


steps:

  uncompress_feature_bc_matrices:
    doc: |
      Extracts the content of TAR file into a folder
    run: ../tools/tar-extract.cwl
    in:
      file_to_extract: filtered_feature_bc_matrix_folder
    out:
    - extracted_folder

  sc_multiome_filter:
    doc: |
      Filters single-cell multiome ATAC and RNA-Seq datasets
      based on the common QC metrics
    run: ../tools/sc-multiome-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/extracted_folder
      aggregation_metadata: aggregation_metadata
      atac_fragments_file: atac_fragments_file
      annotation_gtf_file: annotation_gtf_file
      grouping_data: grouping_data
      blacklist_regions_file: blacklist_regions_file
      barcodes_data: barcodes_data
      rna_minimum_cells:
        default: 1
      minimum_genes:
        source: minimum_genes
        valueFrom: $(split_numbers(self))
      maximum_genes:
        source: maximum_genes
        valueFrom: $(split_numbers(self))
      rna_minimum_umi:
        source: rna_minimum_umi
        valueFrom: $(split_numbers(self))
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      atac_minimum_cells:
        default: 1
      atac_minimum_umi:
        source: atac_minimum_umi
        valueFrom: $(split_numbers(self))
      maximum_nucl_signal:
        source: maximum_nucl_signal
        valueFrom: $(split_numbers(self))
      minimum_tss_enrich:
        source: minimum_tss_enrich
        valueFrom: $(split_numbers(self))
      minimum_frip:
        source: minimum_frip
        valueFrom: $(split_numbers(self))
      maximum_blacklist_fraction:
        source: maximum_blacklist_fraction
        valueFrom: $(split_numbers(self))
      call_by: call_by
      remove_doublets: remove_doublets
      rna_doublet_rate:
        source: rna_doublet_rate
        valueFrom: $(self==""?null:self)                 # safety measure
      rna_doublet_rate_sd:
        source: rna_doublet_rate_sd
        valueFrom: $(self==""?null:self)                 # safety measure
      atac_doublet_rate:
        source: atac_doublet_rate
        valueFrom: $(self==""?null:self)                 # safety measure
      atac_doublet_rate_sd:
        source: atac_doublet_rate_sd
        valueFrom: $(self==""?null:self)                 # safety measure
      verbose:
        default: true
      export_ucsc_cb:
        default: true
      color_theme: color_theme
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
    - raw_1_2_qc_mtrcs_pca_plot_png
    - raw_2_3_qc_mtrcs_pca_plot_png
    - raw_cells_count_plot_png
    - raw_rna_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_dnst_plot_png
    - raw_nvlt_dnst_plot_png
    - raw_atac_umi_dnst_plot_png
    - raw_peak_dnst_plot_png
    - raw_blck_dnst_plot_png
    - raw_rna_atac_umi_corr_plot_png
    - raw_tss_atac_umi_corr_plot_png
    - raw_qc_mtrcs_dnst_plot_png
    - raw_rnadbl_plot_png
    - raw_atacdbl_plot_png
    - raw_vrlpdbl_plot_png
    - raw_tss_nrch_plot_png
    - raw_frgm_hist_png
    - raw_rna_umi_dnst_spl_cnd_plot_png
    - raw_gene_dnst_spl_cnd_plot_png
    - raw_mito_dnst_spl_cnd_plot_png
    - raw_nvlt_dnst_spl_cnd_plot_png
    - raw_atac_umi_dnst_spl_cnd_plot_png
    - raw_peak_dnst_spl_cnd_plot_png
    - raw_blck_dnst_spl_cnd_plot_png
    - mid_fltr_1_2_qc_mtrcs_pca_plot_png
    - mid_fltr_2_3_qc_mtrcs_pca_plot_png
    - mid_fltr_cells_count_plot_png
    - mid_fltr_rna_umi_dnst_plot_png
    - mid_fltr_gene_dnst_plot_png
    - mid_fltr_gene_umi_corr_plot_png
    - mid_fltr_mito_dnst_plot_png
    - mid_fltr_nvlt_dnst_plot_png
    - mid_fltr_atac_umi_dnst_plot_png
    - mid_fltr_peak_dnst_plot_png
    - mid_fltr_blck_dnst_plot_png
    - mid_fltr_rna_atac_umi_corr_plot_png
    - mid_fltr_tss_atac_umi_corr_plot_png
    - mid_fltr_qc_mtrcs_dnst_plot_png
    - mid_fltr_rnadbl_plot_png
    - mid_fltr_atacdbl_plot_png
    - mid_fltr_vrlpdbl_plot_png
    - mid_fltr_tss_nrch_plot_png
    - mid_fltr_frgm_hist_png
    - mid_fltr_rna_umi_dnst_spl_cnd_plot_png
    - mid_fltr_gene_dnst_spl_cnd_plot_png
    - mid_fltr_mito_dnst_spl_cnd_plot_png
    - mid_fltr_nvlt_dnst_spl_cnd_plot_png
    - mid_fltr_atac_umi_dnst_spl_cnd_plot_png
    - mid_fltr_peak_dnst_spl_cnd_plot_png
    - mid_fltr_blck_dnst_spl_cnd_plot_png
    - fltr_1_2_qc_mtrcs_pca_plot_png
    - fltr_2_3_qc_mtrcs_pca_plot_png
    - fltr_cells_count_plot_png
    - fltr_rna_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_dnst_plot_png
    - fltr_nvlt_dnst_plot_png
    - fltr_atac_umi_dnst_plot_png
    - fltr_peak_dnst_plot_png
    - fltr_blck_dnst_plot_png
    - fltr_rna_atac_umi_corr_plot_png
    - fltr_rnadbl_plot_png
    - fltr_atacdbl_plot_png
    - fltr_vrlpdbl_plot_png
    - fltr_tss_atac_umi_corr_plot_png
    - fltr_qc_mtrcs_dnst_plot_png
    - fltr_tss_nrch_plot_png
    - fltr_frgm_hist_png
    - fltr_rna_umi_dnst_spl_cnd_plot_png
    - fltr_gene_dnst_spl_cnd_plot_png
    - fltr_mito_dnst_spl_cnd_plot_png
    - fltr_nvlt_dnst_spl_cnd_plot_png
    - fltr_atac_umi_dnst_spl_cnd_plot_png
    - fltr_peak_dnst_spl_cnd_plot_png
    - fltr_blck_dnst_spl_cnd_plot_png
    - ucsc_cb_config_data
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: sc_multiome_filter/ucsc_cb_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Multiome ATAC and RNA-Seq Filtering Analysis"
s:name: "Single-cell Multiome ATAC and RNA-Seq Filtering Analysis"
s:alternateName: "Filters single-cell multiome ATAC and RNA-Seq datasets based on the common QC metrics"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-multiome-filter.cwl
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
  Single-cell Multiome ATAC and RNA-Seq Filtering Analysis
  
  Filters single-cell multiome ATAC and RNA-Seq datasets
  based on the common QC metrics.