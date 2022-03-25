cwlVersion: v1.1
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
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
      Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate
      experiment in MEX format. The rows consist of all the gene and peak features
      concatenated together and the columns are restricted to those barcodes that are
      identified as cells.
    'sd:upstreamSource': "sc_arc_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File?
    loadContents: true
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities.
      If --mex points to the Cell Ranger ARC Aggregate outputs, the aggr.csv
      file can be used. If Cell Ranger ARC Count outputs have been used in
      --mex, the file should include at least one column - 'library_id' and
      one row with the alias for Cell Ranger ARC Count experiment.
    'sd:upstreamSource': "sc_arc_sample/aggregation_metadata"
    'sd:localLabel': true

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format. Tbi-index file is required.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"
    'sd:localLabel': true

  annotation_gtf_file:
    type: File
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Path to the genome annotation file in GTF format
    'sd:upstreamSource': "sc_arc_sample/genome_indices/genome_indices/annotation_gtf"
    'sd:localLabel': true

  grouping_data:
    type: File?
    label: "TSV/CSV file to define datasets grouping with 'library_id' and 'condition' columns. Rows order should correspond to the aggregation metadata."
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      Default: each dataset is assigned to a separate group.

  blacklisted_regions_file:
    type: File?
    label: "BED file with blacklisted regions"
    doc: |
      Path to the blacklisted regions file in BED format

  barcodes_data:
    type: File?
    label: "Headerless TSV/CSV file with the list of barcodes to select cells of interest (one barcode per line)"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells.
      Default: use all cells.

  gex_minimum_cells:
    type: int?
    default: 5
    label: "Include only GEX features detected in at least this many cells"
    doc: |
      Include only GEX features detected in at least this many cells.
    'sd:layout':
      advanced: true

  gex_minimum_features:
    type: string?
    default: "250"
    label: "Include cells where at least this many GEX features are detected"
    doc: |
      Include cells where at least this many GEX features are detected.
      If multiple values provided, each of them will be applied to the
      correspondent dataset from the --mex input based on the --identity
      file.
    'sd:layout':
      advanced: true

  gex_maximum_features:
    type: string?
    default: "5000"
    label: "Include cells with the number of GEX features not bigger than this value"
    doc: |
      Include cells with the number of GEX features not bigger than this value.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 5000 (applied to all datasets)
    'sd:layout':
      advanced: true

  gex_minimum_umis:
    type: string?
    default: "500"
    label: "Include cells where at least this many GEX UMIs (transcripts) are detected"
    doc: |
      Include cells where at least this many GEX UMIs (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 500 (applied to all datasets)
    'sd:layout':
      advanced: true

  mito_pattern:
    type: string?
    default: "^Mt-"
    label: "Regex pattern to identify mitochondrial GEX features"
    doc: |
      Regex pattern to identify mitochondrial GEX features.
    'sd:layout':
      advanced: true

  maximum_mito_perc:
    type: float?
    default: 5
    label: "Include cells with the percentage of GEX transcripts mapped to mitochondrial genes not bigger than this value"
    doc: |
      Include cells with the percentage of GEX transcripts mapped to mitochondrial
      genes not bigger than this value.
    'sd:layout':
      advanced: true

  minimum_novelty_score:
    type: string?
    default: "0.8"
    label: "Include cells with the novelty score not lower than this value, calculated for GEX as log10(genes)/log10(UMIs)"
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      GEX as log10(genes)/log10(UMIs). If multiple values provided, each of them will
      be applied to the correspondent dataset from the --mex input based on the
      --identity file.
      Default: 0.8 (applied to all datasets)
    'sd:layout':
      advanced: true

  atac_minimum_cells:
    type: int?
    default: 5
    label: "Include only ATAC features detected in at least this many cells"
    doc: |
      Include only ATAC features detected in at least this many cells.
    'sd:layout':
      advanced: true

  atac_minimum_umis:
    type: string?
    default: "1000"
    label: "Include cells where at least this many ATAC UMIs (transcripts) are detected"
    doc: |
      Include cells where at least this many ATAC UMIs (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
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
      them will be applied to the correspondent dataset from the --mex input
      based on the --identity file
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
      of them will be applied to the correspondent dataset from the --mex input
      based on the --identity file.
      Default: 2 (applied to all datasets)
    'sd:layout':
      advanced: true

  minimum_frip:
    type: string?
    default: "0.15"
    label: "Include cells with the FRiP not lower than this value"
    doc: |
      Include cells with the FRiP not lower than this value. If multiple values
      provided, each of them will be applied to the correspondent dataset from
      the --mex input based on the --identity file.
      Default: 0.15 (applied to all datasets)
    'sd:layout':
      advanced: true

  maximum_blacklisted_ratio:
    type: string?
    default: "0.05"
    label: "Include cells with the ratio of reads in genomic blacklist regions not bigger than this value"
    doc: |
      Include cells with the ratio of reads in genomic blacklist regions
      not bigger than this value. If multiple values provided, each of them
      will be applied to the correspondent dataset from the --mex input based
      on the --identity file.
      Default: 0.05 (applied to all datasets)
    'sd:layout':
      advanced: true

  call_peaks:
    type:
    - type: enum
      symbols:
      - "do not call peaks"
      - "sample"
      - "RNA based cluster"
    default: "do not call peaks"
    label: "Cells grouping when calling peaks with MACS2 (discards Cell Ranger peaks)"
    doc: |
      Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.
      Peaks are called per identity (identity) or per GEX cluster (cluster) after applying
      all GEX related thresholds, maximum nucleosome signal, and minimum TSS enrichment
      score filters. If set to 'cluster' GEX clusters are identified based on the parameters
      set with --resolution, --gexndim, --highvargex, --gexnorm, and --skipgexntrg.
      Default: do not call peaks
    'sd:layout':
      advanced: true

  regress_mito_perc:
    type: boolean?
    default: false
    label: "Regress mitochondrial genes expression as a confounding source of variation. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Regress mitochondrial genes expression as a confounding source of variation
      when identifying GEX based clusters for calling custom MACS2 peaks.
      Ignored if --callpeaks is not provided.
      Default: false
    'sd:layout':
      advanced: true

  gex_normalization:
    type:
    - type: enum
      symbols:
      - "sct"
      - "log"
      - "sctglm"
    default: "sct"
    label: "Normalization method to be used when identifying GEX based clusters. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Normalization method to be used when identifying GEX based clusters for
      calling custom MACS2 peaks. Ignored if --callpeaks is not set to 'cluster'.
      Default: sct
    'sd:layout':
      advanced: true

  gex_high_var_features_count:
    type: int?
    default: 3000
    label: "Number of highly variable GEX features to detect. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Number of highly variable GEX features to detect. Used for GEX datasets
      integration, scaling, and dimensional reduction when identifying GEX based
      clusters for calling custom MACS2 peaks. Ignored if --callpeaks is not set
      to 'cluster'.
      Default: 3000
    'sd:layout':
      advanced: true

  integration_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "seurat"
      - "none"
    default: "seurat"
    label: "Integration method for GEX datasets when identifying GEX based clusters. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Integration method for GEX datasets when identifying GEX based clusters
      for calling custom MACS2 peaks. Automatically set to 'none' if --mex points
      to the Cell Ranger ARC Count outputs (single, not aggregated dataset that
      doesn't require any integration). Ignored if --callpeaks is not set to
      'cluster'.
      Default: seurat
    'sd:layout':
      advanced: true

  gex_dimensionality:
    type: int?
    default: 10
    label: "Number of principal components to use in GEX UMAP. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Dimensionality to use in GEX UMAP projection and clustering when identifying
      GEX based clusters for calling custom MACS2 peaks (from 1 to 50). If single
      number N is provided, use from 1 to N PCs. If multiple numbers are provided,
      subset to only selected PCs. Ignored if --callpeaks is not set to 'cluster'.
      Default: from 1 to 10
    'sd:layout':
      advanced: true

  resolution:
    type: float?
    default: 0.3
    label: "Resolution to be used when identifying GEX based clusters. Ignored if not calling MACS2 peaks by RNA based cluster."
    doc: |
      Resolution to be used when identifying GEX based clusters for calling
      custom MACS2 peaks. Ignored if --callpeaks is not set to 'cluster'.
      Default: 0.3
    'sd:layout':
      advanced: true

  low_memory:
    type: boolean?
    default: false
    label: "Minimize RAM usage when integrating multiple datasets (slows down processing)"
    doc: |
      Attempts to minimize RAM usage when integrating multiple datasets
      with SCTransform algorithm (slows down the computation).
      Ignored if --callpeaks is not set to 'cluster', if --ntgr is not set
      to 'seurat', if --gexnorm is not set to either 'sct' or 'glm'.
      Default: false
    'sd:layout':
      advanced: true

  parallel_memory_limit:
    type:
    - type: enum
      symbols:
      - "8"
      - "16"
      - "32"
    default: "32"
    label: "Maximum memory in GB allowed to be shared between the workers when using multiple CPUs"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32
    'sd:layout':
      advanced: true

  vector_memory_limit:
    type:
    - type: enum
      symbols:
      - "32"
      - "64"
      - "96"
    default: "64"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 64
    'sd:layout':
      advanced: true

  threads:
    type:
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
    'sd:layout':
      advanced: true


outputs:

  raw_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_peak_dnst_plot_png
    label: "Peak density per cell (not filtered)"
    doc: |
      Peak density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Peak density per cell (not filtered)'

  raw_bl_cnts_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (not filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (not filtered)'

  raw_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_pca_1_2_qc_mtrcs_plot_png
    label: "PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered)"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered)'

  raw_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered)'

  raw_cell_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_cell_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Number of cells per dataset (not filtered)'
  
  raw_gex_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gex_umi_dnst_plot_png
    label: "GEX UMI density per cell (not filtered)"
    doc: |
      GEX UMI density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'GEX UMI density per cell (not filtered)'

  raw_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_atac_umi_dnst_plot_png
    label: "ATAC UMI density per cell (not filtered)"
    doc: |
      ATAC UMI density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'ATAC UMI density per cell (not filtered)'

  raw_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_dnst_plot_png
    label: "Gene density per cell (not filtered)"
    doc: |
      Gene density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Gene density per cell (not filtered)'

  raw_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gex_atac_umi_corr_plot_png
    label: "GEX vs ATAC UMIs per cell correlation (not filtered)"
    doc: |
      GEX vs ATAC UMIs per cell correlation (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'GEX vs ATAC UMIs per cell correlation (not filtered)'

  raw_tss_enrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_tss_enrch_plot_png
    label: "TSS Enrichment Score (not filtered)"
    doc: |
      TSS Enrichment Score (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'TSS Enrichment Score (not filtered)'

  raw_frg_len_hist_png:
    type: File?
    outputSource: sc_multiome_filter/raw_frg_len_hist_png
    label: "Fragments Length Histogram (not filtered)"
    doc: |
      Fragments Length Histogram (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Fragments Length Histogram (not filtered)'

  raw_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_gene_umi_corr_plot_png
    label: "Genes vs GEX UMIs per cell correlation (not filtered)"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Genes vs GEX UMIs per cell correlation (not filtered)'

  raw_mito_perc_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_mito_perc_dnst_plot_png
    label: "Density of transcripts mapped to mitochondrial genes per cell (not filtered)"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (not filtered)'

  raw_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_nvlt_score_dnst_plot_png
    label: "Novelty score density per cell (not filtered)"
    doc: |
      Novelty score density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Novelty score density per cell (not filtered)'

  raw_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/raw_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'QC metrics densities per cell (not filtered)'


  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_peak_dnst_plot_png
    label: "Peak density per cell (intermediate filtered)"
    doc: |
      Peak density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Peak density per cell (intermediate filtered)'

  mid_fltr_bl_cnts_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (intermediate filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (intermediate filtered)'

  mid_fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_pca_1_2_qc_mtrcs_plot_png
    label: "PC1 and PC2 of ORQ-transformed QC metrics PCA (intermediate filtered)"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'PC1 and PC2 of ORQ-transformed QC metrics PCA (intermediate filtered)'

  mid_fltr_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered)'

  mid_fltr_cell_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_cell_count_plot_png
    label: "Number of cells per dataset (intermediate filtered)"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Number of cells per dataset (intermediate filtered)'

  mid_fltr_gex_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gex_umi_dnst_plot_png
    label: "GEX UMI density per cell (intermediate filtered)"
    doc: |
      GEX UMI density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'GEX UMI density per cell (intermediate filtered)'

  mid_fltr_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_atac_umi_dnst_plot_png
    label: "ATAC UMI density per cell (intermediate filtered)"
    doc: |
      ATAC UMI density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'ATAC UMI density per cell (intermediate filtered)'

  mid_fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gene_dnst_plot_png
    label: "Gene density per cell (intermediate filtered)"
    doc: |
      Gene density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Gene density per cell (intermediate filtered)'

  mid_fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gex_atac_umi_corr_plot_png
    label: "GEX vs ATAC UMIs per cell correlation (intermediate filtered)"
    doc: |
      GEX vs ATAC UMIs per cell correlation (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'GEX vs ATAC UMIs per cell correlation (intermediate filtered)'

  mid_fltr_tss_enrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_tss_enrch_plot_png
    label: "TSS Enrichment Score (intermediate filtered)"
    doc: |
      TSS Enrichment Score (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'TSS Enrichment Score (intermediate filtered)'

  mid_fltr_frg_len_hist_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_frg_len_hist_png
    label: "Fragments Length Histogram (intermediate filtered)"
    doc: |
      Fragments Length Histogram (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Fragments Length Histogram (intermediate filtered)'

  mid_fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_gene_umi_corr_plot_png
    label: "Genes vs GEX UMIs per cell correlation (intermediate filtered)"
    doc: |
      Genes vs GEX UMIs per cell correlation (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Genes vs GEX UMIs per cell correlation (intermediate filtered)'

  mid_fltr_mito_perc_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_mito_perc_dnst_plot_png
    label: "Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered)"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered)'

  mid_fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_nvlt_score_dnst_plot_png
    label: "Novelty score density per cell (intermediate filtered)"
    doc: |
      Novelty score density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Novelty score density per cell (intermediate filtered)'

  mid_fltr_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (intermediate filtered)"
    doc: |
      QC metrics densities per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'QC metrics densities per cell (intermediate filtered)'

  mid_ntgr_gex_elbow_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_ntgr_gex_elbow_plot_png
    label: "Elbow plot from GEX PCA of filtered integrated/scaled datasets"
    doc: |
      Elbow plot from GEX PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Elbow plot from GEX PCA of filtered integrated/scaled datasets'

  mid_ntgr_gex_qc_dim_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/mid_ntgr_gex_qc_dim_corr_plot_png
    label: "Correlation plots between main QC metrics and PCA reduction on GEX assay"
    doc: |
      Correlation plots between main QC metrics and PCA reduction on GEX assay.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Correlation plots between main QC metrics and PCA reduction on GEX assay'

  mid_clst_gex_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_multiome_filter/mid_clst_gex_umap_res_plot_png
    label: "Clustered UMAP projected PCA of filtered GEX datasets"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Clustered UMAP projected PCA of filtered GEX datasets'

  mid_clst_gex_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_multiome_filter/mid_clst_gex_umap_spl_by_cond_res_plot_png
    label: "Split by condition clustered UMAP projected PCA of filtered GEX datasets"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered GEX datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Split by condition clustered UMAP projected PCA of filtered GEX datasets'

  mid_fltr_gex_umap_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: sc_multiome_filter/mid_fltr_gex_umap_qc_mtrcs_res_plot_png
    label: "QC metrics for clustered UMAP projected PCA of filtered GEX datasets"
    doc: |
      QC metrics for clustered UMAP projected PCA of filtered GEX datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'QC metrics for clustered UMAP projected PCA of filtered GEX datasets'


  fltr_peak_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_peak_dnst_plot_png
    label: "Peak density per cell (filtered)"
    doc: |
      Peak density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Peak density per cell (filtered)'

  fltr_bl_cnts_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (filtered)'

  fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_pca_1_2_qc_mtrcs_plot_png
    label: "PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered)"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered)'

  fltr_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered)'

  fltr_cell_count_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_cell_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Number of cells per dataset (filtered)'

  fltr_gex_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gex_umi_dnst_plot_png
    label: "GEX UMI density per cell (filtered)"
    doc: |
      GEX UMI density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'GEX UMI density per cell (filtered)'

  fltr_atac_umi_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_atac_umi_dnst_plot_png
    label: "ATAC UMI density per cell (filtered)"
    doc: |
      ATAC UMI density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'ATAC UMI density per cell (filtered)'

  fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_dnst_plot_png
    label: "Gene density per cell (filtered)"
    doc: |
      Gene density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Gene density per cell (filtered)'

  fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gex_atac_umi_corr_plot_png
    label: "GEX vs ATAC UMIs per cell correlation (filtered)"
    doc: |
      GEX vs ATAC UMIs per cell correlation (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'GEX vs ATAC UMIs per cell correlation (filtered)'

  fltr_tss_enrch_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_tss_enrch_plot_png
    label: "TSS Enrichment Score (filtered)"
    doc: |
      TSS Enrichment Score (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'TSS Enrichment Score (filtered)'

  fltr_frg_len_hist_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_frg_len_hist_png
    label: "Fragments Length Histogram (filtered)"
    doc: |
      Fragments Length Histogram (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Fragments Length Histogram (filtered)'

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_gene_umi_corr_plot_png
    label: "Genes vs GEX UMIs per cell correlation (filtered)"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Genes vs GEX UMIs per cell correlation (filtered)'

  fltr_mito_perc_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_mito_perc_dnst_plot_png
    label: "Density of transcripts mapped to mitochondrial genes per cell (filtered)"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (filtered)'

  fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_nvlt_score_dnst_plot_png
    label: "Novelty score density per cell (filtered)"
    doc: |
      Novelty score density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Novelty score density per cell (filtered)'

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_multiome_filter/fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'QC metrics densities per cell (filtered)'


  seurat_filtered_data_rds:
    type: File?
    outputSource: sc_multiome_filter/seurat_filtered_data_rds
    label: "Filtered Seurat data in RDS format"
    doc: |
      Filtered Seurat data in RDS format
      RDS format

  sc_multiome_filter_stdout_log:
    type: File
    outputSource: sc_multiome_filter/stdout_log
    label: stdout log generated by Single-cell Multiome Filter
    doc: |
      stdout log generated by Single-cell Multiome Filter

  sc_multiome_filter_stderr_log:
    type: File
    outputSource: sc_multiome_filter/stderr_log
    label: stderr log generated by Single-cell Multiome Filter
    doc: |
      stderr log generated by Single-cell Multiome Filter


steps:

  uncompress_feature_bc_matrices:
    in:
      compressed: filtered_feature_bc_matrix_folder
    out:
    - uncompressed
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        compressed:
          type: File
          inputBinding:
            position: 1
      outputs:
        uncompressed:
          type: Directory
          outputBinding:
            glob: "*"
      baseCommand: ["tar", "xzf"]

  sc_multiome_filter:
    run: ../tools/sc-multiome-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      aggregation_metadata: aggregation_metadata
      atac_fragments_file: atac_fragments_file
      annotation_gtf_file: annotation_gtf_file
      grouping_data: grouping_data
      blacklisted_regions_file: blacklisted_regions_file
      barcodes_data: barcodes_data
      gex_minimum_cells: gex_minimum_cells
      gex_minimum_features:
        source: gex_minimum_features
        valueFrom: $(split_numbers(self))
      gex_maximum_features:
        source: gex_maximum_features
        valueFrom: $(split_numbers(self))
      gex_minimum_umis:
        source: gex_minimum_umis
        valueFrom: $(split_numbers(self))
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      regress_mito_perc: regress_mito_perc
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      atac_minimum_cells: atac_minimum_cells
      atac_minimum_umis:
        source: atac_minimum_umis
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
      maximum_blacklisted_ratio:
        source: maximum_blacklisted_ratio
        valueFrom: $(split_numbers(self))
      call_peaks:
        source: call_peaks
        valueFrom: |
          ${
            if (self == "do not call peaks") {
              return null;
            } else if (self == "sample") {
              return "identity";
            } else if (self == "RNA based cluster") {
              return "cluster";
            }
          }
      gex_normalization: gex_normalization
      gex_high_var_features_count: gex_high_var_features_count
      integration_method: integration_method
      gex_dimensionality: gex_dimensionality
      resolution: resolution
      low_memory: low_memory
      parallel_memory_limit:
        source: parallel_memory_limit
        valueFrom: $(parseInt(self))
      threads:
        source: threads
        valueFrom: $(parseInt(self))
      export_pdf_plots:
        default: false
      export_h5seurat_data:
        default: false
      vector_memory_limit:
        source: vector_memory_limit
        valueFrom: $(parseInt(self))
      verbose:
        default: true
    out:
    - raw_peak_dnst_plot_png
    - raw_bl_cnts_dnst_plot_png
    - raw_pca_1_2_qc_mtrcs_plot_png
    - raw_pca_2_3_qc_mtrcs_plot_png
    - raw_cell_count_plot_png
    - raw_gex_umi_dnst_plot_png
    - raw_atac_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_gex_atac_umi_corr_plot_png
    - raw_tss_enrch_plot_png
    - raw_frg_len_hist_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_perc_dnst_plot_png
    - raw_nvlt_score_dnst_plot_png
    - raw_qc_mtrcs_plot_png
    - mid_fltr_peak_dnst_plot_png
    - mid_fltr_bl_cnts_dnst_plot_png
    - mid_fltr_pca_1_2_qc_mtrcs_plot_png
    - mid_fltr_pca_2_3_qc_mtrcs_plot_png
    - mid_fltr_cell_count_plot_png
    - mid_fltr_gex_umi_dnst_plot_png
    - mid_fltr_atac_umi_dnst_plot_png
    - mid_fltr_gene_dnst_plot_png
    - mid_fltr_gex_atac_umi_corr_plot_png
    - mid_fltr_tss_enrch_plot_png
    - mid_fltr_frg_len_hist_png
    - mid_fltr_gene_umi_corr_plot_png
    - mid_fltr_mito_perc_dnst_plot_png
    - mid_fltr_nvlt_score_dnst_plot_png
    - mid_fltr_qc_mtrcs_plot_png
    - mid_ntgr_gex_elbow_plot_png 
    - mid_ntgr_gex_qc_dim_corr_plot_png
    - mid_clst_gex_umap_res_plot_png
    - mid_clst_gex_umap_spl_by_cond_res_plot_png
    - mid_fltr_gex_umap_qc_mtrcs_res_plot_png
    - fltr_peak_dnst_plot_png
    - fltr_bl_cnts_dnst_plot_png
    - fltr_pca_1_2_qc_mtrcs_plot_png
    - fltr_pca_2_3_qc_mtrcs_plot_png
    - fltr_cell_count_plot_png
    - fltr_gex_umi_dnst_plot_png
    - fltr_atac_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_gex_atac_umi_corr_plot_png
    - fltr_tss_enrch_plot_png
    - fltr_frg_len_hist_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_perc_dnst_plot_png
    - fltr_nvlt_score_dnst_plot_png
    - fltr_qc_mtrcs_plot_png
    - seurat_filtered_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Multiome Filter"
s:name: "Single-cell Multiome Filter"
s:alternateName: "Filters single-cell multiome datasets based on the common QC metrics"

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
  Single-cell Multiome Filter
  =====================================================================
  Filters single-cell multiome datasets based on the common QC metrics.