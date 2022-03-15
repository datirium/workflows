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
    - var get_correct_order = function(cluster_files=null, cluster_aliases=null, aggregation_metadata=null) {
        if (cluster_files === null || cluster_aliases === null){
          return [null, null];
        } else if (aggregation_metadata === null) {
          if (cluster_files.length == 1 && cluster_aliases.length == 1){
            return [cluster_files, cluster_aliases];
          } else {
            return [null, null];
          }
        } else {
          let updated_cluster_aliases = cluster_aliases.map(alias => alias.replace(/\t|\s|\[|\]|\>|\<|,|\./g, "_"));
          let ordered_cluster_aliases = aggregation_metadata.contents.split("\n").map(line => line.split(",")[0]).filter(value => value != "" && value != "library_id");
          let ordered_cluster_files = ordered_cluster_aliases.map(alias => cluster_files[updated_cluster_aliases.indexOf(alias)]);
          if (ordered_cluster_files.includes(undefined)){
            return [null, null];
          } else {
            return [ordered_cluster_files, ordered_cluster_aliases];
          }
        }
      };


'sd:upstream':
  sc_arc_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-count.cwl"
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-aggr.cwl"
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  genotype_sample:
  - "souporcell.cwl"


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

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format.
    'sd:upstreamSource': "sc_arc_sample/atac_fragments_file"
    'sd:localLabel': true

  annotation_gtf_file:
    type: File
    label: "Cell Ranger ARC Count/Aggregate Experiment"
    doc: |
      GTF annotation file that includes refGene and mitochondrial DNA annotations.
    'sd:upstreamSource': "sc_arc_sample/genome_indices/genome_indices/annotation_gtf"
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

  gex_genotype_cluster_tsv_file:
    type:
      - "null"
      - File[]
    label: "Souporcell Cluster by Genotype Experiment(s)"
    doc: |
      Cellurar barcodes file(s) clustered by genotype (GEX) generated in Souporcell
    'sd:upstreamSource': "genotype_sample/gex_genotype_cluster_tsv_file"
    'sd:localLabel': true

  atac_genotype_cluster_tsv_file:
    type:
      - "null"
      - File[]
    label: "Souporcell Cluster by Genotype Experiment(s)"
    doc: |
      Cellurar barcodes file(s) clustered by genotype (ATAC) generated in Souporcell
    'sd:upstreamSource': "genotype_sample/atac_genotype_cluster_tsv_file"
    'sd:localLabel': true

  genotype_cluster_names:
    type:
      - "null"
      - string[]
    label: "Souporcell Cluster by Genotype Experiment(s)"
    doc: |
      Aliases of the scRNA-Seq Cell Ranger ARC Experiments used in Souporcell analysis.
      The order of samples will be defined based on the rows from aggregation_metadata
      input if it was provided.
    'sd:upstreamSource': "genotype_sample/sc_rnaseq_sample/alias"
    'sd:localLabel': true

  conditions_data:
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
      Path to the blacklisted regions file in BED format.

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

  regress_mito_perc:
    type: boolean?
    default: false
    label: "Regress mitochondrial genes expression as a confounding source of variation"
    doc: |
      Regress mitochondrial genes expression as a confounding source of variation.
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

  gex_selected_features:
    type: string?
    default: null
    label: "GEX features of interest to evaluate expression"
    doc: |
      GEX features of interest to evaluate expression.
      Default: do not highlight any features
    'sd:layout':
      advanced: true

  gex_high_var_features_count:
    type: int?
    default: 3000
    label: "Number of highly variable GEX features to detect. Used for datasets integration, scaling, and dimensional reduction"
    doc: |
      Number of highly variable GEX features to detect. Used for GEX datasets
      integration, scaling, and dimensional reduction.
    'sd:layout':
      advanced: true

  gex_dimensionality:
    type: int?
    default: 10
    label: "Number of principal components to use in GEX UMAP projection and clustering (from 1 to 50)"
    doc: |
      Number of principal components to use in GEX UMAP projection and clustering
      (from 1 to 50).
    'sd:layout':
      advanced: true

  gex_minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those GEX features that on average have log fold change difference in expression between every tested pair of clusters not lower than this value"
    doc: |
      For putative gene markers identification include only those GEX features that
      on average have log fold change difference in expression between every tested
      pair of clusters not lower than this value.
    'sd:layout':
      advanced: true

  gex_minimum_pct:
    type: float?
    default: 0.1
    label: "Include only those GEX features that are detected in not lower than this fraction of cells in either of the two tested clusters"
    doc: |
      For putative gene markers identification include only those GEX features that
      are detected in not lower than this fraction of cells in either of the two
      tested clusters.
    'sd:layout':
      advanced: true

  gex_only_positive_markers:
    type: boolean?
    default: false
    label: "For putative gene markers identification return only positive markers"
    doc: |
      For putative gene markers identification return only positive markers.
    'sd:layout':
      advanced: true

  gex_test_use:
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
    default: "wilcox"
    label: "Statistical test to use for putative gene markers identification"
    doc: |
      Statistical test to use for putative gene markers identification.
    'sd:layout':
      advanced: true

  no_sct:
    type: boolean?
    default: false
    label: "Do not use SCTransform when running RNA datasets integration"
    doc: |
      Do not use SCTransform when running RNA datasets integration.
      Use LogNormalize instead. Ignored when --mex points to the
      Cell Ranger ARC Count outputs (single, not aggregated dataset
      that doesn't require any integration) or --skipgexntrg parameter
      was applied.
    'sd:layout':
      advanced: true

  # skip_miqc:
  #   type: boolean?
  #   default: true
  #   label: "Skip miQC threshold prediction for the percentage of transcripts mapped to mitochondrial genes"
  #   doc: |
  #     Skip threshold prediction for the percentage of transcripts
  #     mapped to mitochondrial genes (do not run MiQC).
  #     Default: false
  #   'sd:layout':
  #     advanced: true

  skip_gex_ntrg:
    type: boolean?
    default: false
    label: "Do not integrate RNA datasets, use merged data instead"
    doc: |
      Do not integrate RNA datasets, use merged data instead.
      Applied by default if --mex points to the Cell Ranger ARC Count
      outputs (single, not aggregated dataset that doesn't require any
      integration).
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
      - "identity"
      - "cluster"
    default: "do not call peaks"
    label: "Cells grouping when calling peaks with MACS2 (use instead of Cell Ranger peaks)"
    doc: |
      Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.
      Peaks are called per GEX cluster (cluster) or per identity (identity) after
      applying all GEX related thresholds, maximum nucleosome signal, and minimum TSS
      enrichment score filters. If set to 'cluster' GEX clusters are identified based on the
      first value from the --resolution array using --gexndim principal components when
      building nearest-neighbour graph.
      Default: do not call peaks
    'sd:layout':
      advanced: true

  skip_atac_ntrg:
    type: boolean?
    default: false
    label: "Do not integrate ATAC datasets, use merged data instead"
    doc: |
      Do not integrate ATAC datasets, use merged data instead.
      Applied by default if --mex pointed to the Cell Ranger ARC Count
      outputs (single, not aggregated dataset that doesn't require any
      integration).
    'sd:layout':
      advanced: true

  atac_dimensionality:
    type: int?
    default: 10
    label: "Number of principal components to use in ATAC UMAP projection and clustering (from 2 to 50)"
    doc: |
      Number of principal components to use in ATAC UMAP projection and clustering
      (from 2 to 50).
    'sd:layout':
      advanced: true

  atac_high_var_features_perc:
    type: int?
    default: 75
    label: "Minimum percentile to set the top most common ATAC features as highly variable"
    doc: |
      Minimum percentile to set the top most common ATAC features as highly variable.
      For example, setting to 5 will use the the top 95% most common among all cells
      ATAC features as highly variable. Used for ATAC datasets integration, scaling,
      and dimensional reduction.
      Default: 75 (use only the top 25% of all common peaks)
    'sd:layout':
      advanced: true

  resolution:
    type: string?
    default: "0.3"
    label: "Comma or space separated list of clustering resolutions"
    doc: |
      Comma or space separated list of clustering resolutions
    'sd:layout':
      advanced: true

  memory:
    type: int?
    default: 32
    label: "Maximum memory in GB allowed to be shared between the workers when using multiple CPUs"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
    'sd:layout':
      advanced: true


outputs:

  raw_cell_count_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_cell_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Number of cells per dataset (not filtered)'
  
  raw_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_pca_1_2_qc_mtrcs_plot_png
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
    outputSource: seurat_wnn_cluster/raw_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered)'

  raw_gex_umi_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_gex_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/raw_atac_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/raw_gene_dnst_plot_png
    label: "Gene density per cell (not filtered)"
    doc: |
      Gene density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Gene density per cell (not filtered)'

  raw_peak_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_peak_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/raw_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (not filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (not filtered)'

  raw_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_gex_atac_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/raw_tss_enrch_plot_png
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
    outputSource: seurat_wnn_cluster/raw_frg_len_hist_png
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
    outputSource: seurat_wnn_cluster/raw_gene_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/raw_mito_perc_dnst_plot_png
    label: "Density of transcripts mapped to mitochondrial genes per cell (not filtered)"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (not filtered)'

  raw_miqc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_miqc_mtrcs_plot_png
    label: "MiQC prediction of the compromised cells level (not filtered)"
    doc: |
      MiQC prediction of the compromised cells level (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'MiQC prediction of the compromised cells level (not filtered)'

  raw_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_nvlt_score_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/raw_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'QC metrics densities per cell (not filtered)'

  # raw_qc_mtrcs_tsv:
  #   type: File?
  #   outputSource: seurat_wnn_cluster/raw_qc_mtrcs_tsv
  #   label: "QC metrics densities per cell (not filtered)"
  #   doc: |
  #     QC metrics densities per cell (not filtered).
  #     TSV format
  #   'sd:visualPlugins':
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'GEX UMI density per cell (not filtered)'
  #       xAxisTitle: 'GEX UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#2e8b57"]
  #       height: 300
  #       data: [$1, $2]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Gene density per cell (not filtered)'
  #       xAxisTitle: 'Genes per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#191970"]
  #       height: 300
  #       data: [$3, $4]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of transcripts mapped to mitochondrial genes per cell (not filtered)'
  #       xAxisTitle: 'Percentage of transcripts mapped to mitochondrial genes per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ff0000"]
  #       height: 300
  #       data: [$5, $6]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Novelty score density per cell (not filtered)'
  #       xAxisTitle: 'log10 Genes / log10 UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ffd700"]
  #       height: 300
  #       data: [$7, $8]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'ATAC UMI density per cell (not filtered)'
  #       xAxisTitle: 'ATAC UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#7cfc00"]
  #       height: 300
  #       data: [$9, $10]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Peak density per cell (not filtered)'
  #       xAxisTitle: 'Peaks per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#e9967a"]
  #       height: 300
  #       data: [$11, $12]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of mononucleosomal to nucleosome-free fragments ratio per cell (not filtered)'
  #       xAxisTitle: 'Mononucleosomal to nucleosome-free fragments ratio per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#00bfff"]
  #       height: 300
  #       data: [$13, $14]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of FRiP per cell (not filtered)'
  #       xAxisTitle: 'Fraction of reads in peaks per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#0000ff"]
  #       height: 300
  #       data: [$15, $16]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of fraction of reads within blacklisted regions per cell (not filtered)'
  #       xAxisTitle: 'Fraction of reads within blacklisted regions per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ff1493"]
  #       height: 300
  #       data: [$17, $18]

  mid_fltr_cell_count_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_cell_count_plot_png
    label: "Number of cells per dataset (intermediate filtered)"
    doc: |
      Number of cells per dataset (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Number of cells per dataset (intermediate filtered)'

  mid_fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_pca_1_2_qc_mtrcs_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (intermediate filtered)'

  mid_fltr_gex_umi_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_gex_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_atac_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_gene_dnst_plot_png
    label: "Gene density per cell (intermediate filtered)"
    doc: |
      Gene density per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Gene density per cell (intermediate filtered)'

  mid_fltr_peak_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_peak_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (intermediate filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (intermediate filtered)'

  mid_fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_gex_atac_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_tss_enrch_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_frg_len_hist_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_gene_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_mito_perc_dnst_plot_png
    label: "Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered)"
    doc: |
      Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (intermediate filtered)'

  mid_fltr_miqc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_miqc_mtrcs_plot_png
    label: "MiQC prediction of the compromised cells level (intermediate filtered)"
    doc: |
      MiQC prediction of the compromised cells level (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'MiQC prediction of the compromised cells level (intermediate filtered)'

  mid_fltr_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/mid_fltr_nvlt_score_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/mid_fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (intermediate filtered)"
    doc: |
      QC metrics densities per cell (intermediate filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1.5. Intermediate filtered QC'
        Caption: 'QC metrics densities per cell (intermediate filtered)'

  fltr_cell_count_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_cell_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Number of cells per dataset (filtered)'

  fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_pca_1_2_qc_mtrcs_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_pca_2_3_qc_mtrcs_plot_png
    label: "PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered)"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered)'

  fltr_gex_umi_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_gex_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_atac_umi_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_gene_dnst_plot_png
    label: "Gene density per cell (filtered)"
    doc: |
      Gene density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Gene density per cell (filtered)'

  fltr_peak_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_peak_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_bl_cnts_dnst_plot_png
    label: "Density of fraction of reads within blacklisted regions per cell (filtered)"
    doc: |
      Density of fraction of reads within blacklisted regions per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Density of fraction of reads within blacklisted regions per cell (filtered)'

  fltr_gex_atac_umi_corr_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_gex_atac_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_tss_enrch_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_frg_len_hist_png
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
    outputSource: seurat_wnn_cluster/fltr_gene_umi_corr_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_mito_perc_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_nvlt_score_dnst_plot_png
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
    outputSource: seurat_wnn_cluster/fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'QC metrics densities per cell (filtered)'

  # fltr_qc_mtrcs_tsv:
  #   type: File?
  #   outputSource: seurat_wnn_cluster/fltr_qc_mtrcs_tsv
  #   label: "QC metrics densities per cell (filtered)"
  #   doc: |
  #     QC metrics densities per cell (filtered).
  #     TSV format
  #   'sd:visualPlugins':
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'GEX UMI density per cell (filtered)'
  #       xAxisTitle: 'GEX UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#2e8b57"]
  #       height: 300
  #       data: [$1, $2]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Gene density per cell (filtered)'
  #       xAxisTitle: 'Genes per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#191970"]
  #       height: 300
  #       data: [$3, $4]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of transcripts mapped to mitochondrial genes per cell (filtered)'
  #       xAxisTitle: 'Percentage of transcripts mapped to mitochondrial genes per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ff0000"]
  #       height: 300
  #       data: [$5, $6]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Novelty score density per cell (filtered)'
  #       xAxisTitle: 'log10 Genes / log10 UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ffd700"]
  #       height: 300
  #       data: [$7, $8]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'ATAC UMI density per cell (filtered)'
  #       xAxisTitle: 'ATAC UMIs per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#7cfc00"]
  #       height: 300
  #       data: [$9, $10]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Peak density per cell (filtered)'
  #       xAxisTitle: 'Peaks per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#e9967a"]
  #       height: 300
  #       data: [$11, $12]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of mononucleosomal to nucleosome-free fragments ratio per cell (filtered)'
  #       xAxisTitle: 'Mononucleosomal to nucleosome-free fragments ratio per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#00bfff"]
  #       height: 300
  #       data: [$13, $14]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of FRiP per cell (filtered)'
  #       xAxisTitle: 'Fraction of reads in peaks per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#0000ff"]
  #       height: 300
  #       data: [$15, $16]
  #   - scatter:
  #       tab: 'QC (interactive)'
  #       Title: 'Density of fraction of reads within blacklisted regions per cell (filtered)'
  #       xAxisTitle: 'Fraction of reads within blacklisted regions per cell'
  #       yAxisTitle: 'Density'
  #       colors: ["#ff1493"]
  #       height: 300
  #       data: [$17, $18]

  ntgr_gex_elbow_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_gex_elbow_plot_png
    label: "Elbow plot from GEX PCA of filtered integrated/scaled datasets"
    doc: |
      Elbow plot from GEX PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'Elbow plot from GEX PCA of filtered integrated/scaled datasets'

  ntgr_gex_qc_dim_corr_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_gex_qc_dim_corr_plot_png
    label: "Correlation plots between main QC metrics and PCA reduction on GEX assay"
    doc: |
      Correlation plots between main QC metrics and PCA reduction on GEX assay.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'Correlation plots between main QC metrics and PCA reduction on GEX assay'

  ntgr_atac_qc_dim_corr_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_atac_qc_dim_corr_plot_png
    label: "Correlation plots between main QC metrics and LSI reduction on ATAC assay"
    doc: |
      Correlation plots between main QC metrics and LSI reduction on ATAC assay.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'Correlation plots between main QC metrics and LSI reduction on ATAC assay'

  ntgr_gex_pca_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_gex_pca_plot_png
    label: "GEX PCA of filtered integrated/scaled datasets"
    doc: |
      GEX PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 3. Dimensionality evaluation'
        Caption: 'GEX PCA of filtered integrated/scaled datasets'

  clst_gex_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_gex_umap_res_plot_png
    label: "Clustered UMAP projected PCA of filtered GEX datasets"
    doc: |
      Clustered UMAP projected PCA of filtered GEX datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Clustered UMAP projected PCA of filtered GEX datasets'

  clst_gex_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_gex_umap_spl_by_cond_res_plot_png
    label: "Split by condition clustered UMAP projected PCA of filtered GEX datasets"
    doc: |
      Split by condition clustered UMAP projected PCA of filtered GEX datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Split by condition clustered UMAP projected PCA of filtered GEX datasets'

  clst_atac_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_atac_umap_res_plot_png
    label: "Clustered UMAP projected LSI of filtered ATAC datasets"
    doc: |
      Clustered UMAP projected LSI of filtered ATAC datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Clustered UMAP projected LSI of filtered ATAC datasets'

  clst_atac_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_atac_umap_spl_by_cond_res_plot_png
    label: "Split by condition clustered UMAP projected LSI of filtered ATAC datasets"
    doc: |
      Split by condition clustered UMAP projected LSI of filtered ATAC datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Split by condition clustered UMAP projected LSI of filtered ATAC datasets'

  clst_wnn_umap_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_wnn_umap_res_plot_png
    label: "Clustered UMAP projected WNN"
    doc: |
      Clustered UMAP projected WNN.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Clustered UMAP projected WNN'

  clst_wnn_umap_spl_by_cond_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_wnn_umap_spl_by_cond_res_plot_png
    label: "Split by condition clustered UMAP projected WNN"
    doc: |
      Split by condition clustered UMAP projected WNN.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'Split by condition clustered UMAP projected WNN'

  clst_wnn_qc_mtrcs_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/clst_wnn_qc_mtrcs_res_plot_png
    label: "QC metrics for clustered UMAP projected WNN"
    doc: |
      QC metrics for clustered UMAP projected WNN.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 4. Clustering'
        Caption: 'QC metrics for clustered UMAP projected WNN'

  expr_avg_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/expr_avg_per_clst_res_plot_png
    label: "Scaled average log normalized gene expression per cluster"
    doc: |
      Scaled average log normalized gene expression per cluster.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Gene expression'
        Caption: 'Scaled average log normalized gene expression per cluster'
  
  expr_per_clst_cell_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/expr_per_clst_cell_res_plot_png
    label: "Log normalized gene expression per cell of clustered datasets"
    doc: |
      Log normalized gene expression per cell of clustered datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Gene expression'
        Caption: 'Log normalized gene expression per cell of clustered datasets'

  expr_dnst_per_clst_res_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: seurat_wnn_cluster/expr_dnst_per_clst_res_plot_png
    label: "Log normalized gene expression densities per cluster"
    doc: |
      Log normalized gene expression densities per cluster.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 5. Gene expression'
        Caption: 'Log normalized gene expression densities per cluster'

  compressed_cellbrowser_config_data:
    type: File
    outputSource: compress_cellbrowser_config_data/compressed_folder
    label: "Compressed directory with UCSC Cellbrowser configuration data"
    doc: |
      Compressed directory with UCSC Cellbrowser configuration data
  
  cellbrowser_html_data:
    type: Directory
    outputSource: seurat_wnn_cluster/cellbrowser_html_data
    label: "Directory with UCSC Cellbrowser formatted html data"
    doc: |
      Directory with UCSC Cellbrowser formatted html data

  cellbrowser_html_file:
    type: File
    outputSource: seurat_wnn_cluster/cellbrowser_html_file
    label: "Open in UCSC Cell Browser"
    doc: |
      HTML index file from the directory with UCSC Cellbrowser formatted html data
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  gex_clst_pttv_gene_markers:
    type: File
    outputSource: seurat_wnn_cluster/gex_clst_pttv_gene_markers
    label: "GEX putative gene markers file for all clusters and all resolutions"
    doc: |
      GEX putative gene markers file for all clusters and all resolutions.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Step 6. Putative gene markers'
        Title: 'Putative gene markers'

  seurat_clst_data_rds:
    type: File?
    outputSource: seurat_wnn_cluster/seurat_clst_data_rds
    label: "Clustered filtered integrated/scaled Seurat data"
    doc: |
      Clustered filtered integrated/scaled Seurat data.
      RDS format

  seurat_wnn_stdout_log:
    type: File
    outputSource: seurat_wnn_cluster/stdout_log
    label: stdout log generated by Seurat WNN Analysis
    doc: |
      stdout log generated by Seurat WNN Analysis

  seurat_wnn_stderr_log:
    type: File
    outputSource: seurat_wnn_cluster/stderr_log
    label: stderr log generated by Seurat WNN Analysis
    doc: |
      stderr log generated by Seurat WNN Analysis


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

  prepare_metadata:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string?
          default: |
            set -- "$0" "$@"
            if [ ${#@} -lt 3 ]; then
                echo "Not enough arguments provided"
                exit 0
            else
              echo "Provided arguments:"
              for i in "$@";
                  do echo $i;
              done;
            fi

            FILES=()
            ALIASES=()

            for ITEM in "$@"; do
                if [ -f "$ITEM" ]; then
                    echo "Add ${ITEM} as a file"
                    FILES+=("$ITEM")
                else
                    echo "Add ${ITEM} as an alias"
                    ALIASES+=("$ITEM")
                fi
            done;

            CENTER="${#FILES[@]}/2"
            GEX_FILES=( "${FILES[@]:0:$CENTER}" )
            ATAC_FILES=( "${FILES[@]:$CENTER}" )

            echo "Processing GEX clusters"
            echo -e "barcode\tgex_genotype" > gex_metadata.tsv
            let CNT=0
            for FILE in ${GEX_FILES[@]}; do
                ALIAS=${ALIASES[CNT]}
                echo "Saving rows from ${FILE} with label ${ALIAS} and barcode $((CNT+1))"
                cat "$FILE" | grep -v "barcode" | cut -f 1-3 | awk -v i="$((CNT+1))" -v c="${ALIAS}" '{split($1,b,"-"); if ($2 == "singlet") print b[1]"-"i"\t"c"__singlet_"$3; else print b[1]"-"i"\t"c"__"$2}' >> gex_metadata.tsv
                (( CNT++ ))
            done;

            echo "Processing ATAC clusters"
            echo "atac_genotype" > atac_metadata.tsv
            let CNT=0
            for FILE in ${ATAC_FILES[@]}; do
                ALIAS=${ALIASES[CNT]}
                echo "Saving rows from ${FILE} with label ${ALIAS} and barcode $((CNT+1))"
                cat "$FILE" | grep -v "barcode" | cut -f 1-3 | awk -v c="${ALIAS}" '{if ($2 == "singlet") print c"__singlet_"$3; else print c"__"$2}' >> atac_metadata.tsv
                (( CNT++ ))
            done;

            echo "Combining results"
            paste gex_metadata.tsv atac_metadata.tsv > extra_metadata.tsv
            rm -f gex_metadata.tsv atac_metadata.tsv
            head extra_metadata.tsv
          inputBinding:
            position: 5
        gex_genotype_cluster_tsv_file:
          type:
            - "null"
            - File[]
          inputBinding:
            position: 6
        atac_genotype_cluster_tsv_file:
          type:
            - "null"
            - File[]
          inputBinding:
            position: 7
        genotype_cluster_names:
          type:
            - "null"
            - string[]
          inputBinding:
            position: 8
      outputs:
        metadata_file:
          type: File?
          outputBinding:
            glob: "extra_metadata.tsv"
      baseCommand: ["bash", "-c"]
    in:
      gex_genotype_cluster_tsv_file:
        source: [gex_genotype_cluster_tsv_file, genotype_cluster_names, aggregation_metadata]
        valueFrom: $(get_correct_order(self[0], self[1], self[2])[0])
      atac_genotype_cluster_tsv_file:
        source: [atac_genotype_cluster_tsv_file, genotype_cluster_names, aggregation_metadata]
        valueFrom: $(get_correct_order(self[0], self[1], self[2])[0])
      genotype_cluster_names:
        source: [gex_genotype_cluster_tsv_file, genotype_cluster_names, aggregation_metadata]      # doesn't matter either use GEX or ATAC, as the size and order should be the same
        valueFrom: $(get_correct_order(self[0], self[1], self[2])[1])
    out:
    - metadata_file

  seurat_wnn_cluster:
    run: ../tools/seurat-wnn-cluster.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      atac_fragments_file: atac_fragments_file
      annotation_gtf_file: annotation_gtf_file
      aggregation_metadata: aggregation_metadata
      conditions_data: conditions_data
      blacklisted_regions_file: blacklisted_regions_file
      barcodes_data: barcodes_data
      metadata_file: prepare_metadata/metadata_file
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
        valueFrom: $(self=="do not call peaks"?null:self)
      gex_selected_features:
        source: gex_selected_features
        valueFrom: $(split_features(self))
      gex_high_var_features_count: gex_high_var_features_count
      no_sct: no_sct
      skip_gex_ntrg: skip_gex_ntrg
      skip_atac_ntrg: skip_atac_ntrg
      # skip_miqc: skip_miqc
      skip_miqc:
        default: true                                            # running miQC dramatically increases memory footprint and running time
      gex_dimensionality: gex_dimensionality
      gex_minimum_logfc: gex_minimum_logfc
      gex_minimum_pct: gex_minimum_pct
      gex_only_positive_markers: gex_only_positive_markers
      gex_test_use: gex_test_use
      atac_dimensionality: atac_dimensionality
      atac_high_var_features_perc: atac_high_var_features_perc
      resolution:
        source: resolution
        valueFrom: $(split_numbers(self))
      verbose:
        default: true
      export_pdf_plots:
        default: false
      export_rds_data:
        default: true
      memory: memory
      threads: threads
    out:
    - raw_cell_count_plot_png
    - raw_pca_1_2_qc_mtrcs_plot_png
    - raw_pca_2_3_qc_mtrcs_plot_png
    - raw_gex_umi_dnst_plot_png
    - raw_atac_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_peak_dnst_plot_png
    - raw_bl_cnts_dnst_plot_png
    - raw_gex_atac_umi_corr_plot_png
    - raw_tss_enrch_plot_png
    - raw_frg_len_hist_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_perc_dnst_plot_png
    - raw_miqc_mtrcs_plot_png
    - raw_nvlt_score_dnst_plot_png
    - raw_qc_mtrcs_plot_png
    # - raw_qc_mtrcs_tsv
    - mid_fltr_cell_count_plot_png
    - mid_fltr_pca_1_2_qc_mtrcs_plot_png
    - mid_fltr_pca_2_3_qc_mtrcs_plot_png
    - mid_fltr_gex_umi_dnst_plot_png
    - mid_fltr_atac_umi_dnst_plot_png
    - mid_fltr_gene_dnst_plot_png
    - mid_fltr_peak_dnst_plot_png
    - mid_fltr_bl_cnts_dnst_plot_png
    - mid_fltr_gex_atac_umi_corr_plot_png
    - mid_fltr_tss_enrch_plot_png
    - mid_fltr_frg_len_hist_png
    - mid_fltr_gene_umi_corr_plot_png
    - mid_fltr_mito_perc_dnst_plot_png
    - mid_fltr_miqc_mtrcs_plot_png
    - mid_fltr_nvlt_score_dnst_plot_png
    - mid_fltr_qc_mtrcs_plot_png
    - fltr_cell_count_plot_png
    - fltr_pca_1_2_qc_mtrcs_plot_png
    - fltr_pca_2_3_qc_mtrcs_plot_png
    - fltr_gex_umi_dnst_plot_png
    - fltr_atac_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_peak_dnst_plot_png
    - fltr_bl_cnts_dnst_plot_png
    - fltr_gex_atac_umi_corr_plot_png
    - fltr_tss_enrch_plot_png
    - fltr_frg_len_hist_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_perc_dnst_plot_png
    - fltr_nvlt_score_dnst_plot_png
    - fltr_qc_mtrcs_plot_png
    # - fltr_qc_mtrcs_tsv
    - ntgr_gex_elbow_plot_png
    - ntgr_gex_qc_dim_corr_plot_png
    - ntgr_atac_qc_dim_corr_plot_png
    - ntgr_gex_pca_plot_png
    - clst_gex_umap_res_plot_png
    - clst_gex_umap_spl_by_cond_res_plot_png
    - clst_atac_umap_res_plot_png
    - clst_atac_umap_spl_by_cond_res_plot_png
    - clst_wnn_umap_res_plot_png
    - clst_wnn_umap_spl_by_cond_res_plot_png
    - clst_wnn_qc_mtrcs_res_plot_png
    - gex_clst_pttv_gene_markers
    - seurat_clst_data_rds
    - expr_avg_per_clst_res_plot_png
    - expr_per_clst_cell_res_plot_png
    - expr_dnst_per_clst_res_plot_png
    - cellbrowser_config_data
    - cellbrowser_html_data
    - cellbrowser_html_file
    - stdout_log
    - stderr_log

  compress_cellbrowser_config_data:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: seurat_wnn_cluster/cellbrowser_config_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Seurat WNN Analysis"
s:name: "Seurat WNN Analysis"
s:alternateName: "Runs Seurat Weighted Nearest Neighbor Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/seurat-wnn-cluster.cwl
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
  Seurat WNN Analysis
  ===================
  Runs Seurat Weighted Nearest Neighbor Analysis.