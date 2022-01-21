cwlVersion: v1.0
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
  sc_arc_count_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-count.cwl"
  - "cellranger-arc-count.cwl"
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
    label: "Cell Ranger ARC Count Experiment"
    doc: |
      Filtered feature barcode matrix stored as a CSC sparse matrix in MEX format.
      The rows consist of all the gene and peak features concatenated together
      (identical to raw feature barcode matrix) and the columns are restricted to
      those barcodes that are identified as cells.
    'sd:upstreamSource': "sc_arc_count_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ARC Count Experiment"
    doc: |
      Count and barcode information for every ATAC fragment observed in
      the experiment in TSV format.
    'sd:upstreamSource': "sc_arc_count_sample/atac_fragments_file"
    'sd:localLabel': true

  annotation_gtf_file:
    type: File
    label: "Cell Ranger ARC Count Experiment"
    doc: |
      GTF annotation file that includes refGene and mitochondrial DNA annotations.
    'sd:upstreamSource': "sc_arc_count_sample/genome_indices/genome_indices/annotation_gtf"
    'sd:localLabel': true

  gex_genotype_cluster_tsv_file:
    type: File?
    label: "Souporcell Cluster by Genotype Experiment"
    doc: |
      Cellurar barcodes file clustered by genotype (GEX) generated in Souporcell
    'sd:upstreamSource': "genotype_sample/gex_genotype_cluster_tsv_file"
    'sd:localLabel': true

  atac_genotype_cluster_tsv_file:
    type: File?
    label: "Souporcell Cluster by Genotype Experiment"
    doc: |
      Cellurar barcodes file clustered by genotype (ATAC) generated in Souporcell
    'sd:upstreamSource': "genotype_sample/atac_genotype_cluster_tsv_file"
    'sd:localLabel': true

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
    type: int?
    default: 250
    label: "Include cells where at least this many GEX features are detected"
    doc: |
      Include cells where at least this many GEX features are detected.
    'sd:layout':
      advanced: true

  gex_maximum_features:
    type: int?
    default: 5000
    label: "Include cells with the number of GEX features not bigger than this value"
    doc: |
      Include cells with the number of GEX features not bigger than this value.
    'sd:layout':
      advanced: true

  gex_minimum_umis:
    type: int?
    default: 500
    label: "Include cells where at least this many GEX UMIs (transcripts) are detected"
    doc: |
      Include cells where at least this many GEX UMIs (transcripts) are detected.
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
    type: float?
    default: 0.8
    label: "Include cells with the novelty score not lower than this value, calculated for GEX as log10(genes)/log10(UMIs)"
    doc: |
      Include cells with the novelty score not lower than this value
      calculated for GEX as log10(genes)/log10(UMIs).
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
    label: "Number of highly variable features to detect. Used for datasets integration, scaling, and dimensional reduction"
    doc: |
      Number of highly variable features to detect. Used for datasets integration,
      scaling, and dimensional reduction.
    'sd:layout':
      advanced: true

  gex_dimensionality:
    type: int?
    default: 50
    label: "Number of principal components to use in GEX UMAP projection and clustering (from 1 to 50)"
    doc: |
      Number of principal components to use in GEX UMAP projection and clustering
      (from 1 to 50). Use Elbow plot to adjust this parameter.
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
    type: int?
    default: 1000
    label: "Include cells where at least this many ATAC UMIs (transcripts) are detected"
    doc: |
      Include cells where at least this many ATAC UMIs (transcripts) are detected.
    'sd:layout':
      advanced: true

  maximum_nucl_signal:
    type: float?
    default: 4
    label: "Include cells with the nucleosome signal not bigger than this value"
    doc: |
      Include cells with the nucleosome signal not bigger than this value.
      Nucleosome signal quantifies the approximate ratio of mononucleosomal
      to nucleosome-free fragments.
    'sd:layout':
      advanced: true

  minimum_frip:
    type: float?
    default: 0.15
    label: "Include cells with the FRiP not lower than this value"
    doc: |
      Include cells with the FRiP not lower than this value.
    'sd:layout':
      advanced: true

  maximum_blacklisted_ratio:
    type: float?
    default: 0.05
    label: "Include cells with the ratio of reads in genomic blacklist regions not bigger than this value"
    doc: |
      Include cells with the ratio of reads in genomic blacklist regions
      not bigger than this value.
    'sd:layout':
      advanced: true

  call_peaks:
    type: boolean?
    default: false
    label: "Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count"
    doc: |
      Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.
    'sd:layout':
      advanced: true

  atac_dimensionality:
    type: int?
    default: 50
    label: "Number of principal components to use in ATAC UMAP projection and clustering (from 1 to 50)"
    doc: |
      Number of principal components to use in ATAC UMAP projection and clustering
      (from 1 to 50).
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

  threads:
    type: int?
    default: 4
    label: "Threads number to use"
    doc: |
      Threads number
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
        tab: 'QC (not filtered)'
        Caption: 'Number of cells per dataset (not filtered)'
  
  raw_gex_umi_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_gex_umi_dnst_plot_png
    label: "GEX UMI density per cell (not filtered)"
    doc: |
      GEX UMI density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
        Caption: 'GEX vs ATAC UMIs per cell correlation (not filtered)'

  raw_frg_len_hist_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_frg_len_hist_png
    label: "Fragments Length Histogram (not filtered)"
    doc: |
      Fragments Length Histogram (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
        Caption: 'Density of transcripts mapped to mitochondrial genes per cell (not filtered)'

  raw_nvlt_score_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/raw_nvlt_score_dnst_plot_png
    label: "Novelty score density per cell (not filtered)"
    doc: |
      Novelty score density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (not filtered)'
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
        tab: 'QC (not filtered)'
        Caption: 'QC metrics densities per cell (not filtered)'

  raw_qc_mtrcs_tsv:
    type: File?
    outputSource: seurat_wnn_cluster/raw_qc_mtrcs_tsv
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      TSV format
    'sd:visualPlugins':
    - scatter:
        tab: 'QC (interactive)'
        Title: 'GEX UMI density per cell (not filtered)'
        xAxisTitle: 'GEX UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#2e8b57"]
        height: 300
        data: [$1, $2]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Gene density per cell (not filtered)'
        xAxisTitle: 'Genes per cell'
        yAxisTitle: 'Density'
        colors: ["#191970"]
        height: 300
        data: [$3, $4]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of transcripts mapped to mitochondrial genes per cell (not filtered)'
        xAxisTitle: 'Percentage of transcripts mapped to mitochondrial genes per cell'
        yAxisTitle: 'Density'
        colors: ["#ff0000"]
        height: 300
        data: [$5, $6]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Novelty score density per cell (not filtered)'
        xAxisTitle: 'log10 Genes / log10 UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#ffd700"]
        height: 300
        data: [$7, $8]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'ATAC UMI density per cell (not filtered)'
        xAxisTitle: 'ATAC UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#7cfc00"]
        height: 300
        data: [$9, $10]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Peak density per cell (not filtered)'
        xAxisTitle: 'Peaks per cell'
        yAxisTitle: 'Density'
        colors: ["#e9967a"]
        height: 300
        data: [$11, $12]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of mononucleosomal to nucleosome-free fragments ratio per cell (not filtered)'
        xAxisTitle: 'Mononucleosomal to nucleosome-free fragments ratio per cell'
        yAxisTitle: 'Density'
        colors: ["#00bfff"]
        height: 300
        data: [$13, $14]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of FRiP per cell (not filtered)'
        xAxisTitle: 'Fraction of reads in peaks per cell'
        yAxisTitle: 'Density'
        colors: ["#0000ff"]
        height: 300
        data: [$15, $16]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of fraction of reads within blacklisted regions per cell (not filtered)'
        xAxisTitle: 'Fraction of reads within blacklisted regions per cell'
        yAxisTitle: 'Density'
        colors: ["#ff1493"]
        height: 300
        data: [$17, $18]

  fltr_cell_count_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_cell_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (filtered)'
        Caption: 'Number of cells per dataset (filtered)'
  
  fltr_gex_umi_dnst_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_gex_umi_dnst_plot_png
    label: "GEX UMI density per cell (filtered)"
    doc: |
      GEX UMI density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
        Caption: 'GEX vs ATAC UMIs per cell correlation (filtered)'

  fltr_frg_len_hist_png:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_frg_len_hist_png
    label: "Fragments Length Histogram (filtered)"
    doc: |
      Fragments Length Histogram (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
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
        tab: 'QC (filtered)'
        Caption: 'QC metrics densities per cell (filtered)'

  fltr_qc_mtrcs_tsv:
    type: File?
    outputSource: seurat_wnn_cluster/fltr_qc_mtrcs_tsv
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      TSV format
    'sd:visualPlugins':
    - scatter:
        tab: 'QC (interactive)'
        Title: 'GEX UMI density per cell (filtered)'
        xAxisTitle: 'GEX UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#2e8b57"]
        height: 300
        data: [$1, $2]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Gene density per cell (filtered)'
        xAxisTitle: 'Genes per cell'
        yAxisTitle: 'Density'
        colors: ["#191970"]
        height: 300
        data: [$3, $4]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of transcripts mapped to mitochondrial genes per cell (filtered)'
        xAxisTitle: 'Percentage of transcripts mapped to mitochondrial genes per cell'
        yAxisTitle: 'Density'
        colors: ["#ff0000"]
        height: 300
        data: [$5, $6]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Novelty score density per cell (filtered)'
        xAxisTitle: 'log10 Genes / log10 UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#ffd700"]
        height: 300
        data: [$7, $8]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'ATAC UMI density per cell (filtered)'
        xAxisTitle: 'ATAC UMIs per cell'
        yAxisTitle: 'Density'
        colors: ["#7cfc00"]
        height: 300
        data: [$9, $10]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Peak density per cell (filtered)'
        xAxisTitle: 'Peaks per cell'
        yAxisTitle: 'Density'
        colors: ["#e9967a"]
        height: 300
        data: [$11, $12]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of mononucleosomal to nucleosome-free fragments ratio per cell (filtered)'
        xAxisTitle: 'Mononucleosomal to nucleosome-free fragments ratio per cell'
        yAxisTitle: 'Density'
        colors: ["#00bfff"]
        height: 300
        data: [$13, $14]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of FRiP per cell (filtered)'
        xAxisTitle: 'Fraction of reads in peaks per cell'
        yAxisTitle: 'Density'
        colors: ["#0000ff"]
        height: 300
        data: [$15, $16]
    - scatter:
        tab: 'QC (interactive)'
        Title: 'Density of fraction of reads within blacklisted regions per cell (filtered)'
        xAxisTitle: 'Fraction of reads within blacklisted regions per cell'
        yAxisTitle: 'Density'
        colors: ["#ff1493"]
        height: 300
        data: [$17, $18]

  ntgr_elbow_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_elbow_plot_png
    label: "Elbow plot from PCA of filtered integrated/scaled datasets"
    doc: |
      Elbow plot from PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality evaluation'
        Caption: 'Elbow plot from PCA of filtered integrated/scaled datasets'

  ntgr_pca_plot_png:
    type: File?
    outputSource: seurat_wnn_cluster/ntgr_pca_plot_png
    label: "PCA of filtered integrated/scaled datasets"
    doc: |
      PCA of filtered integrated/scaled datasets.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Dimensionality evaluation'
        Caption: 'PCA of filtered integrated/scaled datasets'
  
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
        tab: 'Clustering'
        Caption: 'Clustered UMAP projected PCA of filtered GEX datasets'

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
        tab: 'Clustering'
        Caption: 'Clustered UMAP projected LSI of filtered ATAC datasets'

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
        tab: 'Clustering'
        Caption: 'Clustered UMAP projected WNN'

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
        tab: 'Gene expression'
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
        tab: 'Gene expression'
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
        tab: 'Gene expression'
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
            #!/bin/bash
            if [ -f "$0" ] && [ -f "$1" ]; then
              echo -e "barcode\tgex_genotype" > gex_metadata.tsv
              cat "$0" | grep -v "barcode" | awk 'BEGIN {OFS = "\t"} ; {if ($2 == "singlet") print $1, $3; else print $1, $2}' >> gex_metadata.tsv
              echo "atac_genotype" > atac_metadata.tsv
              cat "$1" | grep -v "barcode" | awk 'BEGIN {OFS = "\t"} ; {if ($2 == "singlet") print $3; else print $2}' >> atac_metadata.tsv
              paste gex_metadata.tsv atac_metadata.tsv > extra_metadata.tsv
              rm -f gex_metadata.tsv atac_metadata.tsv
              head extra_metadata.tsv
            else
              exit 0
            fi
          inputBinding:
            position: 5
        gex_genotype_cluster_tsv_file:
          type: File?
          inputBinding:
            position: 6
        atac_genotype_cluster_tsv_file:
          type: File?
          inputBinding:
            position: 7
      outputs:
        metadata_file:
          type: File?
          outputBinding:
            glob: "extra_metadata.tsv"
      baseCommand: ["bash", "-c"]
    in:
      gex_genotype_cluster_tsv_file: gex_genotype_cluster_tsv_file
      atac_genotype_cluster_tsv_file: atac_genotype_cluster_tsv_file
    out:
    - metadata_file

  seurat_wnn_cluster:
    run: ../tools/seurat-wnn-cluster.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      atac_fragments_file: atac_fragments_file
      annotation_gtf_file: annotation_gtf_file
      blacklisted_regions_file: blacklisted_regions_file
      barcodes_data: barcodes_data
      metadata_file: prepare_metadata/metadata_file
      gex_minimum_cells: gex_minimum_cells
      gex_minimum_features: gex_minimum_features
      gex_maximum_features: gex_maximum_features
      gex_minimum_umis: gex_minimum_umis
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      minimum_novelty_score: minimum_novelty_score
      atac_minimum_cells: atac_minimum_cells
      atac_minimum_umis: atac_minimum_umis
      maximum_nucl_signal: maximum_nucl_signal
      minimum_frip: minimum_frip
      maximum_blacklisted_ratio: maximum_blacklisted_ratio
      call_peaks: call_peaks
      gex_selected_features:
        source: gex_selected_features
        valueFrom: $(split_features(self))
      gex_high_var_features_count: gex_high_var_features_count
      gex_dimensionality: gex_dimensionality
      atac_dimensionality: atac_dimensionality
      resolution:
        source: resolution
        valueFrom: $(split_numbers(self))
      export_pdf_plots:
        default: false
      export_rds_data:
        default: true
      threads: threads
    out:
    - raw_cell_count_plot_png
    - raw_gex_umi_dnst_plot_png
    - raw_atac_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_peak_dnst_plot_png
    - raw_bl_cnts_dnst_plot_png
    - raw_gex_atac_umi_corr_plot_png
    - raw_frg_len_hist_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_perc_dnst_plot_png
    - raw_nvlt_score_dnst_plot_png
    - raw_qc_mtrcs_plot_png
    - raw_qc_mtrcs_tsv
    - fltr_cell_count_plot_png
    - fltr_gex_umi_dnst_plot_png
    - fltr_atac_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_peak_dnst_plot_png
    - fltr_bl_cnts_dnst_plot_png
    - fltr_gex_atac_umi_corr_plot_png
    - fltr_frg_len_hist_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_perc_dnst_plot_png
    - fltr_nvlt_score_dnst_plot_png
    - fltr_qc_mtrcs_plot_png
    - fltr_qc_mtrcs_tsv
    - ntgr_elbow_plot_png
    - ntgr_pca_plot_png
    - clst_gex_umap_res_plot_png
    - clst_atac_umap_res_plot_png
    - clst_wnn_umap_res_plot_png
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