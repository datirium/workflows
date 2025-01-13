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
          var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var get_source_column = function(prefix, reduction, resolution) {
          if (reduction == "RNA" && resolution != null) {
            return prefix + "rna_res." + resolution;
          } else if (reduction == "WNN" && resolution != null) {
            return prefix + "wsnn_res." + resolution;
          } else {
            return null;
          }
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-rna-filter.cwl"
  - "sc-multiome-filter.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-atac-reduce.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-ctype-assign.cwl"
  - "sc-rna-load-rhapsody.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  sc_reference_model:
  - "sc-ctype-assign.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with Filtered RNA-Seq Datasets"
    doc: |
      Analysis that includes filtered
      single-cell data and was run through
      "Single-Cell RNA-Seq Filtering Analysis"
      or "Single-Cell Multiome ATAC-Seq and
      RNA-Seq Filtering Analysis" at any of
      the processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File?
    secondaryFiles:
    - .tbi
    label: "Cell Ranger RNA+ATAC Sample (optional)"
    doc: |
      Any "Cell Ranger RNA+ATAC Sample" for 
      generating ATAC fragments coverage plots
      over the genes of interest. This sample
      can be obtained from the "Cell Ranger
      Count (RNA+ATAC)" or "Cell Ranger
      Aggregate (RNA+ATAC)" pipelines.
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

  reference_data_rds:
    type: File?
    label: "Reference Single-cell Analysis with Assigned Cell Types (for a custom reference attach files below)"
    doc: |
      Analysis that includes single-cell
      RNA-Seq datasets run through the
      "Single-Cell Manual Cell Type
      Assignment" pipeline based on the
      RNA or WNN clustering results.
    "sd:upstreamSource": "sc_reference_model/reference_data_rds"
    "sd:localLabel": true

  reference_data_index:
    type: File?
    "sd:upstreamSource": "sc_reference_model/reference_data_index"

  query_reduction:
    type:
    - "null"
    - type: enum
      symbols:
      - "RNA"
      - "ATAC"
      - "WNN"
    default: "RNA"
    "sd:upstreamSource": "sc_reference_model/query_reduction"

  query_resolution:
    type: float?
    "sd:upstreamSource": "sc_reference_model/query_resolution"

  custom_reference_data_rds:
    type: File?
    label: "Custom reference Seurat Object (optional)"
    doc: |
      RDS file to load the reference Seurat object from.
      This file can be downloaded as ref.Rds from the
      https://azimuth.hubmapconsortium.org/references/

  custom_reference_data_index:
    type: File?
    label: "Custom reference Annoy Index (optional)"
    doc: |
      Annoy index file for the provided reference RDS file.
      This file can be downloaded as idx.annoy from the
      https://azimuth.hubmapconsortium.org/references/

  custom_reference_source_column:
    type: string?
    label: "Custom reference annotation column to select cell types (optional)"
    doc: |
      Column from the metadata of the custom reference
      Seurat object to select the reference annotations.

  minimum_confidence_score:
    type: float?
    default: 0.75
    label: "Minimum prediction confidence score"
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
    default: 0.75
    label: "Minimum prediction mapping score"
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
    default: true
    label: "Find gene markers"
    doc: |
      Identify upregulated genes in each
      predicted cell type compared to all
      other cells. Include only genes that
      are expressed in at least 10% of the
      cells coming from either current cell
      type or from all other cell types
      together. Exclude cells with
      log2FoldChange values less than 0.25.
      Use Wilcoxon Rank Sum test to
      calculate P-values. Keep only genes
      with P-values lower than 0.01. Adjust
      P-values for multiple comparisons
      using Bonferroni correction.
      Default: true

  identify_diff_peaks:
    type: boolean?
    default: false
    label: "Find peak markers"
    doc: |
      Identify differentially accessible
      peaks in each predicted cell type
      compared to all other cells. Include
      only peaks that are present in at
      least 5% of the cells coming from
      either current cell type or from all
      other cell types together. Exclude
      cells with log2FoldChange values less
      than 0.25. Use logistic regression
      framework to calculate P-values. Keep
      only peaks with P-values lower than
      0.01. Adjust P-values for multiple
      comparisons using Bonferroni
      correction.
      Default: false

  genes_of_interest:
    type: string?
    default: null
    label: "Genes of interest"
    doc: |
      Comma or space separated list of genes
      of interest to visualize expression and
      to generate ATAC fragments coverage plots.
      Ignored if "Cell Ranger RNA+ATAC Sample
      (optional)" input is not provided.
      Default: None

  genesets_data:
    type: File?
    label: "GMT file for calculating average expression levels per gene set (optional)"
    doc: |
      Path to the GMT file for calculating average expression levels
      (module scores) per gene set. This file can be downloaded from
      the Molecular Signatures Database (MSigDB) following the link
      https://www.gsea-msigdb.org/gsea/msigdb.
      Default: do not calculate gene set expression scores.

  export_loupe_data:
    type: boolean?
    default: false
    label: "Save raw counts to Loupe file. I confirm that data is generated by 10x technology and accept the EULA available at https://10xgen.com/EULA"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false
    "sd:layout":
      advanced: true

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
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
    label: "Plots color theme"
    doc: |
      Color theme for all plots saved
      as PNG files.
      Default: classic
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "4"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  ref_cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/ref_cell_cnts_gr_ctyp_plot_png
    label: "Number of cells per cell type (all reference cells)"
    doc: |
      Number of cells per cell type.
      All reference cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Reference"
        Caption: "Number of cells per cell type (all reference cells)"

  ref_umap_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/ref_umap_gr_ctyp_plot_png
    label: "Reference UMAP colored by cell type (all reference cells)"
    doc: |
      Reference UMAP colored by cell type.
      All reference cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Reference"
        Caption: "Reference UMAP colored by cell type (all reference cells)"

  cell_cnts_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/cell_cnts_gr_ctyp_plot_png
    label: "Number of cells per cell type (all query cells)"
    doc: |
      Number of cells per cell type.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Number of cells per cell type (all query cells)"

  qc_mtrcs_dnst_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/qc_mtrcs_dnst_gr_ctyp_plot_png
    label: "Distribution of QC metrics per cell colored by cell type (all query cells)"
    doc: |
      Distribution of QC metrics per
      cell colored by cell type.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Distribution of QC metrics per cell colored by cell type (all query cells)"

  umap_cnf_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_cnf_plot_png
    label: "Projected UMAP colored by prediction confidence score (all query cells)"
    doc: |
      Projected UMAP colored by
      prediction confidence score.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Projected UMAP colored by prediction confidence score (all query cells)"

  umap_map_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_map_plot_png
    label: "Projected UMAP colored by prediction mapping score (all query cells)"
    doc: |
      Projected UMAP colored by
      prediction mapping score.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Projected UMAP colored by prediction mapping score (all query cells)"

  gene_umi_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/gene_umi_gr_cnf_spl_ctyp_plot_png
    label: "Genes vs RNA reads per cell colored by prediction confidence score (split by cell type, all query cells)"
    doc: |
      Genes vs RNA reads per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Genes vs RNA reads per cell colored by prediction confidence score (split by cell type, all query cells)"

  gene_umi_gr_map_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/gene_umi_gr_map_spl_ctyp_plot_png
    label: "Genes vs RNA reads per cell colored by prediction mapping score (split by cell type, all query cells)"
    doc: |
      Genes vs RNA reads per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Genes vs RNA reads per cell colored by prediction mapping score (split by cell type, all query cells)"

  umi_mito_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/umi_mito_gr_cnf_spl_ctyp_plot_png
    label: "RNA reads vs mitochondrial % per cell colored by prediction confidence score (split by cell type, all query cells)"
    doc: |
      RNA reads vs mitochondrial % per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs mitochondrial % per cell colored by prediction confidence score (split by cell type, all query cells)"

  umi_mito_gr_map_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/umi_mito_gr_map_spl_ctyp_plot_png
    label: "RNA reads vs mitochondrial % per cell colored by prediction mapping score (split by cell type, all query cells)"
    doc: |
      RNA reads vs mitochondrial % per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs mitochondrial % per cell colored by prediction mapping score (split by cell type, all query cells)"

  tss_frgm_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/tss_frgm_gr_cnf_spl_ctyp_plot_png
    label: "TSS enrichment score vs ATAC fragments in peaks per cell colored by prediction confidence score (split by cell type, all query cells)"
    doc: |
      TSS enrichment score vs ATAC fragments
      in peaks per cell.
      All query cells; split by cell type;
      colored by prediction confidence score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "TSS enrichment score vs ATAC fragments in peaks per cell colored by prediction confidence score (split by cell type, all query cells)"

  tss_frgm_gr_map_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/tss_frgm_gr_map_spl_ctyp_plot_png
    label: "TSS enrichment score vs ATAC fragments in peaks per cell colored by prediction mapping score (split by cell type, all query cells)"
    doc: |
      TSS enrichment score vs ATAC fragments
      in peaks per cell.
      All query cells; split by cell type;
      colored by prediction mapping score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "TSS enrichment score vs ATAC fragments in peaks per cell colored by prediction mapping score (split by cell type, all query cells)"

  rna_atac_cnts_gr_cnf_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/rna_atac_cnts_gr_cnf_spl_ctyp_plot_png
    label: "RNA reads vs ATAC fragments in peaks per cell colored by prediction confidence score (split by cell type, all query cells)"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      All query cells; split by cell type; colored
      by prediction confidence score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs ATAC fragments in peaks per cell colored by prediction confidence score (split by cell type, all query cells)"

  rna_atac_cnts_gr_map_spl_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/rna_atac_cnts_gr_map_spl_ctyp_plot_png
    label: "RNA reads vs ATAC fragments in peaks per cell colored by prediction mapping score (split by cell type, all query cells)"
    doc: |
      RNA reads vs ATAC fragments in peaks per cell.
      All query cells; split by cell type; colored
      by prediction mapping score.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "RNA reads vs ATAC fragments in peaks per cell colored by prediction mapping score (split by cell type, all query cells)"

  rnadbl_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/rnadbl_gr_ctyp_plot_png
    label: "Percentage of RNA doublets per cell type (all query cells)"
    doc: |
      Percentage of RNA doublets
      per cell type.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA doublets per cell type (all query cells)"

  atacdbl_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/atacdbl_gr_ctyp_plot_png
    label: "Percentage of ATAC doublets per cell type (all query cells)"
    doc: |
      Percentage of ATAC doublets
      per cell type.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of ATAC doublets per cell type (all query cells)"

  vrlpdbl_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/vrlpdbl_gr_ctyp_plot_png
    label: "Percentage of RNA and ATAC doublets per cell type (all query cells)"
    doc: |
      Percentage of RNA and ATAC
      doublets per cell type.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "Percentage of RNA and ATAC doublets per cell type (all query cells)"

  umap_gr_ctyp_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ctyp_plot_png
    label: "Projected UMAP colored by cell type (filtered query cells)"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "Projected UMAP colored by cell type (filtered query cells)"

  umap_gr_ctyp_spl_ph_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ctyp_spl_ph_png
    label: "Projected UMAP colored by cell type (split by cell cycle phase, optionally downsampled filtered query cells)"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by cell
      cycle phase; downsampled to the
      smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "Projected UMAP colored by cell type (split by cell cycle phase, optionally downsampled filtered query cells)"

  cmp_gr_ph_spl_ctyp_png:
    type: File?
    outputSource: rna_azimuth/cmp_gr_ph_spl_ctyp_png
    label: "Composition plot colored by cell cycle phase (split by cell type, optionally downsampled filtered query cells)"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by cell type;
      downsampled to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by cell type"
        Caption: "Composition plot colored by cell cycle phase (split by cell type, optionally downsampled filtered query cells)"

  umap_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ctyp_spl_idnt_plot_png
    label: "Projected UMAP colored by cell type (split by dataset, downsampled filtered query cells)"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Projected UMAP colored by cell type (split by dataset, downsampled filtered query cells)"

  cmp_gr_ctyp_spl_idnt_plot_png:
    type: File?
    outputSource: rna_azimuth/cmp_gr_ctyp_spl_idnt_plot_png
    label: "Composition plot colored by cell type (split by dataset, downsampled filtered query cells)"
    doc: |
      Composition plot colored by cell type.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cell type (split by dataset, downsampled filtered query cells)"

  umap_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ph_spl_idnt_plot_png
    label: "Projected UMAP colored by cell cycle phase (split by dataset, downsampled filtered query cells)"
    doc: |
      Projected UMAP colored by cell cycle phase.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Projected UMAP colored by cell cycle phase (split by dataset, downsampled filtered query cells)"

  cmp_gr_ph_spl_idnt_plot_png:
    type: File?
    outputSource: rna_azimuth/cmp_gr_ph_spl_idnt_plot_png
    label: "Composition plot colored by cell cycle phase (split by dataset, downsampled filtered query cells)"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by dataset;
      downsampled to the smallest dataset.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by dataset"
        Caption: "Composition plot colored by cell cycle phase (split by dataset, downsampled filtered query cells)"

  umap_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ctyp_spl_cnd_plot_png
    label: "Projected UMAP colored by cell type (split by grouping condition, downsampled filtered query cells)"
    doc: |
      Projected UMAP colored by cell type.
      Filtered query cells; split by grouping
      condition; first downsampled to the
      smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Projected UMAP colored by cell type (split by grouping condition, downsampled filtered query cells)"

  cmp_gr_ctyp_spl_cnd_plot_png:
    type: File?
    outputSource: rna_azimuth/cmp_gr_ctyp_spl_cnd_plot_png
    label: "Composition plot colored by cell type (split by grouping condition, downsampled filtered query cells)"
    doc: |
      Composition plot colored by cell type.
      Filtered query cells; split by grouping
      condition; first downsampled to the
      smallest dataset, then downsampled to
      the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cell type (split by grouping condition, downsampled filtered query cells)"

  umap_gr_ph_spl_cnd_plot_png:
    type: File?
    outputSource: rna_azimuth/umap_gr_ph_spl_cnd_plot_png
    label: "Projected UMAP colored by cell cycle phase (split by grouping condition, downsampled filtered query cells)"
    doc: |
      Projected UMAP colored by cell cycle phase.
      Filtered query cells; split by grouping
      condition; first downsampled to the smallest
      dataset, then downsampled to the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Projected UMAP colored by cell cycle phase (split by grouping condition, downsampled filtered query cells)"

  cmp_gr_ph_spl_cnd_plot_png:
    type: File?
    outputSource: rna_azimuth/cmp_gr_ph_spl_cnd_plot_png
    label: "Composition plot colored by cell cycle phase (split by grouping condition, downsampled filtered query cells)"
    doc: |
      Composition plot colored by cell cycle phase.
      Filtered query cells; split by grouping condition;
      first downsampled to the smallest dataset, then
      downsampled to the smallest group.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Split by group"
        Caption: "Composition plot colored by cell cycle phase (split by grouping condition, downsampled filtered query cells)"

  xpr_avg_plot_png:
    type: File?
    outputSource: rna_azimuth/xpr_avg_plot_png
    label: "Average gene expression (filtered query cells)"
    doc: |
      Average gene expression.
      Filtered query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Average gene expression (filtered query cells)"

  xpr_dnst_plot_png:
    type: File?
    outputSource: rna_azimuth/xpr_dnst_plot_png
    label: "Gene expression density (filtered query cells)"
    doc: |
      Gene expression density.
      Filtered query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Gene expression density (filtered query cells)"

  xpr_per_cell_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: rna_azimuth/xpr_per_cell_plot_png
    label: "Projected UMAP colored by gene expression (per gene, filtered query cells)"
    doc: |
      Projected UMAP colored by gene expression.
      Filtered query cells; all genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (expression)"
        Caption: "Projected UMAP colored by gene expression (per gene, filtered query cells)"

  cvrg_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputSource: rna_azimuth/cvrg_plot_png
    label: "ATAC fragment coverage (per gene, filtered query cells)"
    doc: |
      ATAC fragment coverage.
      Filtered query cells;
      all genes of interest.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Genes of interest (coverage)"
        Caption: "ATAC fragment coverage (per gene, filtered query cells)"

  gse_per_cell_plot_png:
    type: File?
    outputSource: rna_azimuth/gse_per_cell_plot_png
    label: "Projected UMAP colored by gene set expression score (all query cells)"
    doc: |
      Projected UMAP colored by
      gene set expression score.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Projected UMAP colored by gene set expression score (all query cells)"

  gse_avg_plot_png:
    type: File?
    outputSource: rna_azimuth/gse_avg_plot_png
    label: "Average gene set expression score (all query cells)"
    doc: |
      Average gene set expression score.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Average gene set expression score (all query cells)"

  gse_dnst_plot_png:
    type: File?
    outputSource: rna_azimuth/gse_dnst_plot_png
    label: "Gene set expression score density (all query cells)"
    doc: |
      Gene set expression score density.
      All query cells.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene sets of interest (expression)"
        Caption: "Gene set expression score density (all query cells)"

  xpr_htmp_plot_png:
    type: File?
    outputSource: rna_azimuth/xpr_htmp_plot_png
    label: "Gene expression heatmap (top gene markers, filtered query cells)"
    doc: |
      Gene expression heatmap from
      the filtered query cells.
      Top gene markers.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Gene markers heatmap"
        Caption: "Gene expression heatmap (top gene markers, filtered query cells)"

  xpr_htmp_tsv:
    type: File?
    outputSource: rna_azimuth/xpr_htmp_tsv
    label: "Gene expression heatmap (top gene markers, filtered query cells)"
    doc: |
      Gene expression heatmap from
      the filtered query cells.
      Top gene markers.
      TSV format.

  gene_markers_tsv:
    type: File?
    outputSource: rna_azimuth/gene_markers_tsv
    label: "Gene markers (filtered query cells)"
    doc: |
      Gene markers from the filtered
      query cells.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Gene markers table"
        Title: "Gene markers (filtered query cells)"

  peak_markers_tsv:
    type: File?
    outputSource: rna_azimuth/peak_markers_tsv
    label: "Peak markers (filtered query cells)"
    doc: |
      Peak markers from the filtered
      query cells.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Peak markers table"
        Title: "Peak markers (filtered query cells)"

  ucsc_cb_html_data:
    type: Directory?
    outputSource: rna_azimuth/ucsc_cb_html_data
    label: "UCSC Cell Browser (data)"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputSource: rna_azimuth/ucsc_cb_html_file
    label: "UCSC Cell Browser"
    doc: |
      UCSC Cell Browser html index.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  seurat_data_rds:
    type: File
    outputSource: rna_azimuth/seurat_data_rds
    label: "Seurat object in RDS format"
    doc: |
      Seurat object.
      RDS format.

  seurat_data_scope:
    type: File?
    outputSource: rna_azimuth/seurat_data_scope
    label: "Seurat object in SCope compatible loom format"
    doc: |
      Seurat object.
      SCope compatible.
      Loom format.

  seurat_rna_data_cloupe:
    type: File?
    outputSource: rna_azimuth/seurat_rna_data_cloupe
    label: "Seurat object in Loupe format"
    doc: |
      Seurat object.
      RNA counts.
      Loupe format.

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_report_html_file:
    type: File?
    outputSource: rna_azimuth/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  rna_azimuth_stdout_log:
    type: File
    outputSource: rna_azimuth/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the rna_azimuth step.

  rna_azimuth_stderr_log:
    type: File
    outputSource: rna_azimuth/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the rna_azimuth step.


steps:

  rna_azimuth:
    run: ../tools/sc-rna-azimuth.cwl
    in:
      query_data_rds: query_data_rds
      reference_data_rds:
        source: [reference_data_rds, custom_reference_data_rds, custom_reference_data_index, custom_reference_source_column]
        valueFrom: |
          ${
            if (
              self[1] != null && self[1].class == "File" &&
              self[2] != null && self[2].class == "File" &&
              self[3] != null
            ){
              return self[1];
            } else {
              return self[0];
            }
          }
      reference_data_index:
        source: [reference_data_index, custom_reference_data_rds, custom_reference_data_index, custom_reference_source_column]
        valueFrom: |
          ${
            if (
              self[1] != null && self[1].class == "File" &&
              self[2] != null && self[2].class == "File" &&
              self[3] != null
            ){
              return self[2];
            } else {
              return self[0];
            }
          }
      reference_source_column:
        source: [query_reduction, query_resolution, custom_reference_data_rds, custom_reference_data_index, custom_reference_source_column]
        valueFrom: |
          ${
            if (
              self[2] != null && self[2].class == "File" &&
              self[3] != null && self[3].class == "File" &&
              self[4] != null
            ){
              return self[4];
            } else {
              return get_source_column("custom_", self[0], self[1]);
            }
          }
      minimum_confidence_score: minimum_confidence_score
      minimum_mapping_score: minimum_mapping_score
      atac_fragments_file: atac_fragments_file
      genes_of_interest:
        source: genes_of_interest
        valueFrom: $(split_features(self))
      genesets_data: genesets_data
      identify_diff_genes: identify_diff_genes
      identify_diff_peaks: identify_diff_peaks
      rna_minimum_logfc:
        default: 0.25
      rna_minimum_pct:
        default: 0.1
      atac_minimum_logfc:
        default: 0.25
      atac_minimum_pct:
        default: 0.05
      only_positive_diff_genes:
        default: true
      rna_test_to_use: 
        default: wilcox
      atac_test_to_use:
        default: LR
      verbose:
        default: true
      export_ucsc_cb:
        default: true
      export_scope_data:
        default: true
      export_loupe_data: export_loupe_data
      export_pdf_plots:
        default: true
      color_theme: color_theme
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - ref_cell_cnts_gr_ctyp_plot_png
    - ref_umap_gr_ctyp_plot_png
    - cell_cnts_gr_ctyp_plot_png
    - qc_mtrcs_dnst_gr_ctyp_plot_png
    - umap_cnf_plot_png
    - umap_map_plot_png
    - gene_umi_gr_cnf_spl_ctyp_plot_png
    - gene_umi_gr_map_spl_ctyp_plot_png
    - umi_mito_gr_cnf_spl_ctyp_plot_png
    - umi_mito_gr_map_spl_ctyp_plot_png
    - tss_frgm_gr_cnf_spl_ctyp_plot_png
    - tss_frgm_gr_map_spl_ctyp_plot_png
    - rna_atac_cnts_gr_cnf_spl_ctyp_plot_png
    - rna_atac_cnts_gr_map_spl_ctyp_plot_png
    - rnadbl_gr_ctyp_plot_png
    - atacdbl_gr_ctyp_plot_png
    - vrlpdbl_gr_ctyp_plot_png
    - umap_gr_ctyp_plot_png
    - umap_gr_ctyp_spl_ph_png
    - cmp_gr_ph_spl_ctyp_png
    - umap_gr_ctyp_spl_idnt_plot_png
    - cmp_gr_ctyp_spl_idnt_plot_png
    - umap_gr_ph_spl_idnt_plot_png
    - cmp_gr_ph_spl_idnt_plot_png
    - umap_gr_ctyp_spl_cnd_plot_png
    - cmp_gr_ctyp_spl_cnd_plot_png
    - umap_gr_ph_spl_cnd_plot_png
    - cmp_gr_ph_spl_cnd_plot_png
    - gse_per_cell_plot_png
    - gse_avg_plot_png
    - gse_dnst_plot_png
    - xpr_avg_plot_png
    - xpr_dnst_plot_png
    - xpr_per_cell_plot_png
    - cvrg_plot_png
    - xpr_htmp_plot_png
    - xpr_htmp_tsv
    - gene_markers_tsv
    - peak_markers_tsv
    - all_plots_pdf
    - ucsc_cb_html_data
    - ucsc_cb_html_file
    - seurat_data_rds
    - seurat_data_scope
    - seurat_rna_data_cloupe
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - rna_azimuth/all_plots_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell RNA-Seq Reference Mapping"
s:name: "Single-Cell RNA-Seq Reference Mapping"
s:alternateName: "Single-Cell RNA-Seq Reference Mapping"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-azimuth.cwl
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
  Single-Cell RNA-Seq Reference Mapping

  Uses Azimuth R package to assign identities to cells based on the
  reference annotation from the results of the "Single-Cell Manual
  Cell Type Assignment" pipeline. Alternatively, custom reference
  models can be downloaded from the
  https://azimuth.hubmapconsortium.org/ website. This workflow can
  be run with the outputs of the following pipelines: "Single-Cell
  RNA-Seq Filtering Analysis", "Single-Cell Multiome ATAC-Seq and
  RNA-Seq Filtering Analysis", "Single-Cell RNA-Seq Dimensionality
  Reduction Analysis", "Single-Cell RNA-Seq Cluster Analysis", and
  "Single-Cell WNN Cluster Analysis". It can also be used with the
  outputs of: "Single-Cell ATAC-Seq Dimensionality Reduction Analysis",
  "Single-Cell ATAC-Seq Cluster Analysis", "Single-Cell Manual Cell
  Type Assignment" pipelines if these were part of the multiome data
  analysis. The results of this workflow are compatible with any
  single cell pipeline normally used after the "Single-Cell RNA-Seq
  Filtering Analysis" or "Single-Cell Multiome ATAC-Seq and RNA-Seq
  Filtering Analysis" pipelines, depending on the preceding analysis
  step. In other words, this pipeline predicts cell types for
  high-quality cells without impacting subsequent data analysis steps.