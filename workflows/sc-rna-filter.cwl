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
  sc_rnaseq_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-aggr.cwl"
  - "https://github.com/datirium/workflows/workflows/single-cell-preprocess-cellranger.cwl"
  - "cellranger-aggr.cwl"
  - "single-cell-preprocess-cellranger.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  filtered_feature_bc_matrix_folder:
    type: File
    label: "Cell Ranger Count/Aggregate Experiment"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger Count/Aggregate
      experiment in MEX format. If multiple locations provided data is assumed to be not
      aggregated (outputs from multiple Cell Ranger Count experiments) and will be merged
      before the analysis.
    'sd:upstreamSource': "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File?
    label: "Cell Ranger Count/Aggregate Experiment"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If --mex points to
      the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. In case of
      using feature-barcode matrices from a single or multiple Cell Ranger Count experiments
      the file with identities should include at least one column - 'library_id', and a row
      with aliases per each experiment from the --mex input. The order of rows should correspond
      to the order of feature-barcode matrices provided in the --mex parameter.
    'sd:upstreamSource': "sc_rnaseq_sample/aggregation_metadata"
    'sd:localLabel': true

  grouping_data:
    type: File?
    label: "TSV/CSV file to define datasets grouping with 'library_id' and 'condition' columns. Rows order should correspond to the aggregation metadata."
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      Default: each dataset is assigned to a separate group.

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
      Include only GEX features detected in at least this many cells. Ignored when
      --mex points to the feature-barcode matrices from the multiple Cell Ranger
      Count experiments.
      Default: 5 (applied to all datasets)
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
      Default: 250 (applied to all datasets)
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

  raw_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_pca_1_2_qc_mtrcs_plot_png
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
    outputSource: sc_gex_filter/raw_pca_2_3_qc_mtrcs_plot_png
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
    outputSource: sc_gex_filter/raw_cell_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Number of cells per dataset (not filtered)'
  
  raw_gex_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_gex_umi_dnst_spl_by_cond_plot_png
    label: "Split by condition GEX UMI density per cell (not filtered)"
    doc: |
      Split by condition GEX UMI density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition GEX UMI density per cell (not filtered)'

  raw_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_gene_dnst_spl_by_cond_plot_png
    label: "Split by condition gene density per cell (not filtered)"
    doc: |
      Split by condition gene density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition gene density per cell (not filtered)'

  raw_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_gene_umi_corr_plot_png
    label: "Genes vs GEX UMIs per cell correlation (not filtered)"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Genes vs GEX UMIs per cell correlation (not filtered)'

  raw_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_mito_perc_dnst_spl_by_cond_plot_png
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered)'

  raw_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_nvlt_score_dnst_spl_by_cond_plot_png
    label: "Split by condition novelty score density per cell (not filtered)"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'Split by condition novelty score density per cell (not filtered)'

  raw_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_gex_filter/raw_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (not filtered)"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 1. Not filtered QC'
        Caption: 'QC metrics densities per cell (not filtered)'

  fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_pca_1_2_qc_mtrcs_plot_png
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
    outputSource: sc_gex_filter/fltr_pca_2_3_qc_mtrcs_plot_png
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
    outputSource: sc_gex_filter/fltr_cell_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Number of cells per dataset (filtered)'
  
  fltr_gex_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_gex_umi_dnst_spl_by_cond_plot_png
    label: "Split by condition GEX UMI density per cell (filtered)"
    doc: |
      Split by condition GEX UMI density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition GEX UMI density per cell (filtered)'

  fltr_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_gene_dnst_spl_by_cond_plot_png
    label: "Split by condition gene density per cell (filtered)"
    doc: |
      Split by condition gene density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition gene density per cell (filtered)'

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_gene_umi_corr_plot_png
    label: "Genes vs GEX UMIs per cell correlation (filtered)"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Genes vs GEX UMIs per cell correlation (filtered)'

  fltr_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_mito_perc_dnst_spl_by_cond_plot_png
    label: "Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered)"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered)'

  fltr_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_nvlt_score_dnst_spl_by_cond_plot_png
    label: "Split by condition novelty score density per cell (filtered)"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'Split by condition novelty score density per cell (filtered)'

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputSource: sc_gex_filter/fltr_qc_mtrcs_plot_png
    label: "QC metrics densities per cell (filtered)"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Step 2. Filtered QC'
        Caption: 'QC metrics densities per cell (filtered)'

  seurat_filtered_data_rds:
    type: File
    outputSource: sc_gex_filter/seurat_filtered_data_rds
    label: "Filtered Seurat data in RDS format"
    doc: |
      Filtered Seurat data in RDS format
      RDS format

  sc_gex_filter_stdout_log:
    type: File
    outputSource: sc_gex_filter/stdout_log
    label: stdout log generated by Single-cell GEX Filter
    doc: |
      stdout log generated by Single-cell GEX Filter

  sc_gex_filter_stderr_log:
    type: File
    outputSource: sc_gex_filter/stderr_log
    label: stderr log generated by Single-cell GEX Filter
    doc: |
      stderr log generated by Single-cell GEX Filter


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

  sc_gex_filter:
    run: ../tools/sc-gex-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/uncompressed
      aggregation_metadata: aggregation_metadata
      grouping_data: grouping_data
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
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      export_pdf_plots:
        default: false
      verbose:
        default: true
      export_h5seurat_data:
        default: false
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
    - raw_pca_1_2_qc_mtrcs_plot_png
    - raw_pca_2_3_qc_mtrcs_plot_png
    - raw_cell_count_plot_png
    - raw_gex_umi_dnst_spl_by_cond_plot_png
    - raw_gene_dnst_spl_by_cond_plot_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_perc_dnst_spl_by_cond_plot_png
    - raw_nvlt_score_dnst_spl_by_cond_plot_png
    - raw_qc_mtrcs_plot_png
    - fltr_pca_1_2_qc_mtrcs_plot_png
    - fltr_pca_2_3_qc_mtrcs_plot_png
    - fltr_cell_count_plot_png
    - fltr_gex_umi_dnst_spl_by_cond_plot_png
    - fltr_gene_dnst_spl_by_cond_plot_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_perc_dnst_spl_by_cond_plot_png
    - fltr_nvlt_score_dnst_spl_by_cond_plot_png
    - fltr_qc_mtrcs_plot_png
    - seurat_filtered_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell GEX Filter"
s:name: "Single-cell GEX Filter"
s:alternateName: "Filters single-cell GEX datasets based on the common QC metrics"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-rna-filter.cwl
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
  Single-cell GEX Filter
  ================================================================
  Filters single-cell GEX datasets based on the common QC metrics.