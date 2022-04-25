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
      Path to the compressed folder with feature-barcode matrix from Cell Ranger Count/Aggregate
      experiment in MEX format.
    'sd:upstreamSource': "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  aggregation_metadata:
    type: File?
    label: "Cell Ranger Count/Aggregate Experiment"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.
    'sd:upstreamSource': "sc_rnaseq_sample/aggregation_metadata"
    'sd:localLabel': true

  grouping_data:
    type: File?
    label: "Optional TSV/CSV file to define datasets grouping with 'library_id' and 'condition' columns. Rows order should correspond to the aggregation metadata."
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column - 'library_id'
      with the values and order that correspond to the 'library_id' column from the
      '--identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  barcodes_data:
    type: File?
    label: "Optional headerless TSV/CSV file with the list of barcodes to select cells of interest (one barcode per line)"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells.
      Default: use all cells.

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
    label: "Include cells where at least this many UMI (transcripts) are detected"
    doc: |
      Include cells where at least this many UMI (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
      Default: 500 (applied to all datasets)
    'sd:layout':
      advanced: true

  minimum_novelty_score:
    type: string?
    default: "0.8"
    label: "Include cells with the novelty score not lower than this value, calculated as log10(genes)/log10(UMI)"
    doc: |
      Include cells with the novelty score not lower than this value, calculated
      as log10(genes)/log10(UMI). If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)
    'sd:layout':
      advanced: true

  mito_pattern:
    type: string?
    default: "^Mt-"
    label: "Regex pattern to identify mitochondrial genes"
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^Mt-'
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

  parallel_memory_limit:
    type:
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
    - type: enum
      symbols:
      - "32"
    default: "32"
    label: "Maximum vector memory in GB allowed to be used by R"
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Forced to 32 GB
    'sd:layout':
      advanced: true

  threads:
    type:
    - type: enum
      symbols:
      - "1"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
      Forced to 1
    'sd:layout':
      advanced: true


outputs:

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_1_2_qc_mtrcs_pca_plot_png
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
    outputSource: sc_rna_filter/raw_2_3_qc_mtrcs_pca_plot_png
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
    outputSource: sc_rna_filter/raw_cells_count_plot_png
    label: "Number of cells per dataset (not filtered)"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Number of cells per dataset'

  raw_umi_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_umi_dnst_plot_png
    label: "UMI per cell density (not filtered)"
    doc: |
      UMI per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'UMI per cell density'

  raw_gene_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_gene_dnst_plot_png
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
    outputSource: sc_rna_filter/raw_gene_umi_corr_plot_png
    label: "Genes vs UMI per cell correlation (not filtered)"
    doc: |
      Genes vs UMI per cell correlation (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Genes vs UMI per cell correlation'

  raw_mito_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_mito_dnst_plot_png
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
    outputSource: sc_rna_filter/raw_nvlt_dnst_plot_png
    label: "Novelty score per cell density (not filtered)"
    doc: |
      Novelty score per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Novelty score per cell density'

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_qc_mtrcs_dnst_plot_png
    label: "QC metrics per cell density (not filtered)"
    doc: |
      QC metrics per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'QC metrics per cell density'

  raw_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density (not filtered)"
    doc: |
      Split by grouping condition UMI per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition UMI per cell density'

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/raw_gene_dnst_spl_cnd_plot_png
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
    outputSource: sc_rna_filter/raw_mito_dnst_spl_cnd_plot_png
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
    outputSource: sc_rna_filter/raw_nvlt_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the novelty score per cell density (not filtered)"
    doc: |
      Split by grouping condition the novelty score per cell density (not filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Not filtered QC'
        Caption: 'Split by grouping condition the novelty score per cell density'

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_1_2_qc_mtrcs_pca_plot_png
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
    outputSource: sc_rna_filter/fltr_2_3_qc_mtrcs_pca_plot_png
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
    outputSource: sc_rna_filter/fltr_cells_count_plot_png
    label: "Number of cells per dataset (filtered)"
    doc: |
      Number of cells per dataset (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Number of cells per dataset'

  fltr_umi_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_umi_dnst_plot_png
    label: "UMI per cell density (filtered)"
    doc: |
      UMI per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'UMI per cell density'

  fltr_gene_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_gene_dnst_plot_png
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
    outputSource: sc_rna_filter/fltr_gene_umi_corr_plot_png
    label: "Genes vs UMI per cell correlation (filtered)"
    doc: |
      Genes vs UMI per cell correlation (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Genes vs UMI per cell correlation'

  fltr_mito_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_mito_dnst_plot_png
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
    outputSource: sc_rna_filter/fltr_nvlt_dnst_plot_png
    label: "Novelty score per cell density (filtered)"
    doc: |
      Novelty score per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Novelty score per cell density'

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_qc_mtrcs_dnst_plot_png
    label: "QC metrics per cell density (filtered)"
    doc: |
      QC metrics per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'QC metrics per cell density'

  fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_umi_dnst_spl_cnd_plot_png
    label: "Split by grouping condition UMI per cell density (filtered)"
    doc: |
      Split by grouping condition UMI per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition UMI per cell density'

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputSource: sc_rna_filter/fltr_gene_dnst_spl_cnd_plot_png
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
    outputSource: sc_rna_filter/fltr_mito_dnst_spl_cnd_plot_png
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
    outputSource: sc_rna_filter/fltr_nvlt_dnst_spl_cnd_plot_png
    label: "Split by grouping condition the novelty score per cell density (filtered)"
    doc: |
      Split by grouping condition the novelty score per cell density (filtered).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Filtered QC'
        Caption: 'Split by grouping condition the novelty score per cell density'

  seurat_data_rds:
    type: File
    outputSource: sc_rna_filter/seurat_data_rds
    label: "Processed Seurat data in RDS format"
    doc: |
      Processed Seurat data in RDS format

  sc_rna_filter_stdout_log:
    type: File
    outputSource: sc_rna_filter/stdout_log
    label: "stdout log generated by sc_rna_filter step"
    doc: |
      stdout log generated by sc_rna_filter step

  sc_rna_filter_stderr_log:
    type: File
    outputSource: sc_rna_filter/stderr_log
    label: "stderr log generated by sc_rna_filter step"
    doc: |
      stderr log generated by sc_rna_filter step


steps:

  uncompress_feature_bc_matrices:
    doc: |
      Extracts the content of TAR file into a folder
    run: ../tools/tar-extract.cwl
    in:
      file_to_extract: filtered_feature_bc_matrix_folder
    out:
    - extracted_folder

  sc_rna_filter:
    doc: |
      Filters single-cell RNA-Seq datasets based on the common QC metrics
    run: ../tools/sc-rna-filter.cwl
    in:
      feature_bc_matrices_folder: uncompress_feature_bc_matrices/extracted_folder
      aggregation_metadata: aggregation_metadata
      grouping_data: grouping_data
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
      minimum_novelty_score:
        source: minimum_novelty_score
        valueFrom: $(split_numbers(self))
      mito_pattern: mito_pattern
      maximum_mito_perc: maximum_mito_perc
      verbose:
        default: true
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
    - raw_umi_dnst_plot_png
    - raw_gene_dnst_plot_png
    - raw_gene_umi_corr_plot_png
    - raw_mito_dnst_plot_png
    - raw_nvlt_dnst_plot_png
    - raw_qc_mtrcs_dnst_plot_png
    - raw_umi_dnst_spl_cnd_plot_png
    - raw_gene_dnst_spl_cnd_plot_png
    - raw_mito_dnst_spl_cnd_plot_png
    - raw_nvlt_dnst_spl_cnd_plot_png
    - fltr_1_2_qc_mtrcs_pca_plot_png
    - fltr_2_3_qc_mtrcs_pca_plot_png
    - fltr_cells_count_plot_png
    - fltr_umi_dnst_plot_png
    - fltr_gene_dnst_plot_png
    - fltr_gene_umi_corr_plot_png
    - fltr_mito_dnst_plot_png
    - fltr_nvlt_dnst_plot_png
    - fltr_qc_mtrcs_dnst_plot_png
    - fltr_umi_dnst_spl_cnd_plot_png
    - fltr_gene_dnst_spl_cnd_plot_png
    - fltr_mito_dnst_spl_cnd_plot_png
    - fltr_nvlt_dnst_spl_cnd_plot_png
    - seurat_data_rds
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell RNA-Seq Filtering Analysis"
s:name: "Single-cell RNA-Seq Filtering Analysis"
s:alternateName: "Filters single-cell RNA-Seq datasets based on the common QC metrics"

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
  Single-cell RNA-Seq Filtering Analysis
  
  Filters single-cell RNA-Seq datasets based on the common QC metrics.