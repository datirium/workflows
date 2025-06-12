cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
sd:upstream:
  sc_rnaseq_sample:
  - single-cell-preprocess-cellranger.cwl
  - cellranger-aggr.cwl
inputs:
  alias:
    type: string
    label: Experiment short name/alias
    sd:preview:
      position: 1
  raw_feature_bc_matrices_folder:
    type: File
    format: http://edamontology.org/format_3989
    label: scRNA-Seq Cell Ranger Experiment
    doc: Compressed folder with unfiltered feature-barcode matrices
    sd:upstreamSource: sc_rnaseq_sample/raw_feature_bc_matrices_folder
    sd:localLabel: true
  filtered_feature_bc_matrix_folder:
    type: File
    format: http://edamontology.org/format_3989
    label: scRNA-Seq Cell Ranger Experiment
    doc: Compressed folder with filtered feature-barcode matrices
    sd:upstreamSource: sc_rnaseq_sample/filtered_feature_bc_matrix_folder
    sd:localLabel: true
  secondary_analysis_report_folder:
    type: File
    format: http://edamontology.org/format_3989
    label: scRNA-Seq Cell Ranger Experiment
    doc: Compressed folder with secondary analysis results
    sd:upstreamSource: sc_rnaseq_sample/secondary_analysis_report_folder
    sd:localLabel: true
  genelist_file:
    type: File
    format: http://edamontology.org/format_2330
    label: Target genes list. Headerless text file with 1 gene per line
    doc: Target genes list. Headerless text file with 1 gene per line
  expression_threshold:
    type: float?
    default: 0
    label: Expression threshold for displaying target genes on a plot (expression > threshold)
    doc: Expression threshold for displaying target genes on a plot (expression > threshold)
    sd:layout:
      advanced: true
  fdr:
    type: float?
    default: 0.05
    label: FDR cutoff for expression ratio plots
    doc: FDR cutoff for expression ratio plots
    sd:layout:
      advanced: true
  round_counts:
    type: boolean?
    default: false
    label: Round adjusted counts to integers
    doc: Round adjusted counts to integers
    sd:layout:
      advanced: true
outputs:
  adjusted_feature_bc_matrices_folder:
    type: File
    format: http://edamontology.org/format_3989
    outputSource: estimate_contamination/adjusted_feature_bc_matrices_folder
    label: Compressed folder with adjusted feature-barcode matrices in MEX format
    doc: Compressed folder with adjusted feature-barcode matrices in MEX format
  adjusted_feature_bc_matrices_h5:
    type: File
    format: http://edamontology.org/format_3590
    outputSource: estimate_contamination/adjusted_feature_bc_matrices_h5
    label: Adjusted feature-barcode matrices in HDF5 format
    doc: Adjusted feature-barcode matrices in HDF5 format
  contamination_estimation_plot:
    type: File
    format: http://edamontology.org/format_3508
    outputSource: estimate_contamination/contamination_estimation_plot
    label: Contamination estimation plot
    doc: Contamination estimation plot
  raw_gene_expression_plots:
    type: File?
    format: http://edamontology.org/format_3508
    outputSource: estimate_contamination/raw_gene_expression_plots
    label: Raw gene expression plots
    doc: Raw gene expression plots
  adjusted_gene_expression_plots:
    type: File?
    format: http://edamontology.org/format_3508
    outputSource: estimate_contamination/adjusted_gene_expression_plots
    label: Adjusted gene expression plots
    doc: Adjusted gene expression plots
  raw_gene_expression_to_pure_soup_ratio_plots:
    type: File?
    format: http://edamontology.org/format_3508
    outputSource: estimate_contamination/raw_gene_expression_to_pure_soup_ratio_plots
    label: Raw gene expression to pure soup ratio plots
    doc: Raw gene expression to pure soup ratio plots
  raw_to_adjusted_gene_expression_ratio_plots:
    type: File?
    format: http://edamontology.org/format_3508
    outputSource: estimate_contamination/raw_to_adjusted_gene_expression_ratio_plots
    doc: Raw to adjusted gene expression ratio plots
  soupx_stdout_log:
    type: File
    format: http://edamontology.org/format_2330
    outputSource: estimate_contamination/soupx_stdout_log
    label: SoupX stdout log
    doc: SoupX stdout log
  soupx_stderr_log:
    type: File
    format: http://edamontology.org/format_2330
    outputSource: estimate_contamination/soupx_stderr_log
    label: SoupX stderr log
    doc: SoupX stderr log
steps:
  estimate_contamination:
    run: ../tools/soupx-subworkflow.cwl
    in:
      raw_feature_bc_matrices_folder: raw_feature_bc_matrices_folder
      filtered_feature_bc_matrix_folder: filtered_feature_bc_matrix_folder
      secondary_analysis_report_folder: secondary_analysis_report_folder
      genelist_file: genelist_file
      expression_threshold: expression_threshold
      fdr: fdr
      round_counts: round_counts
    out:
    - adjusted_feature_bc_matrices_folder
    - adjusted_feature_bc_matrices_h5
    - contamination_estimation_plot
    - raw_gene_expression_plots
    - adjusted_gene_expression_plots
    - raw_gene_expression_to_pure_soup_ratio_plots
    - raw_to_adjusted_gene_expression_ratio_plots
    - soupx_stdout_log
    - soupx_stderr_log
label: SoupX Estimate
doc: |
  SoupX Estimate
  ==============
sd:version: 100
