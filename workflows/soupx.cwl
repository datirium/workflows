cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  sc_rnaseq_sample:
    - "single-cell-preprocess-cellranger.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  raw_feature_bc_matrices_folder:
    type: File
    format: "http://edamontology.org/format_3989"
    label: "scRNA-Seq Cell Ranger Experiment"
    doc: "Compressed folder with unfiltered feature-barcode matrices"
    'sd:upstreamSource': "sc_rnaseq_sample/raw_feature_bc_matrices_folder"
    'sd:localLabel': true

  filtered_feature_bc_matrix_folder:
    type: File
    format: "http://edamontology.org/format_3989"
    label: "scRNA-Seq Cell Ranger Experiment"
    doc: "Compressed folder with filtered feature-barcode matrices"
    'sd:upstreamSource': "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

  secondary_analysis_report_folder:
    type: File
    format: "http://edamontology.org/format_3989"
    label: "scRNA-Seq Cell Ranger Experiment"
    doc: "Compressed folder with secondary analysis results"
    'sd:upstreamSource': "sc_rnaseq_sample/secondary_analysis_report_folder"
    'sd:localLabel': true

  genelist_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Target genes list. Headerless text file with 1 gene per line"
    doc: "Target genes list. Headerless text file with 1 gene per line"

  expression_threshold:
    type: float?
    default: 0
    label: "Expression threshold for displaying target genes on a plot (expression > threshold)"
    doc: "Expression threshold for displaying target genes on a plot (expression > threshold)"
    'sd:layout':
      advanced: true

  fdr:
    type: float?
    default: 0.05
    label: "FDR cutoff for expression ratio plots"
    doc: "FDR cutoff for expression ratio plots"
    'sd:layout':
      advanced: true

  round_counts:
    type: boolean?
    default: false
    label: "Round adjusted counts to integers"
    doc: "Round adjusted counts to integers"
    'sd:layout':
      advanced: true


outputs:

  adjusted_feature_bc_matrices_folder:
    type: File
    format: "http://edamontology.org/format_3989"
    outputSource: estimate_contamination/adjusted_feature_bc_matrices_folder
    label: "Compressed folder with adjusted feature-barcode matrices in MEX format"
    doc: "Compressed folder with adjusted feature-barcode matrices in MEX format"

  adjusted_feature_bc_matrices_h5:
    type: File
    format: "http://edamontology.org/format_3590"
    outputSource: estimate_contamination/adjusted_feature_bc_matrices_h5
    label: "Adjusted feature-barcode matrices in HDF5 format"
    doc: "Adjusted feature-barcode matrices in HDF5 format"

  contamination_estimation_plot:
    type: File
    format: "http://edamontology.org/format_3508"
    outputSource: estimate_contamination/contamination_estimation_plot
    label: "Contamination estimation plot"
    doc: "Contamination estimation plot"

  raw_gene_expression_plots:
    type: File?
    format: "http://edamontology.org/format_3508"
    outputSource: estimate_contamination/raw_gene_expression_plots
    label: "Raw gene expression plots"
    doc: "Raw gene expression plots"

  adjusted_gene_expression_plots:
    type: File?
    format: "http://edamontology.org/format_3508"
    outputSource: estimate_contamination/adjusted_gene_expression_plots
    label: "Adjusted gene expression plots"
    doc: "Adjusted gene expression plots"

  raw_gene_expression_to_pure_soup_ratio_plots:
    type: File?
    format: "http://edamontology.org/format_3508"
    outputSource: estimate_contamination/raw_gene_expression_to_pure_soup_ratio_plots
    label: "Raw gene expression to pure soup ratio plots"
    doc: "Raw gene expression to pure soup ratio plots"
  
  raw_to_adjusted_gene_expression_ratio_plots:
    type: File?
    format: "http://edamontology.org/format_3508"
    outputSource: estimate_contamination/raw_to_adjusted_gene_expression_ratio_plots
    doc: "Raw to adjusted gene expression ratio plots"

  soupx_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: estimate_contamination/soupx_stdout_log
    label: "SoupX stdout log"
    doc: "SoupX stdout log"

  soupx_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: estimate_contamination/soupx_stderr_log
    label: "SoupX stderr log"
    doc: "SoupX stderr log"


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


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/9.0/schemaorg-current-http.rdf

s:name: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"
label: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"
s:alternateName: "SoupX - an R package for the estimation and removal of cell free mRNA contamination"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/soupx.cwl
s:codeRepository: https://github.com/datirium/workflows
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
  Devel version of Single-Cell Advanced Cell Ranger Pipeline (SoupX)
  =================================================================

