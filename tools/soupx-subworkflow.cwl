cwlVersion: v1.0
class: Workflow


requirements:
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


inputs:

  raw_feature_bc_matrices_folder:
    type: File
    doc: "Compressed folder with unfiltered feature-barcode matrices"

  filtered_feature_bc_matrix_folder:
    type: File
    doc: "Compressed folder with filtered feature-barcode matrices"

  secondary_analysis_report_folder:
    type: File
    doc: "Compressed folder with secondary analysis results"

  genelist_file:
    type: File?
    doc: "Target genes list. Headerless text file with 1 gene per line"

  expression_threshold:
    type: float?
    doc: "Expression threshold for displaying target genes on a plot (expression > threshold)"

  fdr:
    type: float?
    doc: "FDR cutoff for expression ratio plots"

  round_counts:
    type: boolean?
    doc: "Round adjusted counts to integers"

  matrix_format_version:
    type:
    - "null"
    - type: enum
      name: "matrix_format_version"
      symbols: ["2","3"]
    doc: "Output matrix format version. Corresponds to the latest Cell Ranger matrix format"

  output_prefix:
    type: string?
    doc: "Output prefix"


outputs:

  adjusted_feature_bc_matrices_folder:
    type: File
    outputSource: compress_adjusted_feature_bc_matrices_folder/compressed_folder
    doc: "Compressed folder with adjusted feature-barcode matrices in MEX format"

  adjusted_feature_bc_matrices_h5:
    type: File
    outputSource: estimate_contamination/adjusted_feature_bc_matrices_h5
    doc: "Adjusted feature-barcode matrices in HDF5 format"

  contamination_estimation_plot:
    type: File
    outputSource: estimate_contamination/contamination_estimation_plot
    doc: "Contamination estimation plot"

  raw_gene_expression_plots:
    type: File?
    outputSource: estimate_contamination/raw_gene_expression_plots
    doc: "Raw gene expression plots"

  adjusted_gene_expression_plots:
    type: File?
    outputSource: estimate_contamination/adjusted_gene_expression_plots
    doc: "Adjusted gene expression plots"

  raw_gene_expression_to_pure_soup_ratio_plots:
    type: File?
    outputSource: estimate_contamination/raw_gene_expression_to_pure_soup_ratio_plots
    doc: "Raw gene expression to pure soup ratio plots"
  
  raw_to_adjusted_gene_expression_ratio_plots:
    type: File?
    outputSource: estimate_contamination/raw_to_adjusted_gene_expression_ratio_plots
    doc: "Raw to adjusted gene expression ratio plots"

  soupx_stdout_log:
    type: File
    outputSource: estimate_contamination/stdout_log
    doc: "SoupX stdout log"

  soupx_stderr_log:
    type: File
    outputSource: estimate_contamination/stderr_log
    doc: "SoupX stderr log"


steps:

  extract_count_matrices_to_folder:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: scidap/scidap:v0.0.4
      inputs:
        script:
          type: string?
          default: |
            set -- "$0" "$@"
            mkdir counts && cd counts
            for i in "$@";
                do tar xvf $i;
            done;
          inputBinding:
            position: 1
        compressed_files:
          type: File[]
          inputBinding:
            position: 2
      outputs:
        output_folder:
          type: Directory
          outputBinding:
            glob: "counts"
      baseCommand: [bash, '-c']
    in:
      compressed_files:
        source:
        - raw_feature_bc_matrices_folder
        - filtered_feature_bc_matrix_folder
        - secondary_analysis_report_folder
    out:
    - output_folder

  estimate_contamination:
    run: soupx.cwl
    in:
      feature_bc_matrices_folder: extract_count_matrices_to_folder/output_folder
      genelist_file: genelist_file
      expression_threshold: expression_threshold
      fdr: fdr
      round_counts: round_counts
      matrix_format_version: matrix_format_version
      output_prefix: output_prefix
    out:
    - adjusted_feature_bc_matrices_folder
    - adjusted_feature_bc_matrices_h5
    - contamination_estimation_plot
    - raw_gene_expression_plots
    - adjusted_gene_expression_plots
    - raw_gene_expression_to_pure_soup_ratio_plots
    - raw_to_adjusted_gene_expression_ratio_plots
    - stdout_log
    - stderr_log

  compress_adjusted_feature_bc_matrices_folder:
    run: tar-compress.cwl
    in:
      folder_to_compress: estimate_contamination/adjusted_feature_bc_matrices_folder
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "SoupX (workflow) - an R package for the estimation and removal of cell free mRNA contamination"
label: "SoupX (workflow) - an R package for the estimation and removal of cell free mRNA contamination"
s:alternateName: "SoupX (workflow) - an R package for the estimation and removal of cell free mRNA contamination"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/soupx-subworkflow.cwl
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
  Wrapped in a workflow SoupX tool for easy access to Cell Ranger pipeline compressed outputs.


s:about: |
  Wrapped in a workflow SoupX tool for easy access to Cell Ranger pipeline compressed outputs.