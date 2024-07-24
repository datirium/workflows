cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  rnaseq_sample:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  expression_files:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments"
    doc: "Isoform expression files from RNA-Seq experiments"
    'sd:upstreamSource': "rnaseq_sample/rpkm_isoforms"
    'sd:localLabel': true

  expression_files_genes:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments"
    doc: "Gene expression files from RNA-Seq experiments"
    'sd:upstreamSource': "rnaseq_sample/rpkm_genes"

  expression_files_tss:
    type:
    - "null"
    - File[]
    default: null
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments"
    doc: "Common TSS expression files from RNA-Seq experiments"
    'sd:upstreamSource': "rnaseq_sample/rpkm_common_tss"

  expression_aliases:
    type:
    - "null"
    - string[]
    label: "Isoform expression file aliases"
    doc: "Aliases to make the legend for generated plots. Order corresponds to the isoform expression files"
    'sd:upstreamSource': "rnaseq_sample/alias"

  group_by:
    type:
    - "null"
    - type: enum
      symbols: ["isoforms", "genes", "common tss"]
    default: "isoforms"
    label: "Group by"
    doc: "Grouping method for features: isoforms, genes or common tss"

  genelist_file:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "Gene list to filter input genes. Headerless TSV/CSV file with 1 gene per line"
    doc: "Gene list to filter input genes. Headerless TSV/CSV file with 1 gene per line"
    'sd:layout':
      advanced: true


outputs:

  pca1_vs_pca2_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "PCA1 vs PCA2 plot"
    doc: "PCA1 vs PCA2 plot"
    outputSource: pca/pca1_vs_pca2_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'PCA1 vs PCA2'

  pca1_vs_pca2_plot_pdf:
    type: File
    format: "http://edamontology.org/format_3508"
    label: "PCA1 vs PCA2 plot"
    doc: "PCA1 vs PCA2 plot"
    outputSource: pca/pca1_vs_pca2_plot_pdf

  pca2_vs_pca3_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "PCA2 vs PCA3 plot"
    doc: "PCA2 vs PCA3 plot"
    outputSource: pca/pca2_vs_pca3_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'PCA2 vs PCA3'

  pca2_vs_pca3_plot_pdf:
    type: File
    format: "http://edamontology.org/format_3508"
    label: "PCA2 vs PCA3 plot"
    doc: "PCA2 vs PCA3 plot"
    outputSource: pca/pca2_vs_pca3_plot_pdf

  variance_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "Variance plot"
    doc: "Variance plot"
    outputSource: pca/variance_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Variances'

  variance_plot_pdf:
    type: File
    format: "http://edamontology.org/format_3508"
    label: "Variance plot"
    doc: "Variance plot"
    outputSource: pca/variance_plot_pdf

  pca_3d_html:
    type: File
    format: "http://edamontology.org/format_2331"
    label: "Interactive 3D PCA plot"
    doc: "Plotly generated interactive 3D PCA plot (first three components)"
    outputSource: pca/pca_3d_html

  pca_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "PCA analysis scores results"
    doc: "PCA analysis scores results exported as TSV"
    outputSource: pca/pca_scores_file
    'sd:visualPlugins':
    - scatter3d:
        tab: '3D Plots'
        Title: 'PCA'
        xAxisTitle: 'PCA1'
        yAxisTitle: 'PCA2'
        zAxisTitle: 'PCA3'
        colors: ["#b3de69", "#888888", "#fb8072"]
        height: 600
        data: [$1, $2, $3, $4]

  pca_loadings_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "PCA analysis loadings results"
    doc: "PCA analysis loadings results exported as TSV"
    outputSource: pca/pca_loadings_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'PCA loadings'
        Title: 'PCA loadings'

  pca_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "PCA stdout log"
    doc: "PCA stdout log"
    outputSource: pca/stdout_log

  pca_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "PCA stderr log"
    doc: "PCA stderr log"
    outputSource: pca/stderr_log


steps:

  pca:
    run: ../tools/pca.cwl
    in:
      expression_files:
        source: [group_by, expression_files, expression_files_genes, expression_files_tss]
        valueFrom: |
          ${
              if (self[0] == "isoforms") {
                return self[1];
              } else if (self[0] == "genes") {
                return self[2];
              } else {
                return self[3];
              }
          }
      expression_aliases: expression_aliases
      genelist_file: genelist_file
      target_column:
        default: "Rpkm"
    out:
    - pca1_vs_pca2_plot
    - pca1_vs_pca2_plot_pdf
    - pca2_vs_pca3_plot
    - pca2_vs_pca3_plot_pdf
    - variance_plot
    - variance_plot_pdf
    - pca_3d_html
    - pca_scores_file
    - pca_loadings_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "PCA - Principal Component Analysis"
label: "PCA - Principal Component Analysis"
s:alternateName: ""

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/pca.cwl
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
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/pca.md


doc: |
  Principal Component Analysis
  ---------------

  Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert
  a set of observations of possibly correlated variables (entities each of which takes on various numerical values)
  into a set of values of linearly uncorrelated variables called principal components.

  The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix,
  not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy.