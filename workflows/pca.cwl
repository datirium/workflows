cwlVersion: v1.0
class: Workflow


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
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  expression_file:
    type: File[]
    format: "http://edamontology.org/format_3752"
    label: "Isoform expression files"
    doc: "Isoform expression files"
    'sd:upstreamSource': "rnaseq_sample/rpkm_isoforms"

  legend_name:
    type: string[]
    label: "Isoform expression file aliases"
    doc: "Aliases to make the legend for generated plots. Order corresponds to the isoform expression files"
    'sd:upstreamSource': "rnaseq_sample/alias"


outputs:

  pca1_vs_pca2_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "PCA1 vs PCA2 plot"
    doc: "PCA1 vs PCA2 plot"
    outputSource: pca/pca1_vs_pca2_plot
    'sd:visualPlugins':
    - image:
        Caption: 'PCA1 vs PCA2'

  pca2_vs_pca3_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "PCA2 vs PCA3 plot"
    doc: "PCA2 vs PCA3 plot"
    outputSource: pca/pca2_vs_pca3_plot
    'sd:visualPlugins':
    - image:
        Caption: 'PCA2 vs PCA3'

  variance_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "Variance plot"
    doc: "Variance plot"
    outputSource: pca/variance_plot
    'sd:visualPlugins':
    - image:
        Caption: 'Variances'

  pca_3d_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "PCA1,2,3 plot"
    doc: "First three principal components plot"
    outputSource: pca/pca_3d_plot
    'sd:visualPlugins':
    - image:
        Caption: 'First three principal components'

  pca_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "PCA analysis results"
    doc: "PCA analysis results exported as TSV"
    outputSource: pca/pca_file


steps:

  pca:
    in:
      expression_file: expression_file
      legend_name: legend_name
    out:
    - pca1_vs_pca2_plot
    - pca2_vs_pca3_plot
    - variance_plot
    - pca_3d_plot
    - pca_file
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/pca:v0.0.2
      inputs:
        expression_file:
          type: File[]
          inputBinding:
            prefix: "-i"
        legend_name:
          type: string[]
          inputBinding:
            prefix: "-n"
      outputs:
        pca1_vs_pca2_plot:
          type: File
          outputBinding:
            glob: "*001.png"
        pca2_vs_pca3_plot:
          type: File
          outputBinding:
            glob: "*002.png"
        variance_plot:
          type: File
          outputBinding:
            glob: "*003.png"
        pca_3d_plot:
          type: File
          outputBinding:
            glob: "*004.png"
        pca_file:
          type: File
          outputBinding:
            glob: "*.tsv"
      baseCommand: ["run_pca.R"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

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

doc: |
  Principal Component Analysis
  =====================================

  Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert
  a set of observations of possibly correlated variables (entities each of which takes on various numerical values)
  into a set of values of linearly uncorrelated variables called principal components.

  The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix,
  not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy.
