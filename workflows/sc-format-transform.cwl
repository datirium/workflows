cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  compressed_sparse_matrix:
    type: File
    label: "TAR-gzipped folder with the feature-barcode matrix in MEX format"
    doc: |
      Compressed folder with the feature-barcode
      matrix from the Cell Ranger Count/Aggregate
      experiment in MEX format (TAR-gzipped).

  metadata:
    type: File?
    label: "Aggregation metadata in TSV format"
    doc: |
      Aggregation metadata file from the
      Cell Ranger Aggregate experiment


outputs:

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: pipe/filtered_feature_bc_matrix_folder
    label: "TAR-gzipped folder with the feature-barcode matrix in MEX format"
    doc: |
      Compressed folder with the feature-barcode
      matrix from the Cell Ranger Count/Aggregate
      experiment in MEX format (TAR-gzipped).

  aggregation_metadata:
    type: File?
    outputSource: pipe/aggregation_metadata
    label: "Aggregation metadata in TSV format"
    doc: |
      Aggregation metadata file from the
      Cell Ranger Aggregate experiment


steps:

  pipe:
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
            RNDM_PREFIX=$(tr -dc a-z </dev/urandom | head -c 5)
            cp $0 ./${RNDM_PREFIX}_bc_matrix.tar.gz
            if [[ -n $1 ]]; then
              cat $1 | tr "," "\t" > ./${RNDM_PREFIX}_aggr.tsv
            fi
          inputBinding:
            position: 1
        compressed_sparse_matrix:
          type: File
          inputBinding:
            position: 2
        metadata:
          type: File?
          inputBinding:
            position: 3
      outputs:
        filtered_feature_bc_matrix_folder:
          type: File
          outputBinding:
            glob: "*_bc_matrix.tar.gz"
        aggregation_metadata:
          type: File?
          outputBinding:
            glob: "*_aggr.tsv"
      baseCommand: [bash, '-c']
    in:
      compressed_sparse_matrix: compressed_sparse_matrix
      metadata: metadata
    out:
    - filtered_feature_bc_matrix_folder
    - aggregation_metadata


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Format Transform"
s:name: "Single-cell Format Transform"
s:alternateName: "Transforms single-cell sequencing data formats into Cell Ranger like output"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/sc-format-transform.cwl
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
  Single-cell Format Transform
  
  Transforms single-cell sequencing data formats
  into Cell Ranger like output.
