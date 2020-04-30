cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


'sd:metadata':
  - "../metadata/indices-header.cwl"


inputs:

  genome:
    type: string
    label: "Genome type"
    doc: "Genome type, such as mm10, hg19, hg38, etc"

  fasta_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes"


outputs:

  indices_folder:
    type: Directory
    label: "Bowtie indices"
    doc: "Bowtie generated indices folder"
    outputSource: bowtie_build/indices_folder

  stdout_log:
    type: File
    label: "Bowtie stdout log"
    doc: "Bowtie generated stdout log"
    outputSource: bowtie_build/stdout_log

  stderr_log:
    type: File
    label: "Bowtie stderr log"
    doc: "Bowtie generated stderr log"
    outputSource: bowtie_build/stderr_log


steps:

  bowtie_build:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta_file
      index_base_name: genome
    out:
    - indices_folder
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "Build Bowtie indices"
label: "Build Bowtie indices"
s:alternateName: "Build Bowtie indices"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/dummy/bowtie-index.cwl
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
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681


doc:
  $include: ../descriptions/bowtie-index.md


# doc: |
#   Workflow runs [Bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016) to build indices for reference
#   genome provided in a single FASTA file as fasta_file input. Generated indices are saved in a folder with the name that
#   corresponds to the input genome